// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// qd_linear_solver.hpp — the direct solve K u = b of a quasi-dynamic stage
// behind one interface (GPU plan G3).  K does not change during a run, so a
// back end analyses and factors it once, in MakeQDLinearSolver, and every
// stage is one Mult.  Back ends:
//
//   MUMPS, MUMPS_BLR   CPU (MUMPS through MFEM), the production solver.
//   CUDSS              GPU (MFEM's CuDSSSolver, MFEM_USE_CUDSS).  Under MPI
//                      it gives cuDSS each rank's rows of K (row
//                      distribution, cuDSS's multi-GPU mode, also on one
//                      rank) and loads the communication layer from
//                      CUDSS_COMM_LIB or the path MFEM was built with.
//                      Mult reads and writes device vectors.
//   PCG                MFEM CGSolver (UCERF 3-D plan rev5, Q1.1).  No
//                      factorization: setup is the preconditioner, and each
//                      Mult iterates to QDPCGConfig::rel_tol.  Two
//                      preconditioners (Q1 probe, 2026-11-24, p = 1, 1e-10):
//                      "amg", scalar BoomerAMG on K with strength 0.8 (T-plane
//                      30-36 iterations; UCERF h_f 6 km, 3.5 M DOFs, 32), and
//                      "aux_h1", the two-level auxiliary space on the
//                      continuous p = 1 space (QDAuxH1Preconditioner: 27 and
//                      29 iterations, 1.4x faster per solve).  BoomerAMG in
//                      elasticity mode (SetElasticityOptions) on the DG K
//                      needs 450-580 and is not offered.
//
// All are given K as symmetric positive definite, the full matrix.

#ifndef MFEM_SEAS_QD_LINEAR_SOLVER_HPP
#define MFEM_SEAS_QD_LINEAR_SOLVER_HPP

#include "mfem.hpp"

#ifdef MFEM_USE_MPI

#include <iomanip>
#include <memory>
#include <sstream>
#include <string>

namespace mfem
{
namespace seas
{

enum class QDLinearSolverKind { MUMPS, MUMPS_BLR, CUDSS, PCG };

/// Settings of the PCG back end ([cycle.qd] solver = "pcg").
struct QDPCGConfig
{
   double rel_tol = 1e-10;     ///< stop at ||r|| <= rel_tol ||b|| (preconditioned norm)
   int max_it = 1000;          ///< a solve that has not converged by then aborts
   /// "amg": BoomerAMG on K itself.  "aux_h1": two-level auxiliary space,
   /// a smoother on K and one AMG V-cycle on the continuous (H1) p = 1
   /// space of the same mesh (QDAuxH1Preconditioner).
   std::string preconditioner = "amg";
   double amg_strength = 0.8;  ///< BoomerAMG strength threshold ("amg" only)
   /// false: l1 Gauss-Seidel ("amg": MFEM's default AMG smoother; "aux_h1":
   /// one l1-GS sweep on K).  true: l1 Jacobi, the GPU-friendly choice
   /// ("amg": hypre relax type 18; "aux_h1": two l1-Jacobi sweeps on K and
   /// relax type 18 in the coarse AMG).
   bool l1_jacobi = false;
};

/// Two-level auxiliary-space preconditioner for the SIPG matrix K of a DG
/// vector space (UCERF 3-D plan rev5, Q1.1; Q1 probe 2026-11-24):
///   x  = S b                         (smoother on K, from zero)
///   x += P C P^T (b - K x)           (coarse correction)
///   x += S (b - K x)                 (smoother again)
/// P interpolates the continuous H1 p = 1 vector space of the same mesh
/// into the DG space (element by element, at the DG nodes), Kc = P^T K P
/// (Galerkin) and C is one BoomerAMG V-cycle on Kc in systems mode.  S is
/// symmetric (l1 GS or l1 Jacobi), so the cycle is a symmetric
/// preconditioner for CG.  Kc has the size of the continuous p = 1 space:
/// about 1/21 of the DG p = 1 space on tetrahedra.
class QDAuxH1Preconditioner : public Solver
{
public:
   QDAuxH1Preconditioner(const HypreParMatrix &K, ParFiniteElementSpace &dg,
                         bool l1_jacobi)
      : Solver(K.Height(), K.Width()), K_(K),
        h1_fec_(1, dg.GetParMesh()->Dimension())
   {
      ParMesh &pmesh = *dg.GetParMesh();
      const int vdim = dg.GetVDim();
      MFEM_VERIFY(dg.IsDGSpace(), "QD PCG aux_h1: the space is not DG");
      MFEM_VERIFY(dg.GetTrueVSize() == dg.GetVSize() && K.Height() == dg.GetVSize(),
                  "QD PCG aux_h1: K does not match the DG space");
      h1_ = std::make_unique<ParFiniteElementSpace>(&pmesh, &h1_fec_, vdim,
                                                    Ordering::byNODES);
      SparseMatrix Ploc(dg.GetVSize(), h1_->GetVSize());
      Array<int> vd, vh;
      Vector sh;
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
         const FiniteElement &ed = *dg.GetFE(e), &eh = *h1_->GetFE(e);
         const IntegrationRule &nodes = ed.GetNodes();
         dg.GetElementVDofs(e, vd);
         h1_->GetElementVDofs(e, vh);
         const int nd = ed.GetDof(), nh = eh.GetDof();
         sh.SetSize(nh);
         for (int i = 0; i < nd; i++)
         {
            eh.CalcShape(nodes.IntPoint(i), sh);
            for (int j = 0; j < nh; j++)
            {
               if (std::abs(sh(j)) < 1e-15) { continue; }
               for (int c = 0; c < vdim; c++)
               {
                  Ploc.Set(vd[c * nd + i], vh[c * nh + j], sh(j));
               }
            }
         }
      }
      Ploc.Finalize();
      HypreParMatrix Pl(K.GetComm(), K.GetGlobalNumRows(), h1_->GlobalVSize(),
                        K.GetRowStarts(), h1_->GetDofOffsets(), &Ploc);
      P_.reset(ParMult(&Pl, h1_->Dof_TrueDof_Matrix()));
      Kc_.reset(RAP(&K, P_.get()));
      S_.SetType(l1_jacobi ? HypreSmoother::l1Jacobi : HypreSmoother::l1GS,
                 l1_jacobi ? 2 : 1);
      S_.SetOperator(K);
      amg_ = std::make_unique<HypreBoomerAMG>(*Kc_);
      amg_->SetPrintLevel(0);
      amg_->SetSystemsOptions(vdim, true);
      if (l1_jacobi) { amg_->SetRelaxType(18); }
      Vector z(Kc_->Height()), zz(Kc_->Height());
      z = 0.0;
      zz = 0.0;
      amg_->Setup(z, zz);
   }

   /// CGSolver::SetOperator hands its operator on to the preconditioner:
   /// only the K this was built for is accepted.
   void SetOperator(const Operator &op) override
   {
      MFEM_VERIFY(&op == &K_, "QDAuxH1Preconditioner: K is fixed at construction");
   }

   void Mult(const Vector &b, Vector &x) const override
   {
      r_.SetSize(b.Size());
      d_.SetSize(b.Size());
      rc_.SetSize(P_->Width());
      ec_.SetSize(P_->Width());
      x = 0.0;
      S_.Mult(b, x);
      K_.Mult(x, r_);
      r_.Neg();
      r_ += b;
      P_->MultTranspose(r_, rc_);
      ec_ = 0.0;
      amg_->Mult(rc_, ec_);
      P_->Mult(ec_, d_);
      x += d_;
      K_.Mult(x, r_);
      r_.Neg();
      r_ += b;
      d_ = 0.0;
      S_.Mult(r_, d_);
      x += d_;
   }

   /// Global size of the coarse (continuous p = 1) space.
   HYPRE_BigInt CoarseSize() const { return Kc_->GetGlobalNumRows(); }

private:
   const HypreParMatrix &K_;
   H1_FECollection h1_fec_;
   std::unique_ptr<ParFiniteElementSpace> h1_;
   std::unique_ptr<HypreParMatrix> P_, Kc_;
   mutable HypreSmoother S_;
   std::unique_ptr<HypreBoomerAMG> amg_;
   mutable Vector r_, d_, rc_, ec_;
};

/// PCG on a fixed K, preconditioned by BoomerAMG ("amg") or by
/// QDAuxH1Preconditioner ("aux_h1", needs the DG space).  Mult starts from
/// x = 0 and aborts when the solve does not converge within max_it: a stage
/// with an unconverged displacement would hand the friction solve a wrong
/// traction.
class QDPCGSolver : public Solver
{
public:
   QDPCGSolver(const HypreParMatrix &K, const QDPCGConfig &config,
               ParFiniteElementSpace *dg = nullptr)
      : Solver(K.Height(), K.Width()), config_(config), cg_(K.GetComm())
   {
      MFEM_VERIFY(config.rel_tol > 0.0 && config.rel_tol < 1.0,
                  "QD PCG: rel_tol must lie in (0, 1), got " << config.rel_tol);
      MFEM_VERIFY(config.max_it >= 1,
                  "QD PCG: max_it must be at least 1, got " << config.max_it);
      MFEM_VERIFY(config.amg_strength > 0.0 && config.amg_strength < 1.0,
                  "QD PCG: amg_strength must lie in (0, 1), got "
                  << config.amg_strength);
      MFEM_VERIFY(config.preconditioner == "amg" || config.preconditioner == "aux_h1",
                  "QD PCG: preconditioner must be \"amg\" or \"aux_h1\", got \""
                  << config.preconditioner << "\"");
      if (config.preconditioner == "aux_h1")
      {
         MFEM_VERIFY(dg, "QD PCG: preconditioner \"aux_h1\" needs the DG space");
         auto *aux = new QDAuxH1Preconditioner(K, *dg, config.l1_jacobi);
         coarse_size_ = aux->CoarseSize();
         prec_.reset(aux);
      }
      else
      {
         auto *amg = new HypreBoomerAMG(K);
         amg->SetPrintLevel(0);
         amg->SetStrengthThresh(config.amg_strength);
         if (config.l1_jacobi) { amg->SetRelaxType(18); }
         // Build the hierarchy now, so the setup time is the setup's.
         Vector b(K.Height()), x(K.Height());
         b = 0.0;
         x = 0.0;
         amg->Setup(b, x);
         prec_.reset(amg);
      }
      cg_.SetRelTol(config.rel_tol);
      cg_.SetAbsTol(0.0);
      cg_.SetMaxIter(config.max_it);
      cg_.SetPrintLevel(0);
      cg_.SetPreconditioner(*prec_);
      cg_.SetOperator(K);
   }

   void SetOperator(const Operator &) override
   {
      MFEM_ABORT("QDPCGSolver: K is fixed at construction");
   }

   void Mult(const Vector &b, Vector &x) const override
   {
      cg_.Mult(b, x);
      last_iterations_ = cg_.GetNumIterations();
      last_rel_norm_ = cg_.GetFinalRelNorm();
      total_iterations_ += last_iterations_;
      num_solves_++;
      MFEM_VERIFY(cg_.GetConverged(),
                  "QD PCG: no convergence in " << config_.max_it
                  << " iterations (relative residual " << last_rel_norm_
                  << ", rel_tol " << config_.rel_tol << ", solve "
                  << num_solves_ << ")");
   }

   const QDPCGConfig &Config() const { return config_; }
   int LastIterations() const { return last_iterations_; }
   double LastRelNorm() const { return last_rel_norm_; }
   long long TotalIterations() const { return total_iterations_; }
   long long NumSolves() const { return num_solves_; }
   /// "aux_h1": global size of the continuous coarse space (0 for "amg").
   HYPRE_BigInt CoarseSize() const { return coarse_size_; }

private:
   QDPCGConfig config_;
   mutable CGSolver cg_;
   std::unique_ptr<Solver> prec_;
   HYPRE_BigInt coarse_size_ = 0;
   mutable int last_iterations_ = 0;
   mutable double last_rel_norm_ = 0.0;
   mutable long long total_iterations_ = 0;
   mutable long long num_solves_ = 0;
};

/// What setting up a back end cost, for the log.
struct QDLinearSolverSetup
{
   std::string name;
   double seconds = 0.0;      ///< analysis + factorization, slowest rank
   double device_bytes = -1;  ///< device memory taken, largest rank (-1: CPU)
};

/// True when this MFEM build has the back end @a kind.
inline bool QDLinearSolverAvailable(QDLinearSolverKind kind)
{
   switch (kind)
   {
#ifdef MFEM_USE_MUMPS
      case QDLinearSolverKind::MUMPS:
      case QDLinearSolverKind::MUMPS_BLR: return true;
#endif
#ifdef MFEM_USE_CUDSS
      case QDLinearSolverKind::CUDSS: return true;
#endif
      case QDLinearSolverKind::PCG: return true;
      default: return false;
   }
}

namespace detail
{
/// Free device memory in bytes (0 without a GPU backend).
inline double FreeDeviceBytes()
{
#if defined(MFEM_USE_CUDA)
   if (Device::Allows(Backend::CUDA_MASK))
   {
      size_t free_b = 0, total_b = 0;
      MFEM_GPU_CHECK(cudaMemGetInfo(&free_b, &total_b));
      return static_cast<double>(free_b);
   }
#endif
   return 0.0;
}
} // namespace detail

/// Build the back end @a kind for @a K: analysis and factorization happen
/// here.  @a mumps_print_level and @a blr_tol apply to the MUMPS back ends,
/// @a pcg (required) to PCG, where the setup is the preconditioner; @a dg
/// (the DG space of K) is needed by PCG's "aux_h1" preconditioner.
/// Collective over the communicator of @a K.
inline std::unique_ptr<Solver> MakeQDLinearSolver(
   QDLinearSolverKind kind, const HypreParMatrix &K, int mumps_print_level,
   double blr_tol, QDLinearSolverSetup *setup = nullptr,
   const QDPCGConfig *pcg = nullptr, ParFiniteElementSpace *dg = nullptr)
{
   MFEM_VERIFY(QDLinearSolverAvailable(kind),
               "QD linear solver: this MFEM build has no "
               << (kind == QDLinearSolverKind::CUDSS ? "cuDSS (MFEM_USE_CUDSS)"
                                                     : "MUMPS (MFEM_USE_MUMPS)"));
   MFEM_VERIFY(kind != QDLinearSolverKind::PCG || pcg != nullptr,
               "QD linear solver: PCG needs its QDPCGConfig");
   MPI_Comm comm = K.GetComm();
   MPI_Barrier(comm);
   const double free0 = detail::FreeDeviceBytes();
   const double t0 = MPI_Wtime();
   std::unique_ptr<Solver> solver;
   std::string name;
   switch (kind)
   {
#ifdef MFEM_USE_MUMPS
      case QDLinearSolverKind::MUMPS:
      case QDLinearSolverKind::MUMPS_BLR:
      {
         auto *mumps = new MUMPSSolver(comm);
         mumps->SetMatrixSymType(MUMPSSolver::MatType::SYMMETRIC_POSITIVE_DEFINITE);
         mumps->SetPrintLevel(mumps_print_level);
         if (kind == QDLinearSolverKind::MUMPS_BLR) { mumps->SetBLRTol(blr_tol); }
         mumps->SetOperator(K);
         solver.reset(mumps);
         name = (kind == QDLinearSolverKind::MUMPS) ? "MUMPS" : "MUMPS BLR";
         break;
      }
#endif
#ifdef MFEM_USE_CUDSS
      case QDLinearSolverKind::CUDSS:
      {
         auto *cudss = new CuDSSSolver(comm);
         cudss->SetMatrixSymType(CuDSSSolver::SYMMETRIC_POSITIVE_DEFINITE);
         cudss->SetMatrixViewType(CuDSSSolver::FULL);
         cudss->SetOperator(K);
         solver.reset(cudss);
         int nranks;
         MPI_Comm_size(comm, &nranks);
         name = "cuDSS (" + std::to_string(nranks) + " rank"
                + (nranks > 1 ? "s" : "") + ", row distribution)";
         break;
      }
#endif
      case QDLinearSolverKind::PCG:
      {
         auto *p = new QDPCGSolver(K, *pcg, dg);
         solver.reset(p);
         std::ostringstream o;
         if (pcg->preconditioner == "aux_h1")
         {
            o << "PCG + aux_h1 (" << (pcg->l1_jacobi ? "2 l1-Jacobi sweeps"
                                                     : "1 l1-GS sweep")
              << ", BoomerAMG systems on the continuous p = 1 space, "
              << p->CoarseSize() << " DOFs";
         }
         else
         {
            o << "PCG + BoomerAMG (strength " << pcg->amg_strength << ", "
              << (pcg->l1_jacobi ? "l1-Jacobi" : "default") << " smoother";
         }
         o << ", rel_tol " << pcg->rel_tol << ", max_it " << pcg->max_it << ")";
         name = o.str();
         break;
      }
      default: break;
   }
   if (Device::Allows(Backend::CUDA_MASK | Backend::HIP_MASK))
   {
      MFEM_DEVICE_SYNC;
   }
   const double seconds = MPI_Wtime() - t0;
   const double used = (kind == QDLinearSolverKind::CUDSS)
                       ? free0 - detail::FreeDeviceBytes() : -1.0;
   if (setup)
   {
      setup->name = name;
      MPI_Allreduce(&seconds, &setup->seconds, 1, MPI_DOUBLE, MPI_MAX, comm);
      MPI_Allreduce(&used, &setup->device_bytes, 1, MPI_DOUBLE, MPI_MAX, comm);
   }
   return solver;
}

/// One line for the log, e.g. "[QD-SOLVER] cuDSS (1 rank, row distribution):
/// analysis + factorization 1.234 s (slowest rank), device memory 2.10 GiB
/// (largest rank, ...)".
inline std::string QDLinearSolverSummary(const QDLinearSolverSetup &s)
{
   std::ostringstream o;
   const bool iterative = s.name.rfind("PCG", 0) == 0;
   o << "[QD-SOLVER] " << s.name
     << (iterative ? ": AMG setup " : ": analysis + factorization ")
     << std::fixed << std::setprecision(3) << s.seconds << " s (slowest rank)";
   if (s.device_bytes >= 0.0)
   {
      o << ", device memory " << std::setprecision(2)
        << s.device_bytes / (1024.0 * 1024.0 * 1024.0)
        << " GiB (largest rank, drop in free memory over the setup)";
   }
   return o.str();
}

} // namespace seas
} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_SEAS_QD_LINEAR_SOLVER_HPP
