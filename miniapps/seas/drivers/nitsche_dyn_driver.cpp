// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_dyn_driver.cpp — Nitsche dynamic rupture, plan item N-07 (build
// sequence rev 5.2, phase P7): the driver.
//
// Reads the SAFS spatial deck through spatial::LoadSpatialFrictionConfig (the
// SAFS driver's own parser) plus the Nitsche-only keys the SAFS parser does
// not know ([numerics] penalty_factor, the fault-coordinate specification's
// [fault_geometry] keys, [diagnostics] energy_balance), builds the mesh, the
// material, the boundary configuration and the fault setup, constructs the
// NitscheWaveOperator (N-01) and advances it with MFEM's
// CentralDifferenceSolver (N-02).  Nothing benchmark-specific lives here: the
// law comes from a factory keyed on [meta] law, the friction parameters from
// the SAFS resolver, the prestress and the nucleation from the deck as global
// tensors (spec I1), the station list from the [problem] tag.
//
// Outputs (in --output-dir):
//   <prefix>_station_<name>.dat   TPV102StationWriter (column layout of the
//                                 existing runs) through a DOFData adapter;
//   <prefix>_station_conventions.txt  frame, pin and sign map (I3);
//   <prefix>_station_<name>_copies.dat  per station, the min/max over the DG
//                                 copies at the station point (--station-copies
//                                 average, the default);
//   <prefix>_fault_r<rank>.csv    fault-field dump every --dump-every steps
//                                 (with Prakash-Clifton one more column, sigma*);
//   <prefix>_rupture_time_r<rank>.csv  first |V| > 1 mm/s per fault DOF;
//   stdout log every --log-every steps.
//
// Station copies (Chunhui Zhao, 2026-10-14): a station sits on a fault point
// that carries one DOF per adjacent fault face (DG copies).  With
// --station-copies average (default) the station files hold the face-area-
// weighted average over those copies (summed on rank 0 in face-centroid
// order); --station-copies first keeps the writer's single nearest copy.
//
// Station columns (SCEC map, plan "The SCEC conventions, mapped"):
//   slip1 = delta.dhat (v-slip), slip2 = delta.shat (h-slip),
//   V1 = V.dhat, V2 = V.shat,
//   tau1 = -T.dhat, tau2 = -T.shat (SCEC shear tractions, Pa),
//   sigma_n = -T.nhat (compression positive, Pa),
// with T the total traction (Nitsche traction at the end of the step plus
// prestress plus nucleation) on the face with normal nhat.
//
// Optional [friction.prakash_clifton] table (Chunhui Zhao, 2026-10-14): v_star
// (m/s) and length (m), both required and > 0; its presence switches on the
// operator's Prakash-Clifton regularization of the normal stress the law
// sees, seeded with the initial sigma_n per DOF.  Absent = every bit as
// before.  The station files are unchanged (sigma_n there stays sigma_eff);
// sigma* goes to the fault CSV dump and to fault.vtkhdf (normal_stress_pc).
//
// Scope (Chunhui Zhao, 2026-10-14): tetrahedral meshes; constant material
// (the sidecar / depth-profile / halfspace selections are P8);
// rate-and-state (slip weakening comes with TPV205).

#include "mfem.hpp"
#include "../domain/nitsche_wave_operator.hpp"
#include "../domain/nitsche_energy.hpp"
#include "../solver/nitsche_cfl.hpp"
#include "../friction/rate_state_closure.hpp"
#include "../fault/fault_frame.hpp"
#include "../spatial/code/spatial_friction.hpp"
#include "../dynamic/spatial_nucleation.hpp"
#include "../dynamic/tpv102_stations.hpp"
#ifdef MFEM_USE_HDF5
#include "../io/fault_vtkhdf_writer.hpp"
#include "../io/paraview_output.hpp"   // ParaViewOutput::ProbeH5ZZfpPluginOrAbort
#endif

#include <toml.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
namespace ff = mfem::seas::fault_frame;

namespace
{

/// P3's measured stable constants C(p), p = 1..4 (dt = cfl C(p) h/(p^2 c_p)).
const real_t kCeig[5] = {0.0, 0.113, 0.212, 0.282, 0.332};

int g_rank = 0;

std::string ArgStr(int argc, char *argv[], const char *flag, const std::string &def)
{
   for (int i = 1; i + 1 < argc; i++) { if (std::string(argv[i]) == flag) { return argv[i + 1]; } }
   return def;
}
real_t ArgReal(int argc, char *argv[], const char *flag, real_t def)
{
   const std::string v = ArgStr(argc, argv, flag, "");
   return v.empty() ? def : std::atof(v.c_str());
}
int ArgInt(int argc, char *argv[], const char *flag, int def)
{
   const std::string v = ArgStr(argc, argv, flag, "");
   return v.empty() ? def : std::atoi(v.c_str());
}

const toml::value *DeckTable(const toml::value &root, const char *name)
{
   return root.contains(name) ? &root.at(name) : nullptr;
}
real_t TomlReal(const toml::value *t, const char *key, real_t def)
{
   if (!t || !t->contains(key)) { return def; }
   const auto &v = t->at(key);
   if (v.is_floating()) { return static_cast<real_t>(v.as_floating()); }
   if (v.is_integer()) { return static_cast<real_t>(v.as_integer()); }
   MFEM_ABORT("deck key '" << key << "' must be a number");
   return def;
}
std::string TomlStr(const toml::value *t, const char *key, const std::string &def)
{
   if (!t || !t->contains(key)) { return def; }
   return t->at(key).as_string();
}
bool TomlBool(const toml::value *t, const char *key, bool def)
{
   if (!t || !t->contains(key)) { return def; }
   return t->at(key).as_boolean();
}

/// Station copy average (header, "Station copies").  Every owned DOF within
/// kTie of the globally nearest station distance (the writer's own metric)
/// gets the face-area-weighted average over all DOFs at the station point P
/// (P = the lexicographically smallest such position, copies within kTie of
/// P); rank 0 sums in face-centroid order, so the sum does not depend on the
/// partition, and writes the envelope (min/max over copies).  Collective.
class StationCopyAverage
{
public:
   static constexpr real_t kTie = 1e-6;   // m
   StationCopyAverage(NitscheWaveOperator &op, const std::vector<TPV102Station> &st,
                      const std::string &dir, const std::string &prefix)
      : ns_(static_cast<int>(st.size()))
   {
      ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
      const int nd = E.GetNumFaultDOFs(), nbf = E.GetNbfPerFace();
      const int no = op.NumOwnedFaultDOFs();
      // face centroid and area per local DOF, vertices in a canonical order
      Vector cl(4 * nd);
      for (int k = 0; k < nd; k++)
      {
         long long key[3];
         real_t xv[3][3];
         op.Orientation().FaceVertexData(k / nbf, key, xv);
         std::array<std::array<real_t, 3>, 3> v;
         for (int i = 0; i < 3; i++) { for (int c = 0; c < 3; c++) { v[i][c] = xv[i][c]; } }
         std::sort(v.begin(), v.end());
         real_t e1[3], e2[3];
         for (int c = 0; c < 3; c++)
         {
            cl(4 * k + c) = (v[0][c] + v[1][c] + v[2][c]) / 3.0;
            e1[c] = v[1][c] - v[0][c];
            e2[c] = v[2][c] - v[0][c];
         }
         const real_t cx = e1[1] * e2[2] - e1[2] * e2[1], cy = e1[2] * e2[0] - e1[0] * e2[2],
                      cz = e1[0] * e2[1] - e1[1] * e2[0];
         cl(4 * k + 3) = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
      }
      E.RestrictToOwnedFault(cl, ca_, 4);
      const Vector &X = op.FaultPositions();
      auto dist = [&](int s, int j) {
         const real_t dx = X(3 * j) - st[s].along_strike;
         const real_t dz = std::abs(X(3 * j + 2)) - st[s].down_dip;
         return std::sqrt(dx * dx + dz * dz);
      };
      std::vector<real_t> dmin(ns_, std::numeric_limits<real_t>::max());
      for (int s = 0; s < ns_; s++)
      {
         for (int j = 0; j < no; j++) { dmin[s] = std::min(dmin[s], dist(s, j)); }
      }
      MPI_Allreduce(MPI_IN_PLACE, dmin.data(), ns_, MPITypeMap<real_t>::mpi_type, MPI_MIN,
                    MPI_COMM_WORLD);
      cand_.assign(ns_, {});
      copy_.assign(ns_, {});
      for (int s = 0; s < ns_; s++)
      {
         // P: the lexicographically smallest candidate position over all ranks
         real_t P[3] = {std::numeric_limits<real_t>::max(), 0.0, 0.0};
         for (int j = 0; j < no; j++)
         {
            if (dist(s, j) > dmin[s] + kTie) { continue; }
            cand_[s].push_back(j);
            const real_t *x = &X(3 * j);
            if (std::lexicographical_compare(x, x + 3, P, P + 3)) { for (int c = 0; c < 3; c++) { P[c] = x[c]; } }
         }
         int nr = 0;
         MPI_Comm_size(MPI_COMM_WORLD, &nr);
         std::vector<real_t> all(3 * nr);
         MPI_Allgather(P, 3, MPITypeMap<real_t>::mpi_type, all.data(), 3,
                       MPITypeMap<real_t>::mpi_type, MPI_COMM_WORLD);
         for (int r = 0; r < nr; r++)
         {
            if (std::lexicographical_compare(&all[3 * r], &all[3 * r] + 3, P, P + 3))
            {
               for (int c = 0; c < 3; c++) { P[c] = all[3 * r + c]; }
            }
         }
         for (int j = 0; j < no; j++)
         {
            const real_t *x = &X(3 * j);
            if (std::hypot(x[0] - P[0], x[1] - P[1], x[2] - P[2]) <= kTie) { copy_[s].push_back(j); }
         }
         int nc = static_cast<int>(copy_[s].size());
         MPI_Allreduce(MPI_IN_PLACE, &nc, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         MFEM_VERIFY(nc > 0, "station copy average: no DOF at the point of station " << st[s].name);
         if (g_rank == 0)
         {
            env_.emplace_back(dir + "/" + prefix + "_station_" + st[s].name + "_copies.dat");
            env_.back() << "# station " << st[s].name << ": " << nc << " DG copies at ("
                        << std::setprecision(12) << P[0] << ", " << P[1] << ", " << P[2]
                        << "); the station file holds their face-area-weighted average\n"
                        << "# t ncopies sigma_n_min sigma_n_max tau1_min tau1_max tau2_min tau2_max"
                        << " V1_min V1_max V2_min V2_max  (station columns' units and signs)\n";
         }
      }
   }

   /// Replace the station DOFs' data by the copy average; write the envelope.
   void Apply(real_t t, std::vector<DOFData> &dd)
   {
      constexpr int nf = 8;   // slip1 slip2 V1 V2 tau1 tau2 sigma_n psi
      std::vector<real_t> loc;
      for (int s = 0; s < ns_; s++)
      {
         for (int j : copy_[s])
         {
            const DOFData &q = dd[j];
            const real_t row[5 + nf] = {real_t(s), ca_(4 * j), ca_(4 * j + 1), ca_(4 * j + 2),
                                        ca_(4 * j + 3), q.slip1, q.slip2, q.V1, q.V2,
                                        q.tau1_corr, q.tau2_corr, q.sigma_n_corr, q.psi};
            loc.insert(loc.end(), row, row + 5 + nf);
         }
      }
      int nr = 0, n = static_cast<int>(loc.size());
      MPI_Comm_size(MPI_COMM_WORLD, &nr);
      std::vector<int> cnt(nr), dsp(nr, 0);
      MPI_Gather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
      std::vector<real_t> all;
      if (g_rank == 0)
      {
         for (int r = 1; r < nr; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
         all.resize(dsp[nr - 1] + cnt[nr - 1]);
      }
      MPI_Gatherv(loc.data(), n, MPITypeMap<real_t>::mpi_type, all.data(), cnt.data(), dsp.data(),
                  MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
      std::vector<real_t> avg(ns_ * nf, 0.0);
      if (g_rank == 0)
      {
         const int w = 5 + nf, nrow = static_cast<int>(all.size()) / w;
         std::vector<int> ord(nrow);
         for (int i = 0; i < nrow; i++) { ord[i] = i; }
         std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            return std::lexicographical_compare(&all[w * a], &all[w * a] + 4, &all[w * b], &all[w * b] + 4);
         });
         for (int s = 0; s < ns_; s++)
         {
            real_t wsum = 0.0, lo[nf], hi[nf];
            for (int f = 0; f < nf; f++) { lo[f] = std::numeric_limits<real_t>::max(); hi[f] = -lo[f]; }
            int nc = 0;
            for (int i : ord)
            {
               const real_t *r = &all[w * i];
               if (static_cast<int>(r[0]) != s) { continue; }
               nc++;
               wsum += r[4];
               for (int f = 0; f < nf; f++)
               {
                  avg[s * nf + f] += r[4] * r[5 + f];
                  lo[f] = std::min(lo[f], r[5 + f]);
                  hi[f] = std::max(hi[f], r[5 + f]);
               }
            }
            for (int f = 0; f < nf; f++) { avg[s * nf + f] /= wsum; }
            env_[s] << std::scientific << std::setprecision(10) << t << " " << nc << " "
                    << lo[6] << " " << hi[6] << " " << lo[4] << " " << hi[4] << " " << lo[5] << " "
                    << hi[5] << " " << lo[2] << " " << hi[2] << " " << lo[3] << " " << hi[3] << "\n";
            env_[s].flush();
         }
      }
      MPI_Bcast(avg.data(), ns_ * nf, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
      for (int s = 0; s < ns_; s++)
      {
         const real_t *a = &avg[s * nf];
         for (int j : cand_[s])
         {
            DOFData &q = dd[j];
            q.slip1 = a[0]; q.slip2 = a[1]; q.V1 = a[2]; q.V2 = a[3];
            q.tau1_corr = a[4]; q.tau2_corr = a[5]; q.sigma_n_corr = a[6]; q.psi = a[7];
            q.slip_rate = std::hypot(a[2], a[3]);
         }
      }
   }

private:
   int ns_;
   Vector ca_;                                   // owned: centroid (3) + area
   std::vector<std::vector<int>> cand_, copy_;   // owned DOFs overwritten / averaged
   std::vector<std::ofstream> env_;              // rank 0
};

/// Friction-law factory keyed on [meta] law (plan N-07, Generality).
struct ClosureBundle
{
   std::unique_ptr<DieterichRuinaFriction> fr;
   std::unique_ptr<IFaultClosure> closure;
   std::string name;
};

ClosureBundle MakeClosure(const spatial::SpatialFrictionConfig &cfg)
{
   ClosureBundle B;
   if (cfg.law == spatial::FrictionLawKind::RateState)
   {
      MFEM_VERIFY(cfg.rate_state.has_value(), "[meta] law = rate_state needs "
                  "[friction.rate_state]");
      const auto &rs = *cfg.rate_state;
      B.fr = std::make_unique<DieterichRuinaFriction>(
         DieterichRuinaFriction::Constants{rs.V_0_default, rs.f_0_default,
                                        rs.b_default, rs.Dc_default});
      RateStateClosure::StateLaw law = RateStateClosure::StateLaw::Aging;
      if (rs.state_evolution == spatial::StateEvolutionKind::SlipLawStrongRateWeakening)
      {
         MFEM_ABORT("Nitsche driver: state_evolution = slip-law strong rate "
                    "weakening (TPV104) has no IFaultClosure yet.");
      }
      B.closure = std::make_unique<RateStateClosure>(*B.fr, law);
      B.name = "rate_state (aging law)";
      return B;
   }
   MFEM_ABORT("Nitsche driver: [meta] law = slip_weakening has no "
              "IFaultClosure yet (it comes with TPV205).");
   return B;
}

/// Unit vector helpers.
void Normalize(real_t v[3])
{
   const real_t l = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
   MFEM_VERIFY(l > 0.0, "zero-length vector in [fault_geometry]");
   for (int i = 0; i < 3; i++) { v[i] /= l; }
}

void MakeDirs(const std::string &d)
{
   if (g_rank == 0) { ::mkdir(d.c_str(), 0755); }
   MPI_Barrier(MPI_COMM_WORLD);
}

#ifdef MFEM_USE_HDF5
/// Fault ParaView output, [output] paraview_fault = "hdf5": ONE
/// <output_dir>/fault.vtkhdf per run through the SAFS writer's generic
/// entry point vtkhdf::WriteFaultPackHdf (rank-0 gather, one appended step
/// per frame).  Each owned fault face is split at its edge midpoints into 4
/// sub-triangles; every field is evaluated through the face's own nodal
/// basis at the 3 vertices and 3 edge midpoints (at p = 2 these are exactly
/// the 6 nodes, at p = 1 the midpoints are linear interpolants) and each
/// sub-triangle carries the mean of its 3 corner values as cell data.
/// Every sub-triangle has its own 3 vertices (the writer's layout), so the
/// output stays discontinuous like the DG field.  Owned faces only: a
/// shared fault face is written once.  Sub-triangles are wound so that
/// their normal agrees with the fault frame nhat.  Scope: Chunhui Zhao,
/// 2026-10-14 (P7).
class FaultVTKHDFOutput
{
public:
   FaultVTKHDFOutput(NitscheWaveOperator &op, ParMesh &pmesh, int order,
                     const std::string &dir,
                     ParaViewHDFDataCollection::HDFCompression alg, double param)
      : op_(op), E_(op.Elasticity()), dir_(dir), fe_(order, BasisType::GaussLobatto)
   {
      state_.compression_alg = alg;
      state_.compression_param = param;
      nbf_ = E_.GetNbfPerFace();
      // the face basis must be the one the fault DOFs are nodal values of
      const IntegrationRule &nir = E_.GetFaceQuadrature()->GetNodalRule();
      MFEM_VERIFY(nir.GetNPoints() == nbf_ && fe_.GetDof() == nbf_,
                  "fault ParaView output: face basis has " << fe_.GetDof()
                  << " functions, fault DOFs per face " << nbf_);
      for (int k = 0; k < nbf_; k++)
      {
         MFEM_VERIFY(std::abs(nir.IntPoint(k).x - fe_.GetNodes().IntPoint(k).x) < 1e-12 &&
                     std::abs(nir.IntPoint(k).y - fe_.GetNodes().IntPoint(k).y) < 1e-12,
                     "fault ParaView output: face node " << k << " differs from "
                     "the fault DOF nodal rule");
      }
      // 3 vertices, then the midpoints of edges 01, 12, 20
      const real_t ref[6][2] = {{0, 0}, {1, 0}, {0, 1}, {0.5, 0}, {0.5, 0.5}, {0, 0.5}};
      phi_.SetSize(6, nbf_);
      Vector sh(nbf_);
      for (int c = 0; c < 6; c++)
      {
         IntegrationPoint ip;
         ip.Set2(ref[c][0], ref[c][1]);
         fe_.CalcShape(ip, sh);
         for (int k = 0; k < nbf_; k++) { phi_(c, k) = sh(k); }
      }
      // geometry and a check of the face map against the DOF positions
      Vector Xl, Fl;
      E_.ExpandOwnedToLocalFault(op.FaultPositions(), Xl, 3);
      E_.ExpandOwnedToLocalFault(op.FaultFrames(), Fl, 9);
      const Array<int> &omap = E_.GetOwnedFaultFaceMap();
      const Array<int> &fint = E_.GetFaultInteriorFaces();
      const Array<int> &fsh = E_.GetFaultSharedFaces();
      const int nint = fint.Size();
      local_faces_.resize(E_.GetNumOwnedFaultFaces());
      real_t pos_err = 0.0;
      int nflip = 0;
      for (int o = 0; o < E_.GetNumOwnedFaultFaces(); o++)
      {
         const int fi = omap[o];
         local_faces_[o] = fi;
         const int mf = (fi < nint) ? fint[fi] : pmesh.GetSharedFace(fsh[fi - nint]);
         ElementTransformation *T = pmesh.GetFaceTransformation(mf);
         real_t P[6][3];
         for (int c = 0; c < 6; c++)
         {
            IntegrationPoint ip;
            ip.Set2(ref[c][0], ref[c][1]);
            T->SetIntPoint(&ip);
            Vector x(3);
            T->Transform(ip, x);
            for (int i = 0; i < 3; i++) { P[c][i] = x(i); }
         }
         for (int k = 0; k < nbf_; k++)
         {
            Vector x(3);
            T->Transform(nir.IntPoint(k), x);
            for (int i = 0; i < 3; i++)
            {
               pos_err = std::max(pos_err, std::abs(x(i) - Xl(3 * (fi * nbf_ + k) + i)));
            }
         }
         real_t e1[3], e2[3];
         for (int i = 0; i < 3; i++) { e1[i] = P[1][i] - P[0][i]; e2[i] = P[2][i] - P[0][i]; }
         const real_t nx[3] = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                               e1[0] * e2[1] - e1[1] * e2[0]};
         const real_t *nh = &Fl(9 * (fi * nbf_));
         const bool flip = (nx[0] * nh[0] + nx[1] * nh[1] + nx[2] * nh[2]) < 0.0;
         nflip += flip ? 1 : 0;
         for (int s = 0; s < 4; s++)
         {
            int c0 = kSub[s][0], c1 = kSub[s][1], c2 = kSub[s][2];
            if (flip) { std::swap(c1, c2); }
            const int b = static_cast<int>(pack_.vertices.size());
            for (int c : {c0, c1, c2}) { pack_.vertices.push_back({P[c][0], P[c][1], P[c][2]}); }
            pack_.triangles.push_back({b, b + 1, b + 2});
            corners_.push_back({c0, c1, c2});
         }
      }
      MPI_Allreduce(MPI_IN_PLACE, &pos_err, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                    MPI_COMM_WORLD);
      long long cnt[2] = {E_.GetNumOwnedFaultFaces(), nflip};
      MPI_Allreduce(MPI_IN_PLACE, cnt, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      MFEM_VERIFY(pos_err < 1e-6, "fault ParaView output: the face map does not reproduce "
                  "the fault DOF positions (max error " << pos_err << " m)");
      if (g_rank == 0)
      {
         std::cout << "  fault ParaView output: " << dir_ << "/fault.vtkhdf, " << cnt[0]
                   << " fault faces x 4 sub-triangles (" << cnt[1] << " faces re-wound to"
                   << " nhat); face map vs DOF positions max error " << pos_err << " m\n";
      }
      pack_.field_names = {"slip_dip", "slip_strike", "slip_rate_dip", "slip_rate_strike",
                           "slip_rate", "traction_dip", "traction_strike", "state_variable",
                           "normal_stress", "param_a", "param_Dc"};
      if (op_.HasPrakashClifton()) { pack_.field_names.push_back("normal_stress_pc"); nc_ = 11; }
      pack_.field_arrays.assign(pack_.field_names.size(), std::vector<double>());
   }

   /// Collective.  tau = -(X* + gamma V) and sigma_n_eff as the stations;
   /// at t = 0 (no step yet) the projected background traction.
   void Write(int cycle, real_t t, real_t dt)
   {
      const auto t0 = std::chrono::steady_clock::now();
      const int no = op_.FaultDc().Size();
      const Vector &d = op_.Slip(), &V = op_.SlipRate(), &psi = op_.State();
      const Vector &Fr = op_.FaultFrames();
      Vector Xs = op_.LastShear(), sg = op_.LastSigmaEff();
      const bool have_x = Xs.Size() == 2 * no;
      if (!have_x)
      {
         Vector S;
         op_.BackgroundTraction(op_.GetTime(), S);
         Xs.SetSize(2 * no);
         sg.SetSize(no);
         for (int j = 0; j < no; j++)
         {
            fault_frame::ToLocal(&S(3 * j), 0.0, &Fr(9 * j), &Fr(9 * j + 3), &Fr(9 * j + 6),
                                 Xs(2 * j), Xs(2 * j + 1), sg(j));
         }
      }
      const int nc = nc_;   // sd ss vd vs td ts psi sn a Dc [sigma*]
      Vector own(nc * no), loc;
      for (int j = 0; j < no; j++)
      {
         const real_t *dd = &Fr(9 * j + 3), *ss = &Fr(9 * j + 6);
         real_t *q = &own(nc * j);
         q[0] = fault_frame::Dot(&d(3 * j), dd);
         q[1] = fault_frame::Dot(&d(3 * j), ss);
         q[2] = fault_frame::Dot(&V(3 * j), dd);
         q[3] = fault_frame::Dot(&V(3 * j), ss);
         const real_t g = have_x ? 0.5 * op_.FaultSigmaF()(j) * dt : 0.0;
         q[4] = -(Xs(2 * j) + g * q[2]);
         q[5] = -(Xs(2 * j + 1) + g * q[3]);
         q[6] = psi(j);
         q[7] = sg(j);
         q[8] = op_.FaultA()(j);
         q[9] = op_.FaultDc()(j);
         if (nc == 11) { q[10] = op_.PrakashCliftonSigma()(j); }
      }
      E_.ExpandOwnedToLocalFault(own, loc, nc);
      for (auto &a : pack_.field_arrays) { a.clear(); }
      for (std::size_t o = 0; o < local_faces_.size(); o++)
      {
         const int fi = local_faces_[o];
         real_t cv[6][11];
         for (int c = 0; c < 6; c++)
         {
            for (int m = 0; m < nc; m++)
            {
               real_t v = 0.0;
               for (int k = 0; k < nbf_; k++) { v += phi_(c, k) * loc(nc * (fi * nbf_ + k) + m); }
               cv[c][m] = v;
            }
         }
         for (int s = 0; s < 4; s++)
         {
            const auto &cr = corners_[4 * o + s];
            auto mean = [&](int m) { return (cv[cr[0]][m] + cv[cr[1]][m] + cv[cr[2]][m]) / 3.0; };
            auto speed = [&](int c) { return std::hypot(cv[c][2], cv[c][3]); };
            pack_.field_arrays[0].push_back(mean(0));
            pack_.field_arrays[1].push_back(mean(1));
            pack_.field_arrays[2].push_back(mean(2));
            pack_.field_arrays[3].push_back(mean(3));
            pack_.field_arrays[4].push_back((speed(cr[0]) + speed(cr[1]) + speed(cr[2])) / 3.0);
            for (int m = 4; m < nc; m++) { pack_.field_arrays[m + 1].push_back(mean(m)); }
         }
      }
      const auto t1 = std::chrono::steady_clock::now();
      vtkhdf::WriteFaultPackHdf(state_, dir_, cycle, t, pack_, g_rank, Mpi::WorldSize(),
                                MPI_COMM_WORLD);
      frames_++;
      const auto t2 = std::chrono::steady_clock::now();
      pack_seconds_ += std::chrono::duration<double>(t1 - t0).count();
      seconds_ += std::chrono::duration<double>(t2 - t0).count();
   }

   int Frames() const { return frames_; }
   double Seconds() const { return seconds_; }
   double PackSeconds() const { return pack_seconds_; }

private:
   /// sub-triangles of the 6 corners (0,1,2 vertices; 3,4,5 midpoints of
   /// edges 01, 12, 20), each wound like the face (0,1,2)
   static constexpr int kSub[4][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}};

   NitscheWaveOperator &op_;
   ElasticityDomainOperator<ParMesh> &E_;
   std::string dir_;
   H1_TriangleElement fe_;
   int nbf_ = 0;
   int nc_ = 10;   ///< packed fields per DOF (11 with Prakash-Clifton)
   DenseMatrix phi_;
   std::vector<int> local_faces_;
   std::vector<std::array<int, 3>> corners_;
   vtu::LocalFaultPack pack_;
   vtkhdf::FaultHDFState state_;
   int frames_ = 0;
   double seconds_ = 0.0, pack_seconds_ = 0.0;
};
#endif // MFEM_USE_HDF5

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const int nranks = Mpi::WorldSize();
   auto wall0 = std::chrono::steady_clock::now();

   const std::string cfg_path = ArgStr(argc, argv, "--config", "");
   if (cfg_path.empty())
   {
      if (g_rank == 0)
      {
         std::cerr << "usage: seas_nitsche_dyn_driver --config deck.toml [--mesh m.msh]"
                      " [--order p] [--tfinal s] [--output-dir d] [--prefix name]"
                      " [--station-every n] [--dump-every n] [--log-every n]"
                      " [--station-copies average|first]\n";
      }
      return 1;
   }

   // ------------------------------------------------------------ deck
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(cfg_path);
   const toml::value root = toml::parse(cfg_path);
   const toml::value *t_num = DeckTable(root, "numerics");
   const toml::value *t_fg = DeckTable(root, "fault_geometry");
   const toml::value *t_diag = DeckTable(root, "diagnostics");
   const toml::value *t_fric = DeckTable(root, "friction");
   const toml::value *t_pc = (t_fric && t_fric->contains("prakash_clifton"))
                             ? &t_fric->at("prakash_clifton") : nullptr;

   const std::string mesh_path = ArgStr(argc, argv, "--mesh", cfg.mesh.path);
   const int p = ArgInt(argc, argv, "--order", cfg.mesh.order);
   const real_t tfinal = ArgReal(argc, argv, "--tfinal", cfg.time.tfinal);
   const std::string out_dir = ArgStr(argc, argv, "--output-dir", cfg.output.output_dir);
   const std::string prefix = ArgStr(argc, argv, "--prefix",
                                     cfg.problem.tag.empty() ? "nitsche" : cfg.problem.tag);
   const int station_every = ArgInt(argc, argv, "--station-every", 1);
   const int dump_every = ArgInt(argc, argv, "--dump-every", 0);
   const int log_every = ArgInt(argc, argv, "--log-every", 100);
   const std::string station_copies = ArgStr(argc, argv, "--station-copies", "average");
   MFEM_VERIFY(station_copies == "average" || station_copies == "first",
               "--station-copies " << station_copies << ": expected average or first");
   const real_t cfl = cfg.numerics.cfl;
   const real_t penalty_factor = TomlReal(t_num, "penalty_factor", 1.0);
   const bool energy_on = TomlBool(t_diag, "energy_balance", false);
   const bool pc_on = t_pc != nullptr;
   real_t pc_vstar = 0.0, pc_len = 0.0;
   if (pc_on)
   {
      MFEM_VERIFY(t_pc->is_table(), "[friction.prakash_clifton] must be a table");
      MFEM_VERIFY(t_pc->contains("v_star") && t_pc->contains("length"),
                  "[friction.prakash_clifton] needs both v_star (m/s) and length (m); "
                  "there are no defaults");
      pc_vstar = TomlReal(t_pc, "v_star", 0.0);
      pc_len = TomlReal(t_pc, "length", 0.0);
      MFEM_VERIFY(pc_vstar > 0.0 && pc_len > 0.0, "[friction.prakash_clifton] v_star = "
                  << pc_vstar << ", length = " << pc_len << "; both must be > 0");
   }
   MFEM_VERIFY(p >= 1 && p <= 4, "order p must be in 1..4 (P3 constants); got " << p);
   MFEM_VERIFY(cfl > 0.0 && cfl <= 1.0, "[numerics] cfl (fraction of C(p)) must be in (0,1]");
   MakeDirs(out_dir);

   // [output] ParaView keys, SAFS semantics (spatial_dyn_driver): the master
   // gate paraview_enabled switches the ParaView collections; the fault
   // surface is written as ONE fault.vtkhdf per run.  Every key this driver
   // cannot honour aborts instead of being ignored.
   const auto &ocfg = cfg.output;
   const bool pv_fault = ocfg.paraview_enabled && ocfg.paraview_fault != "off";
   MFEM_VERIFY(ocfg.paraview_fault == "off" || ocfg.paraview_fault == "hdf5"
               || ocfg.paraview_fault == "vtu",
               "[output] paraview_fault = '" << ocfg.paraview_fault
               << "': expected \"off\", \"hdf5\" or \"vtu\"");
   MFEM_VERIFY(!(pv_fault && ocfg.paraview_fault == "vtu"),
               "[output] paraview_fault = \"vtu\" writes one file per snapshot; the "
               "Nitsche driver writes the fault as ONE file per run: use \"hdf5\" "
               "(VTKHDF) or \"off\"");
   MFEM_VERIFY(ocfg.paraview_free_surface == "off",
               "[output] paraview_free_surface = '" << ocfg.paraview_free_surface
               << "': the Nitsche driver writes no free-surface slice; set it to \"off\" "
               "(the parser's default is \"vtu\")");
   if (ocfg.paraview_enabled)
   {
      MFEM_VERIFY(ocfg.paraview_volume == "off" && ocfg.paraview_bulk == "off",
                  "[output] paraview_volume / paraview_bulk: the Nitsche driver writes "
                  "only the fault surface; set both to \"off\"");
      MFEM_VERIFY(!ocfg.paraview_fault_legacy_ascii,
                  "[output] paraview_fault_legacy_ascii: not supported by the Nitsche "
                  "driver (single VTKHDF file only)");
      MFEM_VERIFY(ocfg.paraview_coseismic_dt <= 0.0 && ocfg.paraview_nucleation_dt <= 0.0
                  && ocfg.paraview_interseismic_dt <= 0.0,
                  "[output] paraview_{coseismic,nucleation,interseismic}_dt: the Nitsche "
                  "driver writes at the fixed cadence paraview_fault_dt (or "
                  "paraview_every_steps); remove the regime-adaptive keys");
   }
#ifndef MFEM_USE_HDF5
   MFEM_VERIFY(!pv_fault, "[output] paraview_fault = \"hdf5\" needs MFEM built with "
               "MFEM_USE_HDF5=YES");
#endif

   if (g_rank == 0)
   {
      std::cout << "=== seas_nitsche_dyn_driver (plan N-07) on " << nranks << " rank(s) ===\n"
                << "  deck " << cfg_path << "\n  " << cfg.description << "\n";
      // velocity-stress numerics keys: listed once, ignored
      const char *vs_keys[] = {"ader_order", "time_integrator", "mixed_flux",
                               "interior_flux", "fault_iterator", "use_pml",
                               "cfl_safety", "cfl_dg_safety", "lts"};
      std::string ign;
      for (const char *k : vs_keys) { if (t_num && t_num->contains(k)) { ign += std::string(" ") + k; } }
      if (!ign.empty()) { std::cout << "  ignored [numerics] keys (velocity-stress path):" << ign << "\n"; }
   }

   // ------------------------------------------------------------ mesh
   std::unique_ptr<Mesh> smesh = std::make_unique<Mesh>(mesh_path.c_str(), 1, 1);
   MFEM_VERIFY(smesh->Dimension() == 3, "mesh must be 3-D");
   for (int e = 0; e < smesh->GetNE(); e++)
   {
      MFEM_VERIFY(smesh->GetElementBaseGeometry(e) == Geometry::TETRAHEDRON,
                  "Nitsche driver: tetrahedral meshes only; element " << e << " is "
                  << Geometry::Name[smesh->GetElementBaseGeometry(e)]);
   }
   const long long ne_global = smesh->GetNE();
   ParMesh pmesh(MPI_COMM_WORLD, *smesh);
   smesh.reset();

   // ------------------------------------------------------------ material
   MFEM_VERIFY(cfg.material.kind == spatial::MaterialKind::Constant && !cfg.velocity.use_sidecar,
               "Nitsche driver (P7): only the constant material "
               "([velocity] use_sidecar = false, [material] kind constant) is "
               "wired; the sidecar / depth-profile / halfspace selections are P8.");
   MaterialField material = MaterialField::MakeConstant(
      cfg.material_fallback.lambda, cfg.material_fallback.mu, cfg.material_fallback.rho);

   // ------------------------------------------------------------ boundary
   MFEM_VERIFY(cfg.boundary.fault_attr > 0, "[boundary] fault_attr is required");
   BoundaryConfig bc;
   bc.fault_attr = cfg.boundary.fault_attr;
   for (int a : cfg.boundary.natural_attrs) { bc.natural_attrs.insert(a); }
   for (int a : cfg.boundary.absorbing_attrs) { bc.absorbing_attrs.insert(a); }

   // ------------------------------------------------------------ fault geometry
   NitscheFaultSetup fs;
   {
      const std::string orient = TomlStr(t_fg, "orientation", "propagated");
      const std::string pin_by = TomlStr(t_fg, "pin_by", "reference");
      const std::string tangent = TomlStr(t_fg, "tangent", "pointwise");
      const std::string junction = TomlStr(t_fg, "junction_policy", "refuse");
      MFEM_VERIFY(orient == "propagated", "[fault_geometry] orientation = '" << orient
                  << "': only \"propagated\" is implemented (N-08)");
      MFEM_VERIFY(tangent == "pointwise", "[fault_geometry] tangent = '" << tangent
                  << "': only \"pointwise\" is implemented (N-09; F8 transport is not)");
      MFEM_VERIFY(junction == "refuse", "[fault_geometry] junction_policy = '" << junction
                  << "': only \"refuse\" is implemented");
      auto &oc = fs.orientation;
      for (int i = 0; i < 3; i++)
      {
         oc.ref_normal[i] = cfg.fault_geometry.ref_normal[i];
         oc.up[i] = cfg.fault_geometry.up[i];
         fs.up[i] = cfg.fault_geometry.up[i];
      }
      if (pin_by == "reference") { oc.pin = FaultOrientationConfig::Pin::Reference; }
      else
      {
         MFEM_ABORT("[fault_geometry] pin_by = '" << pin_by << "': only \"reference\" "
                    "is wired in the driver (Pin::Point needs a pin_point key)");
      }
      oc.c_min = TomlReal(t_fg, "pin_min_margin", oc.c_min);
      oc.fault_attrs = {bc.fault_attr};
   }

   // ------------------------------------------------------------ closure
   ClosureBundle law = MakeClosure(cfg);
   fs.closure = law.closure.get();
   fs.sigma_min = cfg.stress.pore_pressure.min_sigma_n_pa;
   MFEM_VERIFY(cfg.stress.pore_pressure.P_p_pa == 0.0
               && cfg.stress.pore_pressure.P_p_grad_pa_per_m == 0.0,
               "Nitsche driver (P7): pore pressure is not wired (P_p must be 0)");

   // ------------------------------------------------------------ prestress (I1)
   MFEM_VERIFY(cfg.stress.kind == spatial::StressSourceKind::ConstantTensor,
               "Nitsche driver: [stress] kind must be \"constant_tensor\" (the global "
               "Cauchy tensor, tension positive; spec I1).  fault_local_prestress is "
               "the slot pattern I2 forbids.");
   const real_t s0[9] = {cfg.stress.sigma_xx_pa, cfg.stress.sigma_xy_pa, cfg.stress.sigma_xz_pa,
                         cfg.stress.sigma_xy_pa, cfg.stress.sigma_yy_pa, cfg.stress.sigma_yz_pa,
                         cfg.stress.sigma_xz_pa, cfg.stress.sigma_yz_pa, cfg.stress.sigma_zz_pa};
   fs.prestress = [s0](const real_t[3], real_t s[9]) { for (int i = 0; i < 9; i++) { s[i] = s0[i]; } };

   // pin frame (F1 with the reference normal): nhat, shat, dhat
   real_t pn[3], pu[3], ps[3], pd[3];
   for (int i = 0; i < 3; i++) { pn[i] = cfg.fault_geometry.ref_normal[i]; pu[i] = cfg.fault_geometry.up[i]; }
   Normalize(pn);
   ff::Cross(pu, pn, ps);
   Normalize(ps);
   ff::Cross(ps, pn, pd);

   // ------------------------------------------------------------ nucleation
   if (cfg.nucleation.enabled)
   {
      MFEM_VERIFY(cfg.nucleation.kind == spatial::NucleationKind::GradualOverstressCompactCircular,
                  "Nitsche driver (P7): only [nucleation] kind = "
                  "\"gradual_overstress_compact_circular\" is wired");
      const auto spec = cfg.nucleation.compact_circular;
      MFEM_VERIFY(spec.radius_m > 0.0 && spec.T_nuc_s > 0.0, "compact-circular nucleation: "
                  "radius_m and T_nuc_s must be > 0");
      // strike overstress delta_tau F(r) G(t) as the global tensor
      //   dsigma = -delta_tau F G (shat (x) nhat + nhat (x) shat),
      // so that -(dsigma nhat).shat = +delta_tau F G.  r is the in-plane
      // distance from the centre measured with (dhat, shat), as the SAFS
      // resolver ResolveGradualOverstressCompactCircular measures it.
      std::array<real_t, 3> c = {spec.center_x_m, spec.center_y_m, spec.center_z_m};
      std::array<real_t, 3> sv = {ps[0], ps[1], ps[2]}, nv = {pn[0], pn[1], pn[2]},
                            dv = {pd[0], pd[1], pd[2]};
      fs.nucleation_spatial = [spec, c, sv, nv, dv](const real_t x[3], real_t s[9]) {
         real_t e[3] = {x[0] - c[0], x[1] - c[1], x[2] - c[2]};
         const real_t ed = e[0] * dv[0] + e[1] * dv[1] + e[2] * dv[2];
         const real_t es = e[0] * sv[0] + e[1] * sv[1] + e[2] * sv[2];
         const real_t amp = -spec.delta_tau_pa
                            * spatial::CompactBellFactor(std::sqrt(ed * ed + es * es), spec.radius_m);
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) { s[3 * i + j] = amp * (sv[i] * nv[j] + nv[i] * sv[j]); }
      };
      const real_t T = spec.T_nuc_s;
      fs.nucleation_temporal = [T](real_t t) { return spatial::SmoothStep(t, T); };
   }

   // ------------------------------------------------------------ friction parameters
   spatial::RateStateBlock rs = *cfg.rate_state;
   rs.eta_auto = false;      // eta is quasi-dynamic radiation damping; the
   rs.eta_default = 0.0;     // Nitsche closure uses gamma = 1/2 sigma_F dt
   spatial::SpatialFrictionResolver resolver;
   const real_t b0 = rs.b_default, f00 = rs.f_0_default, V00 = rs.V_0_default;
   fs.parameters = [&](const real_t x[3], FaultPointState &st) {
      Vector xc(3);
      for (int i = 0; i < 3; i++) { xc(i) = x[i]; }
      Array<int> el(1), at(1);
      el[0] = 0;
      at[0] = bc.fault_attr;
      Vector sn0;
      spatial::RateStatePerDOFParams r = resolver.ResolveRateState(
         rs, xc, el, at, material, pmesh, cfg.stress.pore_pressure, sn0);
      MFEM_VERIFY(r.b(0) == b0 && r.f_0(0) == f00 && r.V_0(0) == V00,
                  "Nitsche driver: per-DOF b / f_0 / V_0 differ from the defaults at ("
                  << x[0] << ", " << x[1] << ", " << x[2] << "); the closure takes "
                  "them as constants (spatially varying b is not wired)");
      st.a = r.a(0);
      st.Dc = r.Dc(0);
   };

   // ------------------------------------------------------------ operator
   DomainConfig dcfg;
   dcfg.penalty_factor = penalty_factor;
   NitscheWaveOperator op(pmesh, p, material, bc, dcfg, &fs);
   MFEM_VERIFY(op.Orientation().GetReport().junction_edges == 0,
               "[fault_geometry] junction_policy = \"refuse\": the fault has "
               << op.Orientation().GetReport().junction_edges << " junction edges");
   const real_t dt_cfl = ComputeMaxDt(op.Sample(), pmesh, p, cfl * kCeig[p]);
   const real_t dt = std::min(dt_cfl, cfg.time.dt_max);
   op.SetTimeStep(dt);

   // initial fault state: d = 0, psi from the continuous friction condition
   // at the projected prestress, V = V_init along -X_tan (C-07)
   const int no = op.NumOwnedFaultDOFs();
   {
      Vector d(3 * no), V(3 * no), psi(no), S0, sig0(no);
      d = 0.0;
      op.BackgroundTraction(cfg.time.t_initial, S0);
      Vector xc(3);
      Array<int> el(1), at(1);
      el[0] = 0; at[0] = bc.fault_attr;
      for (int j = 0; j < no; j++)
      {
         const real_t *n = &op.FaultFrames()(9 * j), *dd = &op.FaultFrames()(9 * j + 3),
                       *ss = &op.FaultFrames()(9 * j + 6);
         real_t X[2], sig;
         ff::ToLocal(&S0(3 * j), 0.0, n, dd, ss, X[0], X[1], sig);
         const real_t xm = std::hypot(X[0], X[1]);
         for (int c = 0; c < 3; c++) { xc(c) = op.FaultPositions()(3 * j + c); }
         Vector sn0;
         spatial::RateStatePerDOFParams r = resolver.ResolveRateState(
            rs, xc, el, at, material, pmesh, cfg.stress.pore_pressure, sn0);
         const real_t Vi = r.V_init(0);
         FaultPointState st;
         st.a = op.FaultA()(j);
         st.Dc = op.FaultDc()(j);
         law.closure->SeedEquilibrium(xm, Vi, sig, st);
         psi(j) = st.psi;
         real_t Vl[2] = {0.0, 0.0};
         if (xm > 0.0) { Vl[0] = -Vi * X[0] / xm; Vl[1] = -Vi * X[1] / xm; }
         ff::Embed(Vl[0], Vl[1], dd, ss, &V(3 * j));
         sig0(j) = sig;
      }
      op.SetFaultState(d, V, psi);
      if (pc_on)
      {
         // sigma*^0 = the initial sigma_n per DOF
         op.SetPrakashClifton(pc_vstar, pc_len);
         op.SetPrakashCliftonState(sig0);
      }
      // initial fault state per owned DOF (generic; e.g. to check psi^0
      // against a benchmark's closed form outside the driver)
      std::ofstream ini(out_dir + "/" + prefix + "_initial_r" + std::to_string(g_rank) + ".csv");
      ini << std::setprecision(17) << "# x y z a Dc psi0 tau0_strike(-T.shat) sigma_n0 V0_strike\n";
      for (int j = 0; j < no; j++)
      {
         const real_t *n = &op.FaultFrames()(9 * j), *dd = &op.FaultFrames()(9 * j + 3),
                       *ss = &op.FaultFrames()(9 * j + 6);
         real_t X[2], sig;
         ff::ToLocal(&S0(3 * j), 0.0, n, dd, ss, X[0], X[1], sig);
         const real_t *x = &op.FaultPositions()(3 * j);
         ini << x[0] << " " << x[1] << " " << x[2] << " " << op.FaultA()(j) << " "
             << op.FaultDc()(j) << " " << psi(j) << " " << -X[1] << " " << sig << " "
             << ff::Dot(&V(3 * j), ss) << "\n";
      }
   }

   // ------------------------------------------------------------ nucleation check
   // The projected SCEC along-strike traction -(S.shat) at the fault DOF
   // nearest the deck's hypocentre against the analytic tau_0 + delta_tau
   // F(r) G(t) at that DOF's position (plan P7 decision 9).
   if (cfg.nucleation.enabled)
   {
      const auto &spec = cfg.nucleation.compact_circular;
      const real_t h[3] = {cfg.hypocenter.x_m, cfg.hypocenter.y_m, cfg.hypocenter.z_m};
      struct { double d; int r; } loc = {1e300, g_rank}, glob;
      int jh = -1;
      for (int j = 0; j < no; j++)
      {
         const real_t *x = &op.FaultPositions()(3 * j);
         const real_t dd = std::hypot(x[0] - h[0], x[1] - h[1], x[2] - h[2]);
         if (dd < loc.d) { loc.d = dd; jh = j; }
      }
      MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
      if (glob.r == g_rank && jh >= 0)
      {
         const real_t *x = &op.FaultPositions()(3 * jh);
         const real_t *ss = &op.FaultFrames()(9 * jh + 6), *dd = &op.FaultFrames()(9 * jh + 3);
         real_t e[3] = {x[0] - spec.center_x_m, x[1] - spec.center_y_m, x[2] - spec.center_z_m};
         const real_t r = std::hypot(ff::Dot(e, dd), ff::Dot(e, ss));
         const real_t F = spatial::CompactBellFactor(r, spec.radius_m);
         Vector S0;
         op.BackgroundTraction(0.0, S0);
         const real_t tau0 = -ff::Dot(&S0(3 * jh), ss);
         std::cout << std::setprecision(8) << "  nucleation check at the fault DOF nearest the hypocentre ("
                   << x[0] << ", " << x[1] << ", " << x[2] << "), " << glob.d << " m away, r = "
                   << r << " m, F(r) = " << F << ":\n     t [s]   projected -T.shat [MPa]   analytic "
                   "tau_0 + dtau F G [MPa]   difference [Pa]\n";
         for (real_t tc : {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5})
         {
            Vector S;
            op.BackgroundTraction(tc, S);
            const real_t proj = -ff::Dot(&S(3 * jh), ss);
            const real_t ana = tau0 + spec.delta_tau_pa * F * spatial::SmoothStep(tc, spec.T_nuc_s);
            std::cout << "     " << std::setw(5) << tc << "   " << std::setw(20) << proj / 1e6
                      << "   " << std::setw(26) << ana / 1e6 << "   " << proj - ana << "\n";
         }
         std::cout << std::flush;
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   // ------------------------------------------------------------ outputs
   // station adapter: owned DOF positions -> TPV102StationWriter
   std::vector<Vector> fcoords(no);
   for (int j = 0; j < no; j++)
   {
      fcoords[j].SetSize(3);
      for (int c = 0; c < 3; c++) { fcoords[j](c) = op.FaultPositions()(3 * j + c); }
   }
   std::unique_ptr<TPV102StationWriter> stations;
   std::unique_ptr<StationCopyAverage> station_avg;
   if (cfg.problem.tag == "tpv102")
   {
      stations = std::make_unique<TPV102StationWriter>();
      stations->Open(out_dir, prefix, DefaultStations(), fcoords, no, MPI_COMM_WORLD);
      if (station_copies == "average")
      {
         station_avg = std::make_unique<StationCopyAverage>(op, DefaultStations(), out_dir, prefix);
      }
      if (g_rank == 0)
      {
         std::ofstream h(out_dir + "/" + prefix + "_station_conventions.txt");
         h << "# Nitsche dynamic-rupture station conventions (spec I3)\n"
           << "# frame: nhat = pin reference normal (" << pn[0] << ", " << pn[1] << ", " << pn[2]
           << "), up = (" << pu[0] << ", " << pu[1] << ", " << pu[2] << ")\n"
           << "#        shat = up x nhat / |up x nhat| = (" << ps[0] << ", " << ps[1] << ", " << ps[2]
           << "), dhat = shat x nhat = (" << pd[0] << ", " << pd[1] << ", " << pd[2] << ")\n"
           << "# pin: orientation propagated, pinned by the reference normal; plus side = side nhat points into\n"
           << "# columns (TPV102StationWriter layout): t slip1 slip2 V1 V2 tau1 tau2 sigma_n log10_theta\n"
           << "#   slip1 = delta.dhat (SCEC v-slip)      slip2 = delta.shat (SCEC h-slip)\n"
           << "#   V1    = V.dhat     (SCEC v-slip-rate) V2    = V.shat     (SCEC h-slip-rate)\n"
           << "#   tau1  = -T.dhat    (SCEC v-shear-stress, Pa)  tau2 = -T.shat (SCEC h-shear-stress, Pa)\n"
           << "#   sigma_n = -T.nhat (compression positive, Pa); T = total traction on the face with normal nhat\n"
           << "#   delta = u(+) - u(-), V = d delta/dt (SCEC slip definition); log10_theta from psi\n";
      }
   }
   std::vector<DOFData> dofdata(no);
   auto fill_dofdata = [&]() {
      const Vector &d = op.Slip(), &V = op.SlipRate(), &psi = op.State();
      Vector Xs = op.LastShear(), sg = op.LastSigmaEff();
      const bool have_x = Xs.Size() == 2 * no;
      if (!have_x)
      {
         // before the first step: no bulk traction yet, the total traction
         // is the projected background (prestress + nucleation at t)
         Vector S;
         op.BackgroundTraction(op.GetTime(), S);
         Xs.SetSize(2 * no);
         sg.SetSize(no);
         for (int j = 0; j < no; j++)
         {
            const real_t *n = &op.FaultFrames()(9 * j), *dd = &op.FaultFrames()(9 * j + 3),
                          *ss = &op.FaultFrames()(9 * j + 6);
            ff::ToLocal(&S(3 * j), 0.0, n, dd, ss, Xs(2 * j), Xs(2 * j + 1), sg(j));
         }
      }
      for (int j = 0; j < no; j++)
      {
         const real_t *dd = &op.FaultFrames()(9 * j + 3), *ss = &op.FaultFrames()(9 * j + 6);
         DOFData &q = dofdata[j];
         q.slip1 = ff::Dot(&d(3 * j), dd);
         q.slip2 = ff::Dot(&d(3 * j), ss);
         q.V1 = ff::Dot(&V(3 * j), dd);
         q.V2 = ff::Dot(&V(3 * j), ss);
         const real_t g = have_x ? 0.5 * op.FaultSigmaF()(j) * dt : 0.0;
         // endpoint traction tangential components in the frame: X* + gamma V
         const real_t X1 = Xs(2 * j) + g * q.V1;
         const real_t X2 = Xs(2 * j + 1) + g * q.V2;
         q.tau1_corr = -X1;
         q.tau2_corr = -X2;
         q.sigma_n_corr = sg(j);
         q.psi = psi(j);
         q.Dc = op.FaultDc()(j);
         q.slip_rate = std::hypot(q.V1, q.V2);
      }
   };
   std::unique_ptr<std::ofstream> dump;
   Vector cf;
   if (dump_every > 0)
   {
      dump = std::make_unique<std::ofstream>(out_dir + "/" + prefix + "_fault_r"
                                             + std::to_string(g_rank) + ".csv");
      *dump << std::setprecision(12);
      // face centroid per owned DOF (keys DG copies uniquely)
      ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
      const int nd = E.GetNumFaultDOFs(), nbf = E.GetNbfPerFace();
      Vector cl(3 * nd);
      for (int k = 0; k < nd; k++)
      {
         long long key[3];
         real_t xv[3][3];
         op.Orientation().FaceVertexData(k / nbf, key, xv);
         for (int c = 0; c < 3; c++) { cl(3 * k + c) = (xv[0][c] + xv[1][c] + xv[2][c]) / 3.0; }
      }
      E.RestrictToOwnedFault(cl, cf, 3);
   }
   std::vector<real_t> trupt(no, -1.0);

   // fault ParaView output (one fault.vtkhdf per run)
#ifdef MFEM_USE_HDF5
   std::unique_ptr<FaultVTKHDFOutput> pv;
#endif
   int pv_every = 0;
   bool pv_capped = false;
   if (pv_fault)
   {
#ifdef MFEM_USE_HDF5
      pv_every = ocfg.paraview_every_steps > 0
                 ? ocfg.paraview_every_steps
                 : std::max(1, static_cast<int>(std::lround(ocfg.paraview_fault_dt / dt)));
      auto alg = ParaViewHDFDataCollection::HDFCompression::Deflate;
      double param = 6.0;   // the writer's default (lossless deflate level 6)
      if (ocfg.paraview_fault_zfp_tol > 0.0)
      {
         ParaViewOutput<ParMesh>::ProbeH5ZZfpPluginOrAbort();
         alg = ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy;
         param = ocfg.paraview_fault_zfp_tol;
      }
      else if (ocfg.paraview_fault_deflate_level >= 0)
      {
         param = ocfg.paraview_fault_deflate_level;
      }
      pv = std::make_unique<FaultVTKHDFOutput>(op, pmesh, p, out_dir, alg, param);
      if (g_rank == 0)
      {
         std::cout << "  fault ParaView output every " << pv_every << " steps (" << pv_every * dt
                   << " s), " << (alg == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy
                                  ? "ZFP accuracy " : "deflate level ") << param
                   << ", at most " << ocfg.max_snapshots << " snapshots\n";
      }
#endif
   }

   std::unique_ptr<NitscheEnergyMonitor> mon;

   const long nsteps = std::lround(std::ceil((tfinal - cfg.time.t_initial) / dt - 1e-9));
   if (g_rank == 0)
   {
      long long nloc = op.Height();
      std::cout << "  mesh " << mesh_path << ": " << ne_global << " tets; p = " << p
                << "; law " << law.name << "; fault faces " << op.Orientation().GetReport().faces
                << ", shared fault faces " << op.NumSharedFaultFacesGlobal() << "\n"
                << "  dt = " << dt << " s (CFL walk: cfl " << cfl << " x C(p) " << kCeig[p]
                << " -> " << dt_cfl << "; dt_max " << cfg.time.dt_max << "), " << nsteps
                << " steps to t = " << tfinal << " s; penalty_factor " << penalty_factor
                << "; energy balance " << (energy_on ? "on" : "off") << "\n";
      if (pc_on)
      {
         std::cout << "  Prakash-Clifton regularization of the law's normal stress: v_star = "
                   << pc_vstar << " m/s, length = " << pc_len << " m (L/v_star = "
                   << pc_len / pc_vstar << " s); sigma*^0 = initial sigma_n\n";
      }
      (void)nloc;
   }

   // ------------------------------------------------------------ loop
   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   real_t t = cfg.time.t_initial;
   op.SetTime(t);
   if (energy_on)
   {
      mon = std::make_unique<NitscheEnergyMonitor>(op);
      Vector A0(n);
      op.Mult(U, Ud, A0);
      mon->Start(mon->Evaluate(U, Ud, A0, t));
   }
   CentralDifferenceSolver ode;
   ode.Init(op);
   if (stations)
   {
      fill_dofdata();
      if (station_avg) { station_avg->Apply(t, dofdata); }
      stations->WriteStep(t, dofdata);
   }
   real_t max_ratio_run = 0.0, max_V_run = 0.0, min_sig_run = 1e300, max_sig_run = -1e300;
   real_t min_pc_run = 1e300, max_pc_run = -1e300;
   for (long s = 1; s <= nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      const NitscheStepRecord &rec = op.LastStep();
      max_ratio_run = std::max(max_ratio_run, rec.max_V_dt_over_Dc);
      max_V_run = std::max(max_V_run, rec.max_V);
      min_sig_run = std::min(min_sig_run, rec.min_sigma_eff);
      max_sig_run = std::max(max_sig_run, rec.max_sigma_eff);
      if (pc_on)
      {
         min_pc_run = std::min(min_pc_run, rec.min_sigma_pc);
         max_pc_run = std::max(max_pc_run, rec.max_sigma_pc);
      }
      MFEM_VERIFY(std::isfinite(rec.max_V), "non-finite slip rate at step " << s);
      const Vector &V = op.SlipRate();
      for (int j = 0; j < no; j++)
      {
         if (trupt[j] < 0.0 && std::sqrt(V(3 * j) * V(3 * j) + V(3 * j + 1) * V(3 * j + 1)
                                         + V(3 * j + 2) * V(3 * j + 2)) > 1e-3) { trupt[j] = t; }
      }
      real_t rho = 0.0;
      NitscheEnergyRecord er;
      if (mon)
      {
         er = mon->Evaluate(U, Ud, op.LastAcceleration(), t);
         rho = mon->Step(er);
      }
      if (stations && s % station_every == 0)
      {
         fill_dofdata();
         if (station_avg) { station_avg->Apply(t, dofdata); }
         stations->WriteStep(t, dofdata);
      }
      if (dump && s % dump_every == 0)
      {
         const Vector &d = op.Slip(), &psi = op.State();
         const Vector &Xs = op.LastShear(), &sg = op.LastSigmaEff();
         const Vector &X = op.FaultPositions();
         for (int j = 0; j < no; j++)
         {
            *dump << s << "," << t << "," << X(3 * j) << "," << X(3 * j + 1) << "," << X(3 * j + 2)
                  << "," << cf(3 * j) << "," << cf(3 * j + 1) << "," << cf(3 * j + 2)
                  << "," << d(3 * j) << "," << d(3 * j + 1) << "," << d(3 * j + 2)
                  << "," << V(3 * j) << "," << V(3 * j + 1) << "," << V(3 * j + 2)
                  << "," << psi(j) << "," << Xs(2 * j) << "," << Xs(2 * j + 1) << "," << sg(j);
            if (pc_on) { *dump << "," << op.PrakashCliftonSigma()(j); }
            *dump << "\n";
         }
      }
#ifdef MFEM_USE_HDF5
      if (pv && s % pv_every == 0)
      {
         // Frames() is the same on every rank, so the collective Write is uniform
         if (pv->Frames() < ocfg.max_snapshots) { pv->Write(static_cast<int>(s), t, dt); }
         else if (!pv_capped)
         {
            pv_capped = true;
            if (g_rank == 0)
            {
               std::cout << "  fault ParaView output: max_snapshots = " << ocfg.max_snapshots
                         << " reached at t = " << t << " s; no further frames\n";
            }
         }
      }
#endif
      if (g_rank == 0 && (s % log_every == 0 || s == nsteps))
      {
         const real_t wall = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - wall0).count();
         std::cout << std::setprecision(6) << "  step " << s << " t = " << t << "  max|V| = "
                   << rec.max_V << "  max|V|dt/Dc = " << rec.max_V_dt_over_Dc
                   << "  sigma_n_eff in [" << rec.min_sigma_eff << ", " << rec.max_sigma_eff << "]";
         if (pc_on)
         {
            std::cout << "  sigma_pc in [" << rec.min_sigma_pc << ", " << rec.max_sigma_pc << "]";
         }
         if (mon) { std::cout << "  E_h = " << er.E_h << "  sum rho = " << mon->SumRho(); }
         std::cout << "  wall " << wall << " s\n" << std::flush;
         (void)rho;
      }
   }
   if (stations) { stations->Close(); }

   // rupture times
   {
      std::ofstream rt(out_dir + "/" + prefix + "_rupture_time_r" + std::to_string(g_rank) + ".csv");
      rt << std::setprecision(12) << "# x y z t_rupture (first |V| > 1e-3 m/s; -1 = never)\n";
      const Vector &X = op.FaultPositions();
      for (int j = 0; j < no; j++)
      {
         rt << X(3 * j) << " " << X(3 * j + 1) << " " << X(3 * j + 2) << " " << trupt[j] << "\n";
      }
   }
#ifdef MFEM_USE_HDF5
   if (pv && g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  fault ParaView output: " << pv->Frames()
                << " frames in " << out_dir << "/fault.vtkhdf; on rank 0 "
                << pv->Seconds() << " s ("
                << (pv->Frames() ? pv->Seconds() / pv->Frames() : 0.0) << " s/frame: pack "
                << (pv->Frames() ? pv->PackSeconds() / pv->Frames() : 0.0)
                << ", gather + HDF write "
                << (pv->Frames() ? (pv->Seconds() - pv->PackSeconds()) / pv->Frames() : 0.0)
                << ")\n";
   }
#endif
   const real_t wall = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - wall0).count();
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "=== done: " << nsteps << " steps, t = " << t
                << " s; max|V| " << max_V_run << " m/s; max|V|dt/Dc " << max_ratio_run
                << "; sigma_n_eff in [" << min_sig_run << ", " << max_sig_run << "] Pa";
      if (pc_on) { std::cout << "; sigma_pc in [" << min_pc_run << ", " << max_pc_run << "] Pa"; }
      std::cout << "; wall " << wall << " s ===\n";
   }
   return 0;
}
