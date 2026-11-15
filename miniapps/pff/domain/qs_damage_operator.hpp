// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// qs_damage_operator.hpp — the damage problem of one staggered pass.
//
// MOOSE/RACCOON fracture app: ADPFFDiffusion + ADPFFSource with
// psi = alpha Gc/(c0 l) + g psie_active (ADDerivativeParsedMaterial), bounds
// VariableOldValueBounds (lower) and ConstantBounds 1 (upper).  Residual
//
//   R(d) = int  2 Gc l/c0 grad d . grad v  +  (Gc/(c0 l) alpha'(d) + g'(d) D) v
//
// with D the driving energy (psie_active, element averages by default).  For
// p = 2 the source is linear in d, Source(d, D) = Source(0, D) +
// SourceDerivative(0, D) d, so R(d) = A d - b with
//
//   A = kappa K + M[SourceDerivative(0, D)],   b = -int Source(0, D) v.
//
// Other p need a Newton loop and abort here.

#ifndef MFEM_PFF_QS_DAMAGE_OPERATOR_HPP
#define MFEM_PFF_QS_DAMAGE_OPERATOR_HPP

#include "mfem.hpp"
#include "../phasefield/phase_field_functions.hpp"

#include <memory>

namespace mfem
{
namespace pff
{

/// f(D(x)) for a driving-energy field D (any GridFunction, e.g. L2 order 0).
class DrivingCoefficient : public Coefficient
{
public:
   enum class Kind { Source0, SourceDerivative0 };

   DrivingCoefficient(const PhaseFieldModel &pf, const GridFunction &D, Kind kind)
      : pf_(pf), D_(D), kind_(kind) {}

   real_t Eval(ElementTransformation &T, const IntegrationPoint &ip) override
   {
      const real_t D = D_.GetValue(T, ip);
      return kind_ == Kind::Source0 ? -pf_.Source(0.0, D) : pf_.SourceDerivative(0.0, D);
   }

private:
   const PhaseFieldModel &pf_;
   const GridFunction &D_;
   Kind kind_;
};

class QsDamageOperator
{
public:
   QsDamageOperator(ParFiniteElementSpace &fes, const PhaseFieldModel &pf,
                    const IntegrationRule *ir)
      : fes_(fes), pf_(pf), ir_(ir), kappa_(pf.DiffusionCoefficient())
   {
      MFEM_VERIFY(pf.deg.p == 2.0, "QsDamageOperator: only p = 2 (a linear damage "
                  "residual) is implemented; p = " << pf.deg.p);
      ParBilinearForm k(&fes_);
      ConstantCoefficient kappa(kappa_);
      auto *di = new DiffusionIntegrator(kappa);
      if (ir_) { di->SetIntRule(ir_); }
      k.AddDomainIntegrator(di);
      k.Assemble();
      k.Finalize();
      K_.reset(k.ParallelAssemble());
   }

   /// A and b for the driving field D.
   void Assemble(const GridFunction &D)
   {
      DrivingCoefficient react(pf_, D, DrivingCoefficient::Kind::SourceDerivative0);
      DrivingCoefficient src(pf_, D, DrivingCoefficient::Kind::Source0);
      ParBilinearForm m(&fes_);
      auto *mi = new MassIntegrator(react);
      if (ir_) { mi->SetIntRule(ir_); }
      m.AddDomainIntegrator(mi);
      m.Assemble();
      m.Finalize();
      std::unique_ptr<HypreParMatrix> M(m.ParallelAssemble());
      A_.reset(Add(1.0, *K_, 1.0, *M));

      ParLinearForm f(&fes_);
      auto *li = new DomainLFIntegrator(src);
      if (ir_) { li->SetIntRule(ir_); }
      f.AddDomainIntegrator(li);
      f.Assemble();
      b_.SetSize(fes_.GetTrueVSize());
      f.ParallelAssemble(b_);
   }

   const HypreParMatrix &A() const { return *A_; }
   const Vector &b() const { return b_; }

private:
   ParFiniteElementSpace &fes_;
   const PhaseFieldModel &pf_;
   const IntegrationRule *ir_;
   real_t kappa_;
   std::unique_ptr<HypreParMatrix> K_, A_;
   Vector b_;
};

} // namespace pff
} // namespace mfem

#endif
