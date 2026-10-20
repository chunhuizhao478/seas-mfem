// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// rheology_stack.hpp — the ordered set of rheology models active in a run
// (rheology plan rev 1, sections 4.1 and 6).
//
// The models run in a fixed order at each rheology node, attenuation first,
// then plasticity or damage-breakage (SeisSol orders attenuation and
// plasticity the same way); each model receives the correction of the ones
// before it (IRheologyModel::Update).  The allowed sets are
//
//    {drucker_prager}  {viscoelastic}  {viscoelastic, drucker_prager}
//    {damage_breakage} {viscoelastic, damage_breakage}
//
// and any other combination is refused by name.  The state of the models is
// stored back to back per node; StateOffset(k) is model k's first slot.
// An empty stack is linear elasticity: the coupler is never built.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_RHEOLOGY_STACK_HPP
#define MFEM_SEAS_RHEOLOGY_STACK_HPP

#include "rheology_model.hpp"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// Stack order is the enum order.
enum class RheologyKind { Viscoelastic = 0, DruckerPrager = 1, DamageBreakage = 2 };

/// The model's deck table name, [rheology.<name>].
inline const char *RheologyKindName(RheologyKind k)
{
   switch (k)
   {
      case RheologyKind::Viscoelastic: return "viscoelastic";
      case RheologyKind::DruckerPrager: return "drucker_prager";
      case RheologyKind::DamageBreakage: return "damage_breakage";
   }
   return "unknown";
}

class RheologyStack
{
public:
   /// "" when `kinds` (each at most once) is an allowed set, else the reason.
   static std::string CheckCombination(std::vector<RheologyKind> kinds)
   {
      std::sort(kinds.begin(), kinds.end());
      for (size_t i = 1; i < kinds.size(); i++)
      {
         if (kinds[i] == kinds[i - 1])
         {
            return std::string("[rheology] ") + RheologyKindName(kinds[i]) + " appears twice";
         }
      }
      const bool dp = std::count(kinds.begin(), kinds.end(), RheologyKind::DruckerPrager) > 0;
      const bool db = std::count(kinds.begin(), kinds.end(), RheologyKind::DamageBreakage) > 0;
      if (dp && db)
      {
         return "[rheology] drucker_prager and damage_breakage cannot both be enabled: "
                "the allowed sets are {drucker_prager}, {viscoelastic}, "
                "{viscoelastic, drucker_prager}, {damage_breakage}, "
                "{viscoelastic, damage_breakage}";
      }
      return "";
   }

   /// Models must be added in stack order, forming an allowed set.
   void Add(RheologyKind kind, std::unique_ptr<IRheologyModel> model)
   {
      MFEM_VERIFY(model, "RheologyStack::Add: null model");
      MFEM_VERIFY(kinds_.empty() || kind > kinds_.back(),
                  "RheologyStack::Add: " << RheologyKindName(kind) << " added out of order");
      std::vector<RheologyKind> k = kinds_;
      k.push_back(kind);
      const std::string err = CheckCombination(k);
      MFEM_VERIFY(err.empty(), err);
      offset_.push_back(nstate_);
      nstate_ += model->NumStateVars();
      kinds_.push_back(kind);
      models_.push_back(std::move(model));
   }

   bool Empty() const { return models_.empty(); }
   int Size() const { return static_cast<int>(models_.size()); }
   RheologyKind Kind(int k) const { return kinds_[k]; }
   const IRheologyModel &Model(int k) const { return *models_[k]; }
   int StateOffset(int k) const { return offset_[k]; }

   /// State variables per node, all models.
   int NumStateVars() const { return nstate_; }

   /// True when every model is linear in the strain history.
   bool IsLinearInStrain() const
   {
      for (const auto &m : models_) { if (!m->IsLinearInStrain()) { return false; } }
      return true;
   }

   /// e.g. "viscoelastic + drucker_prager"
   std::string Describe() const
   {
      std::string s;
      for (int k = 0; k < Size(); k++) { s += (k ? " + " : "") + std::string(models_[k]->Name()); }
      return s;
   }

   void InitState(const RheologyPointContext &ctx, real_t *q) const
   {
      for (int k = 0; k < Size(); k++) { models_[k]->InitState(ctx, q + offset_[k]); }
   }

   /// All models in order; dsig is the total correction.  True when any
   /// model was inelastic at this node.
   bool Update(const RheologyPointContext &ctx, const real_t eps[6], real_t dt, real_t *q,
               real_t dsig[6]) const
   {
      for (int i = 0; i < 6; i++) { dsig[i] = 0.0; }
      bool inelastic = false;
      for (int k = 0; k < Size(); k++)
      {
         const bool ik = models_[k]->Update(ctx, eps, dt, q + offset_[k], dsig);
         inelastic = inelastic || ik;
      }
      return inelastic;
   }

   real_t StiffnessBound(const RheologyPointContext &ctx, const real_t *q) const
   {
      real_t b = 0.0;
      for (int k = 0; k < Size(); k++)
      {
         b = std::max(b, models_[k]->StiffnessBound(ctx, q + offset_[k]));
      }
      return b;
   }

private:
   std::vector<RheologyKind> kinds_;
   std::vector<std::unique_ptr<IRheologyModel>> models_;
   std::vector<int> offset_;
   int nstate_ = 0;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RHEOLOGY_STACK_HPP
