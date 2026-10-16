// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_cycle_api_pin.cpp — fully dynamic BP5 cycle, rule M5 (plan rev 3).
//
// Pins every shared function and type the two adapters will call
// (cycle/qd_phase.hpp, cycle/dyn_phase.hpp; plan Sections 3.2-3.3), so that a
// merge which changes one of them fails here, at compile time, on the line
// that names it, instead of deep inside an adapter.  Each pin is an
// unevaluated call (decltype) with the argument types the adapter passes, or
// a constructibility check; nothing runs.  The binary only reports the count.
//
// Hooks are pinned when the cycle starts to call them
// (document/fd_cycle_dev/HOOKS.md).  P1 calls H1a (FaultGeometry 4-argument
// constructor, SetRateStatePerDOF, GetVwValues) and H1b
// (SetInitialFieldsPerDOF), committed on this branch (they are not on the
// base 010b56e).  H2 (SetBackgroundTraction) is pinned in P2.
//
// Build and run (from miniapps/seas):
//   make -f cycle/cycle.mk test-cycle-api-pin

#include "mfem.hpp"

// Adapter 1 (QDPhase): the seas_bp5_full stack.
#include "../../solver/seas_operator.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp5_params.hpp"
#include "../../io/bp5_benchmark_output.hpp"
#include "../../io/checkpoint.hpp"
#include "../../io/petsc_ts_checkpoint.hpp"
#include "../../common/mpi_context.hpp"

// Adapter 2 (DynPhase) and the state bridge: the Nitsche dynamic stack.
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../domain/domain_config.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../friction/rate_state_closure.hpp"
#include "../../fault/fault_frame.hpp"
#include "../../fault/fault_orientation.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../spatial/code/spatial_friction.hpp"

#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

using namespace mfem;
using namespace mfem::seas;

namespace
{

template <class T> T &LV();                // an lvalue of T (never defined)
template <class T> const T &CLV();         // a const lvalue of T

constexpr int kPinBase = __COUNTER__;

// The expression must compile.
#define PIN(...)                                                              \
   static_assert(sizeof(decltype((void)(__VA_ARGS__), 0)) > 0, #__VA_ARGS__); \
   static_assert(__COUNTER__ >= 0, "")
// The expression must compile and have exactly type T.  A data member is
// written unparenthesized, so decltype gives its declared type.
#define PIN_TYPE(T, ...)                                                      \
   static_assert(std::is_same_v<decltype(__VA_ARGS__), T>,                    \
                 "type of " #__VA_ARGS__ " is no longer " #T);                \
   static_assert(__COUNTER__ >= 0, "")
// T must be constructible from the listed argument types.
#define PIN_CTOR(T, ...)                                                      \
   static_assert(std::is_constructible_v<T, __VA_ARGS__>,                     \
                 #T " is no longer constructible from (" #__VA_ARGS__ ")");   \
   static_assert(__COUNTER__ >= 0, "")

// ===================================================================
// Adapter 1, QDPhase (plan 3.3 d; template bp5_verification_full.cpp
// :1549-1652, PETSc :2479-2520)
// ===================================================================
using QDDomain = ElasticityDomainOperator<ParMesh>;
using QDGeom = FaultGeometry<ParMesh>;
using QDFault = RateStateFaultOperator<ParMesh, 2>;
using QDOp = PBP5SEASOp;

static_assert(std::is_same_v<QDOp, SEASQuasiDynamicOperator<ParMesh, QDDomain, QDFault>>,
              "PBP5SEASOp");

// Domain: the legacy constructor (Dirichlet {5}, natural {1}, fault 3).
PIN_CTOR(QDDomain, ParMesh &, int, real_t, real_t, real_t, real_t, real_t,
         DGMethod, SolverType, BCMode, int);
PIN(DGMethod::IP);
PIN(SolverType::MUMPS);
PIN(SolverType::MUMPS_BLR);
PIN(BCMode::FarField);
PIN(BasisType::GaussLobatto);
PIN(LV<QDDomain>().SetPenaltyFactor(real_t()));
PIN(LV<QDDomain>().SetBLRTol(real_t()));
PIN_TYPE(const Array<int> &, CLV<QDDomain>().GetFaultInteriorFaces());
PIN_TYPE(const Array<int> &, CLV<QDDomain>().GetFaultSharedFaces());
PIN(CLV<QDDomain>().RestrictToOwnedFault(CLV<Vector>(), LV<Vector>()));
PIN_TYPE(int, CLV<QDDomain>().GetNumFaultDOFs());
PIN_TYPE(int, CLV<QDDomain>().GetNumOwnedFaultDOFs());
PIN_TYPE(int, CLV<QDDomain>().GetNbfPerFace());
PIN(CLV<QDDomain>().GetFaultDOFBasis(LV<DenseMatrix>()));

// The seed parameters the legacy constructor and the kernel read.
PIN_TYPE(real_t, CLV<BP5Params>().lambda());
PIN_TYPE(real_t, CLV<BP5Params>().mu());
PIN_TYPE(real_t, CLV<BP5Params>().eta());
PIN(CLV<BP5Params>().Vp);
PIN(CLV<BP5Params>().Wf);
PIN(CLV<BP5Params>().lf);
PIN(CLV<BP5Params>().V0);
PIN(CLV<BP5Params>().f0);
PIN(CLV<BP5Params>().b);
PIN(CLV<BP5Params>().L0);
PIN(CLV<BP5Params>().cs);
PIN(CLV<BP5Params>().rho);

// Fault geometry: the BP5 constructor (the parity test's reference) and the
// per-DOF arrays QDPhase and the station writer read.
PIN_CTOR(QDGeom, DomainOperator<ParMesh> &, const BP5Params &, MPIContext *);
PIN_TYPE(int, CLV<QDGeom>().NumFaultDOFs());
PIN_TYPE(int, CLV<QDGeom>().NumGlobalFaultDOFs());
PIN_TYPE(const Vector &, CLV<QDGeom>().fault_dof_coords_3d());
PIN_TYPE(const Vector &, CLV<QDGeom>().GetCoordsX2());
PIN_TYPE(const Vector &, CLV<QDGeom>().GetCoordsX3());
PIN_TYPE(const Vector &, CLV<QDGeom>().GetAValues());
PIN_TYPE(const Vector &, CLV<QDGeom>().GetDcValues());
PIN_TYPE(const Vector &, CLV<QDGeom>().GetEtaValues());
PIN_TYPE(const Vector &, CLV<QDGeom>().GetTauPre());
PIN_TYPE(const Vector &, CLV<QDGeom>().GetVInit());
PIN_TYPE(const Vector &, CLV<QDGeom>().sigma_n_per_dof());

// Hook H1a: the constructor without the analytic BP5 fields, and the per-DOF
// rate-and-state fill from the resolver.
PIN_CTOR(QDGeom, DomainOperator<ParMesh> &, const BP5Params &, MPIContext *, bool);
PIN_TYPE(void, LV<QDGeom>().SetRateStatePerDOF(CLV<spatial::RateStatePerDOFParams>(),
                                               CLV<Vector>()));
PIN_TYPE(const Vector &, CLV<QDGeom>().GetVwValues());
// Hook H1b: per-DOF tau0, V_init and sigma_n.
PIN_TYPE(void, LV<QDGeom>().SetInitialFieldsPerDOF(CLV<Vector>(), CLV<Vector>(),
                                                   CLV<Vector>()));

// Friction law and aging law.
PIN_CTOR(DieterichRuinaFriction, const DieterichRuinaFriction::Constants &);
PIN(CLV<DieterichRuinaFriction::Constants>().V0);
PIN(CLV<DieterichRuinaFriction::Constants>().f0);
PIN(CLV<DieterichRuinaFriction::Constants>().b);
PIN(CLV<DieterichRuinaFriction::Constants>().Dc);
PIN_CTOR(AgingLawPsi, real_t, real_t, real_t);
static_assert(std::is_base_of_v<StateEvolution, AgingLawPsi>, "AgingLawPsi");

// Fault kernel.
PIN_CTOR(QDFault, QDGeom *, DieterichRuinaFriction *, StateEvolution *,
         const BP5Params &, MPIContext *);
PIN(LV<QDFault>().SetScecPsiInit(bool()));
PIN_TYPE(int, CLV<QDFault>().StateSize());
PIN_TYPE(int, CLV<QDFault>().TractionSize());
PIN_TYPE(int, CLV<QDFault>().NumNodes());
static_assert(QDFault::StatePerNode == 3 && QDFault::PsiIndex == 2,
              "QD kernel state layout per node: slip dip, slip strike, psi");
PIN(LV<QDFault>().PreInit(LV<Vector>()));
PIN_TYPE(real_t, LV<QDFault>().Init(CLV<Vector>(), LV<Vector>()));
PIN(CLV<QDFault>().GetSlip(CLV<Vector>(), LV<Vector>()));
PIN(CLV<QDFault>().GetTheta(CLV<Vector>(), LV<Vector>()));
PIN(CLV<QDFault>().RecomputeSlipRate(CLV<Vector>(), CLV<Vector>(), LV<Vector>()));
PIN_TYPE(const Vector &, CLV<QDFault>().GetSlipRate());
PIN(LV<QDFault>().SetSlipRate(CLV<Vector>()));

// Station locator of seas_bp5_full's station writer.
PIN_CTOR(Probe2DInterpolator, const Vector &, const Vector &,
         const std::vector<Probe2DInterpolator::Station> &, int, int, bool);
PIN_TYPE(bool, CLV<Probe2DInterpolator>().HasExactMatch(int()));
PIN_TYPE(real_t, CLV<Probe2DInterpolator>().GetMatchDistance(int()));
PIN_TYPE(int, CLV<Probe2DInterpolator>().GetNearestDOF(int()));
PIN_TYPE(int, CLV<Probe2DInterpolator>().GetFaceStart(int()));
PIN_TYPE(const Vector &, CLV<Probe2DInterpolator>().GetExactWeights(int()));

// Quasi-dynamic operator.
PIN_CTOR(QDOp, QDDomain *, QDFault *, MPIContext *);
static_assert(std::is_base_of_v<TimeDependentOperator, QDOp>, "QD operator");
PIN(LV<QDOp>().SetInitialCondition(LV<Vector>()));
PIN(CLV<QDOp>().Mult(CLV<Vector>(), LV<Vector>()));
PIN(LV<QDOp>().SetElasticSigmaN(bool()));
PIN(LV<QDOp>().SetDisplacement(CLV<Vector>()));
PIN_TYPE(const QDOp::GridFuncType &, CLV<QDOp>().GetDisplacement());
PIN_TYPE(const Vector &, CLV<QDOp>().GetTraction());
PIN_TYPE(const Vector &, CLV<QDOp>().GetNormalTraction());
PIN_TYPE(real_t, CLV<QDOp>().GetMaxSlipRate());
PIN_TYPE(const QDDomain *, CLV<QDOp>().GetDomain());
PIN_TYPE(const QDFault *, CLV<QDOp>().GetFault());

// MPI context.
PIN_CTOR(MPIContext, int *, char ***);
PIN_TYPE(MPI_Comm, CLV<MPIContext>().GetComm());
PIN(CLV<MPIContext>().Rank());
PIN(CLV<MPIContext>().IsRoot());

#ifdef MFEM_USE_PETSC
// PETSc TS (MFEM wrapper; the TS options come from the Tandem options file).
PIN_CTOR(PetscODESolver, MPI_Comm, const std::string &);
PIN(LV<PetscODESolver>().SetAbsTol(real_t()));
PIN(LV<PetscODESolver>().SetRelTol(real_t()));
PIN(LV<PetscODESolver>().SetMaxIter(int()));
PIN(LV<PetscODESolver>().Run(LV<Vector>(), LV<real_t>(), LV<real_t>(), real_t()));
PIN(LV<PetscODESolver>().Init(LV<TimeDependentOperator>(),
                              PetscODESolver::ODE_SOLVER_GENERAL));
PIN_TYPE(petsc::TS, static_cast<petsc::TS>(CLV<PetscODESolver>()));
#endif

// Checkpoints (cycle/checkpoint.hpp calls these unchanged, then appends its
// own block).
PIN_TYPE(void, WriteCheckpoint(CLV<std::string>(), real_t(), real_t(), int(), int(),
                               bool(), CLV<Vector>(), CLV<Vector>(), CLV<Vector>(),
                               CLV<Vector>(), bool(), CLV<Vector>(),
                               static_cast<const MPIContext *>(nullptr)));
PIN_TYPE(bool, ReadCheckpoint(CLV<std::string>(), LV<real_t>(), LV<real_t>(), LV<int>(),
                              LV<int>(), LV<bool>(), LV<Vector>(), LV<Vector>(),
                              LV<Vector>(), LV<Vector>(), LV<bool>(), LV<Vector>(),
                              static_cast<const MPIContext *>(nullptr)));
#ifdef MFEM_USE_PETSC
PIN_TYPE(void, WritePetscTSCheckpoint(CLV<std::string>(), real_t(), real_t(), int(), int(),
                                      int(), real_t(), real_t(), int(), int(), real_t(),
                                      static_cast<const MPIContext *>(nullptr)));
PIN_TYPE(bool, ReadPetscTSCheckpoint(CLV<std::string>(), LV<real_t>(), LV<real_t>(),
                                     LV<int>(), LV<int>(), LV<int>(), LV<real_t>(),
                                     LV<real_t>(), LV<int>(), LV<int>(), LV<real_t>(),
                                     static_cast<const MPIContext *>(nullptr)));
#endif

// ===================================================================
// Adapter 2, DynPhase (plan 3.3 e; template drivers/nitsche_dyn_driver.cpp
// :700-930, :1180)
// ===================================================================
using W = NitscheWaveOperator;

PIN_CTOR(W, ParMesh &, int, const MaterialField &, const BoundaryConfig &,
         const DomainConfig &, const NitscheFaultSetup *);
static_assert(std::is_base_of_v<SecondOrderTimeDependentOperator, W>, "W");
PIN(LV<W>().SetTimeStep(real_t()));
PIN_TYPE(real_t, CLV<W>().GetTimeStep());
PIN(LV<W>().SetTime(real_t()));
PIN(LV<W>().ImplicitSolve(real_t(), real_t(), CLV<Vector>(), CLV<Vector>(), LV<Vector>()));
PIN(CLV<W>().Mult(CLV<Vector>(), CLV<Vector>(), LV<Vector>()));
PIN(LV<W>().SetFaultState(CLV<Vector>(), CLV<Vector>(), CLV<Vector>()));
PIN_TYPE(int, CLV<W>().NumOwnedFaultDOFs());
PIN_TYPE(const Vector &, CLV<W>().Slip());
PIN_TYPE(const Vector &, CLV<W>().SlipRate());
PIN_TYPE(const Vector &, CLV<W>().State());
PIN_TYPE(const Vector &, CLV<W>().FaultPositions());
PIN_TYPE(const Vector &, CLV<W>().FaultFrames());
PIN_TYPE(const Vector &, CLV<W>().FaultA());
PIN_TYPE(const Vector &, CLV<W>().FaultDc());
PIN_TYPE(const Vector &, CLV<W>().PrestressTraction());
PIN(CLV<W>().BackgroundTraction(real_t(), LV<Vector>()));
PIN_TYPE(QDDomain &, LV<W>().Elasticity());
PIN_TYPE(const FaultOrientation &, CLV<W>().Orientation());
PIN(CLV<W>().Orientation().GetReport().junction_edges);
PIN_TYPE(const NitscheMaterialSample &, CLV<W>().Sample());
PIN(CLV<W>().LastStep().max_V);
PIN_TYPE(int, CLV<W>().NumSharedFaultFacesGlobal());
PIN(LV<W>().SetPrakashClifton(real_t(), real_t()));
PIN(LV<W>().SetPrakashCliftonState(CLV<Vector>()));
PIN_TYPE(const Vector &, CLV<W>().PrakashCliftonSigma());
PIN(LV<W>().SetKelvinVoigt(real_t(), W::KVScale::Node, real_t(), W::KVLength::Inscribed,
                           static_cast<const BoundaryConfig *>(nullptr)));

// Fault setup and the pointwise closure state.
PIN_TYPE(const IFaultClosure *, CLV<NitscheFaultSetup>().closure);
PIN_TYPE(FaultOrientationConfig, CLV<NitscheFaultSetup>().orientation);
PIN(CLV<NitscheFaultSetup>().up[2]);
PIN_TYPE(std::function<void(const real_t *, real_t *)>,
         CLV<NitscheFaultSetup>().prestress);
PIN_TYPE(std::function<void(const real_t *, FaultPointState &)>,
         CLV<NitscheFaultSetup>().parameters);
PIN_TYPE(real_t, CLV<NitscheFaultSetup>().sigma_min);
PIN(CLV<FaultOrientationConfig>().ref_normal[2]);
PIN(CLV<FaultOrientationConfig>().up[2]);
PIN(CLV<FaultOrientationConfig>().c_min);
PIN(CLV<FaultOrientationConfig>().fault_attrs);
PIN(FaultOrientationConfig::Pin::Reference);
PIN_TYPE(real_t, CLV<FaultPointState>().psi);
PIN_TYPE(real_t, CLV<FaultPointState>().a);
PIN_TYPE(real_t, CLV<FaultPointState>().Dc);
PIN_CTOR(RateStateClosure, const DieterichRuinaFriction &, RateStateClosure::StateLaw);
PIN(RateStateClosure::StateLaw::Aging);
static_assert(std::is_base_of_v<IFaultClosure, RateStateClosure>, "RateStateClosure");
PIN(CLV<IFaultClosure>().SeedEquilibrium(real_t(), real_t(), real_t(),
                                         LV<FaultPointState>()));
PIN_TYPE(real_t, CLV<IFaultClosure>().StateRate(real_t(), CLV<FaultPointState>()));

// Material, domain configuration, boundaries.
PIN_TYPE(MaterialField, MaterialField::MakeConstant(real_t(), real_t(), real_t()));
PIN(LV<DomainConfig>().penalty_factor = real_t());
PIN_TYPE(std::set<int>, LV<BoundaryConfig>().dirichlet_attrs);
PIN_TYPE(std::set<int>, LV<BoundaryConfig>().natural_attrs);
PIN_TYPE(std::set<int>, LV<BoundaryConfig>().absorbing_attrs);
PIN_TYPE(int, LV<BoundaryConfig>().fault_attr);
using DirichletMap = std::map<int, DirichletFunc>;
PIN_TYPE(DirichletMap, LV<BoundaryConfig>().dirichlet_funcs);
PIN_TYPE(DirichletFunc, LV<BoundaryConfig>().default_dirichlet_func);
static_assert(std::is_same_v<DirichletFunc,
              std::function<void(const Vector &, real_t, Vector &)>>, "DirichletFunc");
PIN_TYPE(DirichletFunc, MakeBP5DirichletFunc(real_t()));

// Explicit time step and stepper.
PIN_TYPE(real_t, ComputeMaxDt(CLV<NitscheMaterialSample>(), LV<ParMesh>(), int(), real_t()));
PIN(LV<CentralDifferenceSolver>().Init(LV<SecondOrderTimeDependentOperator>()));
PIN(LV<CentralDifferenceSolver>().Step(LV<Vector>(), LV<Vector>(), LV<real_t>(),
                                       LV<real_t>()));

// ===================================================================
// State bridge and fault fields (plan 3.3 b, c, f)
// ===================================================================
namespace ff = mfem::seas::fault_frame;
PIN_TYPE(real_t, ff::EpsilonF(static_cast<const real_t *>(nullptr),
                              static_cast<const real_t *>(nullptr)));
PIN_TYPE(void, ff::ToKernel(static_cast<const real_t *>(nullptr), real_t(),
                            static_cast<const real_t *>(nullptr),
                            static_cast<const real_t *>(nullptr),
                            LV<real_t>(), LV<real_t>()));
PIN_TYPE(void, ff::FromKernel(real_t(), real_t(), real_t(), real_t(),
                              static_cast<const real_t *>(nullptr),
                              static_cast<const real_t *>(nullptr),
                              static_cast<const real_t *>(nullptr),
                              static_cast<real_t *>(nullptr)));
PIN_TYPE(void, ff::ToLocal(static_cast<const real_t *>(nullptr), real_t(),
                           static_cast<const real_t *>(nullptr),
                           static_cast<const real_t *>(nullptr),
                           static_cast<const real_t *>(nullptr),
                           LV<real_t>(), LV<real_t>(), LV<real_t>()));
PIN_TYPE(void, ff::Embed(real_t(), real_t(), static_cast<const real_t *>(nullptr),
                         static_cast<const real_t *>(nullptr),
                         static_cast<real_t *>(nullptr)));

using spatial::SpatialFrictionResolver;
using spatial::RateStatePerDOFParams;
PIN(SpatialFrictionResolver());
PIN_TYPE(RateStatePerDOFParams,
         CLV<SpatialFrictionResolver>().ResolveRateState(
            CLV<spatial::RateStateBlock>(), CLV<Vector>(), CLV<Array<int>>(),
            CLV<Array<int>>(), CLV<MaterialField>(), LV<ParMesh>(),
            CLV<spatial::PorePressureSpec>(), CLV<Vector>()));
PIN_TYPE(Vector, LV<RateStatePerDOFParams>().a);
PIN_TYPE(Vector, LV<RateStatePerDOFParams>().b);
PIN_TYPE(Vector, LV<RateStatePerDOFParams>().Dc);
PIN_TYPE(Vector, LV<RateStatePerDOFParams>().V_init);
PIN_TYPE(Vector, LV<RateStatePerDOFParams>().f_0);
PIN_TYPE(Vector, LV<RateStatePerDOFParams>().V_0);
PIN_TYPE(Vector, LV<RateStatePerDOFParams>().eta);
PIN_TYPE(Vector, LV<RateStatePerDOFParams>().sigma_n_eff);
PIN_TYPE(Vector, LV<RateStatePerDOFParams>().V_w);
#ifdef SEAS_USE_TOML
PIN_TYPE(spatial::SpatialFrictionConfig,
         spatial::LoadSpatialFrictionConfig(CLV<std::string>()));
#endif

constexpr int kPinCount = __COUNTER__ - kPinBase - 1;

} // namespace

int main()
{
   std::cout << "PASSED: test_cycle_api_pin: " << kPinCount
             << " API pins compiled (base 010b56e with hooks H1a, H1b; H2 not pinned)\n";
   return 0;
}
