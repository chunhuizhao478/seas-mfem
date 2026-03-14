# SEAS Miniapp Development Instructions

## Reference Materials

Before implementing any feature, ALWAYS check:

1. **Tandem source code**: `/Users/chunhuizhao/projects/tandem/`
   - `examples/tandem/2d/` - BP1, BP3, BP5 configurations
   - `app/localoperator/` - Poisson, Elasticity operators
   - `app/form/` - DG operators, boundary handling

2. **Documentation**: `/Users/chunhuizhao/projects/mfem/miniapps/seas/document/`
   - `tandem_to_mfem_porting_report.md` - Overall porting strategy
   - `fault_interface_implementation_guide.md` - Fault operator design
   - `bp2_implementation_plan.md` - Phased implementation plan
   - `tandem_boundary_condition_analysis.md` - BC implementation details

3. **SCEC SEAS benchmark specifications** for BP1, BP2, BP3, etc.

## Key Implementation Notes

- Boundary conditions controlled by `BCMode` enum (default: `Tandem`)
  - **FarField** (default): Dirichlet on attrs 1-4 (x=+-Lx, y=+-Ly), Natural on attrs 5-6 (z=0, z=Lz)
  - **XOnly**: Dirichlet on attrs 1-2 (x=+-Lx) only
  - **AllDirichlet**: Legacy (wrong), Dirichlet on all attrs
- Tandem's `boundary_linear=true` is an optimization flag (boundary = f(x)*t), NOT a BC-type selector
- Tandem's BP5 geo: Physical Surface(1) = {top,bottom} → Natural; Surface(5) = far-field → Dirichlet (H25 fix)
