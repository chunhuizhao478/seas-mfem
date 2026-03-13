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

- Boundary conditions: See `tandem_boundary_condition_analysis.md` for correct BC approach
- Plate loading is applied on **ALL boundaries** (attrs 1-6): u = (0, sgn(x)·Vp·t/2, 0), matching Tandem's `boundary_linear=true`
- No free surface — consistent with Tandem bp5.lua (H10 fix)
