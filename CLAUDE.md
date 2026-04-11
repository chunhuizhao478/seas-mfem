# Project Instructions

## SEAS Miniapp Development

When implementing the SEAS (Sequences of Earthquakes and Aseismic Slip) miniapp:

- **Always refer to Tandem code** at `/Users/chunhuizhao/projects/tandem` for implementation details
- **Always consult documents** in `/Users/chunhuizhao/projects/mfem/miniapps/seas/document/` for:
  - Porting reports and implementation guides
  - Boundary condition analysis
  - BP2 implementation plans
  - Any other reference documentation

Before implementing any feature, check:
1. How Tandem implements it (source code in `/Users/chunhuizhao/projects/tandem/`)
2. Relevant documentation in the seas/document folder
3. SCEC SEAS benchmark specifications

## Environment Setup

- **Building MFEM/SEAS code (mpicxx, MPI):** `conda activate mfem-dev`
- **Gmsh mesh generation:** `conda activate pythonenv`

## Proposing Fixes

Before proposing any fix:
- **Always read the debug documents** in `miniapps/seas/debug_document/bp5_debug_document/` to understand the full history of what has been tried, what worked, and what failed
- **Never revert a previous fix without careful justification.** Each fix in the debug history was made for a specific reason. If you believe a prior fix should be reverted, you must:
  1. Cite the specific debug document and fix you want to revert
  2. Explain why the original reasoning was wrong
  3. Show evidence (data, math, or code analysis) that reverting improves correctness
  4. Get explicit approval before reverting

## Debugging Behavior

When encountering errors or issues:
- Do NOT silently fall back to simpler approaches or workarounds
- Do NOT assume a simpler solution exists without asking
- REPORT the exact error message and what you tried
- ASK before changing the approach or simplifying the implementation
- If something fails, explain WHY it failed before proposing alternatives

## When Things Don't Work

1. Show me the actual error output
2. Explain what you think is causing it
3. Ask me how I want to proceed rather than autonomously pivoting

Do not assume I want a simpler solution. I may need the complex approach to work.
