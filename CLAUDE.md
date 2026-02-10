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
