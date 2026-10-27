# Version control: GitHub first (MANDATORY, Chunhui Zhao 2026-10-15)

GitHub is the single source of truth for every version that runs on a cluster
(Frontera, Discovery, condo, any HPC).  BP5 fully dynamic and all later work
build on these runs, so a cluster run must always be reproducible from a
GitHub commit.

1. Commit and push FIRST, then pull on the cluster.  Every change to code,
   Makefiles, input decks (`*.toml`) and job scripts (`*.sbatch`) that a
   cluster job will use is committed and pushed to GitHub before the cluster
   sees it.  The cluster checkout gets it with `git fetch` + `git checkout` /
   `git pull` of that pushed commit.
2. Never patch, scp or edit source on the cluster: no `git apply` of patches,
   no copied `.cpp/.hpp/Makefile`, no `git stash` there.  A cluster checkout
   has no uncommitted tracked changes when it builds or runs.
3. Decks and job scripts live in the repo (e.g. `miniapps/seas/<problem>/configs/`)
   and are committed as they run; the job reads them from the checkout.
4. Every build/run job script starts with this guard and logs the commit:
   ```bash
   cd $ROOT
   git fetch -q origin || { echo "ERROR: git fetch failed"; exit 2; }
   [ -z "$(git status --porcelain --untracked-files=no)" ] \
      || { echo "ERROR: uncommitted changes in $ROOT"; git status --short; exit 2; }
   [ -n "$(git branch -r --contains HEAD)" ] \
      || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not on GitHub"; exit 2; }
   echo "commit $(git rev-parse HEAD)"
   ```
   A binary must be built by a job that passed this guard at the same commit.
5. Results record the commit hash (job stdout / driver log); comparisons cite it.
6. Exceptions (e.g. a throw-away debug build) need Chunhui Zhao's explicit
   approval in the conversation, and such runs are never used as results.

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
