# jobs

One folder per cluster (`vista/`, `ls6/`, `discovery/`), sbatch scripts committed as they run.
Every script starts with the commit guard of the repository CLAUDE.md and logs the commit.

Storage rule (project rule, 2026-11-16): build in `$HOME` or `$WORK`; every
result (outputs, logs, checkpoints, the sbatch `-o`/`-e` files) goes to
`$SCRATCH`, never `$WORK` or `$HOME`. Scripts pass `-o $SCRATCH/pff/<run>` to
the drivers and `cd` to a `$SCRATCH` run directory before running.
