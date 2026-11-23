# Archived job scripts

Older orders, resolutions and test runs, moved here on 2026-11-23 so that
`jobs/vista_nitsche/` holds only the script of each case's SCEC submission run (plus
the build and mesh jobs it needs).  Each script is unchanged.  They source
`$SLURM_SUBMIT_DIR/...` (run_nitsche_vista.sh), so submit them from `jobs/vista_nitsche/`:

    cd miniapps/seas/jobs/vista_nitsche && sbatch archive/<script>.sbatch
