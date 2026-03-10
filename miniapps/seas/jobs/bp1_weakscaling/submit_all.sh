#!/bin/bash
# Submit all BP1 weak scaling jobs

echo "Submitting BP1 weak scaling jobs..."

sbatch bp1_weak_np4.sbatch
sbatch bp1_weak_np8.sbatch
sbatch bp1_weak_np16.sbatch
sbatch bp1_weak_np32.sbatch
sbatch bp1_weak_np40.sbatch
sbatch bp1_weak_np80.sbatch

echo "All weak scaling jobs submitted."
