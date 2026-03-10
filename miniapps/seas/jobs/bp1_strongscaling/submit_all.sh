#!/bin/bash
# Submit all BP1 strong scaling jobs

sbatch bp1_strong_np1.sbatch
sbatch bp1_strong_np2.sbatch
sbatch bp1_strong_np4.sbatch
sbatch bp1_strong_np8.sbatch
sbatch bp1_strong_np16.sbatch
sbatch bp1_strong_np32.sbatch
sbatch bp1_strong_np40.sbatch
sbatch bp1_strong_np80.sbatch
