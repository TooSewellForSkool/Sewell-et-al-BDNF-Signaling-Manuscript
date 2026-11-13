#!/bin/bash
#
# 03_Leiden.py on the HPC.
# Run from the folder that contains expression_matrix_analysis.csv and 03_Leiden.py:
#   sbatch Pipeline2_03_Leiden.slurm 2        (2 = cluster round; leave off for round 1)
# To use another email without editing the file:  sbatch --mail-user=you@virginia.edu ...
#
# ----------------------------------------------------------------------------
# EDIT THESE (the #SBATCH lines cannot use variables, so edit them directly)
#SBATCH --mail-user=jms8rp@virginia.edu
#SBATCH -A zunderlab
# ----------------------------------------------------------------------------
# Resources (unchanged from your v2 files)
#SBATCH --partition=standard
#SBATCH -n 1
#SBATCH --mem=50000
#SBATCH --time=72:00:00
# Notifications: start, finish, and failure
#SBATCH --mail-type=BEGIN,END,FAIL
# Log files named after the job and its ID, so jobs never overwrite each other
#SBATCH --job-name=Pipeline2_03_Leiden
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

# Name of your conda environment. The python path below is built from it
# ($HOME is your own home folder, so no user name is typed here).
CONDA_ENV_NAME="Pipeline2_jms8rp_env_v2"
PYTHON_BIN="$HOME/.conda/envs/${CONDA_ENV_NAME}/bin/python"

if [ ! -x "$PYTHON_BIN" ]; then
  echo "ERROR: python not found at $PYTHON_BIN -- check CONDA_ENV_NAME at the top of this file"
  exit 1
fi

module purge
module load miniforge

# Always run in the folder the job was submitted from
cd "$SLURM_SUBMIT_DIR"
export PYTHONUNBUFFERED=1   # show progress in the log as it happens

CLUSTER_ROUND=${1:-1}   # first argument after the .slurm file name; 1 if not given
echo "Job $SLURM_JOB_ID on $(hostname): 03_Leiden.py, cluster round $CLUSTER_ROUND"
"$PYTHON_BIN" 03_Leiden.py --round "$CLUSTER_ROUND"
