#!/bin/bash
#SBATCH --job-name=xca-@NAME@
#SBATCH --partition=@PARTITION@
#SBATCH --nodes=@NODES@
#SBATCH --ntasks=@NTASKS@
#SBATCH --cpus-per-task=@CPUS_PER_TASK@
#SBATCH --time=@TIME@
#SBATCH --chdir=@RUNDIR@
#SBATCH --output=@RUNDIR@/slurm-%j.out
#SBATCH --error=@RUNDIR@/slurm-%j.out
@GPU_LINE@
@MEM_LINE@
@EXTRA_LINE@

set -uo pipefail
echo "run=@NAME@ commit=@COMMIT@ host=$(hostname) start=$(date -u +%FT%TZ)"

module purge 2>/dev/null || true
for m in @MODULES@; do module load "$m" 2>/dev/null || echo "WARN module $m not loaded"; done

# Compute nodes usually have no outbound network. The depot was already
# instantiated by validate.sh on the login node; force offline so a missing
# package fails loudly here instead of hanging on a download.
export JULIA_DEPOT_PATH="@DEPOT@"
export JULIA_PKG_OFFLINE=true
export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

mkdir -p @RUNDIR@/summary
cd "@RUNDIR@/src/@CASE_DIR@" || { echo "FATAL no case dir @CASE_DIR@"; exit 2; }

rc=0
@SOLVER_CMD@ || rc=$?

{
  echo "exit_code=$rc"
  echo "commit=@COMMIT@"
  echo "end=$(date -u +%FT%TZ)"
  echo "elapsed=${SECONDS}s"
  echo "nodes=@NODES@ ntasks=@NTASKS@ cpus_per_task=@CPUS_PER_TASK@"
} >>@RUNDIR@/summary/result.txt

# small artifacts only; fields stay on scratch
for f in residuals.dat forces.dat convergence.csv monitor.csv error_norms.csv; do
  [[ -f "$f" ]] && tail -n 5000 "$f" >"@RUNDIR@/summary/$f"
done
cp -- *.png @RUNDIR@/summary/ 2>/dev/null || true

echo "exit=$rc"
exit $rc
