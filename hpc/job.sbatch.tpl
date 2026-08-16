#!/bin/bash
#SBATCH --job-name=xca-@NAME@
#SBATCH --partition=@PARTITION@
@ACCOUNT_LINE@
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

# Post-process: on success, render the pressure field from the last .vtk
# written (off-screen, via the pyvista venv set up separately on the login
# node — compute nodes have no network to install anything). On failure,
# write a text summary instead. Either way this must not flip the job's
# own exit code.
if (( rc == 0 )); then
  last_vtk=$(ls -t -- *.vtk 2>/dev/null | head -1)
  if [[ -n "$last_vtk" ]]; then
    # Load individually: these two carry different GCCcore versions and a
    # combined load can fail as a pair, taking down the one that would
    # otherwise have worked.
    for pm in python-uoneasy/3.12.3-GCCcore-13.3.0 mesa-uoneasy/23.1.9-GCCcore-13.2.0; do
      module load "$pm" 2>/dev/null || echo "WARN postprocess module $pm not loaded"
    done
    if [[ -x "$HOME/cfd-pipeline/venv-pyvista/bin/python3" ]]; then
      "$HOME/cfd-pipeline/venv-pyvista/bin/python3" \
        "@RUNDIR@/src/hpc/postprocess_vtk.py" "$last_vtk" "@RUNDIR@/summary" \
        --field p --decimate 0.75 \
        || echo "WARN postprocess_vtk.py failed on $last_vtk"
    else
      echo "WARN no pyvista venv at \$HOME/cfd-pipeline/venv-pyvista, skipping postprocess"
    fi
  else
    echo "WARN no .vtk file found for postprocessing"
  fi
else
  {
    echo "run=@NAME@ FAILED exit_code=$rc commit=@COMMIT@"
    echo "--- tail of slurm log ---"
    tail -n 100 "@RUNDIR@/slurm-${SLURM_JOB_ID:-unknown}.out" 2>/dev/null
  } >@RUNDIR@/summary/error.txt
fi

echo "exit=$rc"
exit $rc
