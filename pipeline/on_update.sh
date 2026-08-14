#!/usr/bin/env bash
# TRIGGER_MODE="hook" only. This is your original design: a script in the repo,
# edited from the phone, that decides what to run and calls sbatch itself.
#
# The poller has already: pulled, run validate.sh, and snapshotted this commit
# into $CFD_RUNDIR/src (you are cd'd there). It handed you:
#   CFD_RUNDIR    scratch dir for this commit
#   CFD_COMMIT    short sha
#   CFD_JOB_REPORT  append "name jobid" per submission so the poller tracks it
#   CFD_ACCOUNT CFD_PARTITION CFD_MODULES CFD_DEPOT CFD_DRY_RUN
# Guards (max active jobs, rate limit, DRY_RUN) are enforced by the poller and
# cannot be overridden from here. That is deliberate.
set -euo pipefail

submit() { # submit <name> <sbatch-script>
  local name=$1 script=$2 jid
  if [[ "${CFD_DRY_RUN:-0}" == "1" ]]; then
    echo "DRY_RUN: would submit $name ($script)"; return 0
  fi
  jid=$(sbatch --parsable "$script")
  echo "$name $jid" >>"$CFD_JOB_REPORT"
  echo "submitted $name as $jid"
}

# ---- example: one job per case listed here -------------------------------
CASES=(
  "tgv_re1600:cases/taylor_green:02:00:00:16"
  # "cavity_re1000:cases/cavity:00:30:00:8"
)

for entry in "${CASES[@]}"; do
  IFS=: read -r name dir h m s cpus <<<"$entry"
  walltime="$h:$m:$s"
  [[ -d "$dir" ]] || { echo "skip $name: no $dir"; continue; }

  script="$CFD_RUNDIR/$name.sbatch"
  cat >"$script" <<EOF
#!/bin/bash
#SBATCH --job-name=xca-$name
#SBATCH --partition=$CFD_PARTITION
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$cpus
#SBATCH --time=$walltime
#SBATCH --chdir=$CFD_RUNDIR
#SBATCH --output=$CFD_RUNDIR/slurm-%j.out
set -uo pipefail
module purge 2>/dev/null || true
for mod in $CFD_MODULES; do module load "\$mod" 2>/dev/null || true; done
export JULIA_DEPOT_PATH="$CFD_DEPOT"
export JULIA_PKG_OFFLINE=true
export JULIA_NUM_THREADS=\${SLURM_CPUS_PER_TASK:-1}
mkdir -p $CFD_RUNDIR/summary
cd "$CFD_RUNDIR/src/$dir"
rc=0
julia --project=$CFD_RUNDIR/src run.jl || rc=\$?
echo "exit_code=\$rc commit=$CFD_COMMIT" >>$CFD_RUNDIR/summary/result.txt
exit \$rc
EOF

  submit "$name" "$script"
done
