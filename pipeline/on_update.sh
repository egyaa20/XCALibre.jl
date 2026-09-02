#!/bin/bash
# Hook script for hpc/poller.sh (TRIGGER_MODE=hook or "both"). Runs when a
# commit lands on `queue` with `[run]` in the commit message.
#
# Per the pipeline's original design (see CLAUDE.md "Launching a run
# (TRIGGER_MODE=hook)"): this script decides everything and calls sbatch;
# the poller owns the guards, the snapshot, and the reporting. Invoked from
# INSIDE a frozen git-archive snapshot of `queue` at the triggering commit
# (not a persistent clone), with these env vars set by poller.sh:
#   CFD_RUNDIR CFD_COMMIT CFD_JOB_REPORT CFD_ACCOUNT CFD_PARTITION
#   CFD_MODULES CFD_DEPOT CFD_DRY_RUN
#
# This only SUBMITS jobs and returns -- submit.jl needs no XCALibre
# precompilation to do that (TOML/SHA stdlib only), so it's fast, well
# inside poller.sh's HOOK_TIMEOUT (default 300s). The actual solve happens
# in the separately-submitted SLURM jobs, asynchronously, using the
# Manifest.toml this snapshot already carries (copied in by poller.sh's
# snapshot() step) against the shared, pre-instantiated depot.
#
# To change which case/sweep this triggers: edit the CASE= line below,
# commit, and push to `queue` with `[run]` in the commit message. That's the
# whole workflow -- no SSH, no VPN.

set -euo pipefail

# Cron doesn't source /etc/profile.d, so Lmod's `module` function may not
# exist here even though poller.sh's own top-level bootstrap made it
# available in ITS process -- bash functions aren't inherited by this
# subshell (only exported env vars are). Same pattern as validate.sh.
if ! declare -F module >/dev/null 2>&1; then
  for f in /etc/profile.d/z00_lmod.sh /etc/profile.d/00-modulepath.sh; do
    [[ -r "$f" ]] && source "$f"
  done
fi
if declare -F module >/dev/null 2>&1 && [[ -n "${CFD_MODULES:-}" ]]; then
  for m in $CFD_MODULES; do
    module load "$m" || { echo "FATAL: could not load module $m"; exit 1; }
  done
fi

# sbatch itself is expected already on PATH here: unlike the `module`
# function, PATH is a plain exported env var and poller.sh's own top-level
# bootstrap (module load slurm, if needed) propagates to this subshell
# normally.
command -v sbatch >/dev/null 2>&1 || { echo "FATAL: sbatch not on PATH"; exit 1; }

export JULIA_DEPOT_PATH="${CFD_DEPOT:-$HOME/.julia}"
export JULIA_PKG_OFFLINE=true   # this snapshot's Manifest must already resolve

HPC="cases/tatsumoto/configs/hpc.toml"

CASES=(
  "cases/tatsumoto/configs/supercritical.toml"
)

for CASE in "${CASES[@]}"; do
  echo "on_update.sh: submitting $CASE for commit ${CFD_COMMIT:-unknown}"

  # No --reprep: that was a one-time flag for the 2026-08-22 relaunch, to
  # stop the live-property-fix sweep from silently reusing warmup
  # checkpoints computed under the pre-fix (frozen-property) solver.
  # Every checkpoint on disk now postdates that fix, so the default
  # (reuse an existing checkpoint if present) is correct again -- hardcoding
  # --reprep here permanently would mean every future push redoes warmup
  # from scratch even when nothing warmup-relevant changed, silently
  # burning hours on every single hook-triggered run.
  args=(cases/tatsumoto/pipeline/submit.jl "$CASE" --hpc "$HPC")
  [[ -n "${CFD_JOB_REPORT:-}" ]] && args+=(--job-report "$CFD_JOB_REPORT")
  [[ "${CFD_DRY_RUN:-0}" == "1" ]] && args+=(--dry-run)

  julia --project=. "${args[@]}"
done
