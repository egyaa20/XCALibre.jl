#!/usr/bin/env bash
# Gate between a phone-typed instruction and the queue. Runs on the LOGIN node
# before any sbatch. Non-zero exit = nothing is submitted.
# Keep it under VALIDATE_TIMEOUT (900 s) and off the compute nodes.

# Cron (and any non-interactive `bash -c`, which is how poller.sh invokes
# this) doesn't source /etc/profile.d, so Lmod's `module` function may not
# exist here even though it's fine in an interactive ssh shell. Bootstrap it
# before turning on -u below, since Lmod's own init scripts aren't
# nounset-clean.
if ! declare -F module >/dev/null 2>&1; then
  for f in /etc/profile.d/z00_lmod.sh /etc/profile.d/00-modulepath.sh; do
    [[ -r "$f" ]] && source "$f"
  done
fi

set -euo pipefail

if declare -F module >/dev/null 2>&1 && [[ -n "${CFD_MODULES:-}" ]]; then
  for m in $CFD_MODULES; do
    module load "$m" || { echo "FAIL: could not load module $m"; exit 1; }
  done
fi

export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-$HOME/.julia}"

echo "== instantiate =="
# Do this here, on the login node, where there IS network. Compute nodes run
# with JULIA_PKG_OFFLINE=true and cannot download anything.
julia --project=. -e '
  using Pkg
  Pkg.instantiate()
  Pkg.precompile()
' || { echo "FAIL: instantiate/precompile"; exit 1; }

echo "== load =="
julia --project=. -e 'using XCALibre' || { echo "FAIL: package does not load"; exit 1; }

echo "== load (offline) =="
# The actual compute-node job runs with JULIA_PKG_OFFLINE=true (no network
# there). Instantiating above can silently paper over a Manifest that's
# stale relative to the installed Julia version or registry state, since
# online mode is free to fix things up quietly. Re-check under the same
# offline constraint the real job faces, so that gap gets caught here
# instead of after a real sbatch submission.
JULIA_PKG_OFFLINE=true julia --project=. -e 'using XCALibre' \
  || { echo "FAIL: package does not load offline (JULIA_PKG_OFFLINE=true) — Manifest.toml is likely out of sync with this Julia version; try Pkg.resolve()"; exit 1; }

echo "== fast tests =="
# Swap for your real smoke test. Keep it short — this runs on every submission.
# Good candidates: an MMS convergence check on a tiny mesh, or a single-step
# solve on a 8x8 case. Bad candidate: the full test suite.
if [[ -f test/smoke.jl ]]; then
  julia --project=. test/smoke.jl || { echo "FAIL: smoke test"; exit 1; }
else
  echo "  (no test/smoke.jl — add one, this is your only real guard)"
fi

echo "validation passed"
