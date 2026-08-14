#!/usr/bin/env bash
# Gate between a phone-typed instruction and the queue. Runs on the LOGIN node
# before any sbatch. Non-zero exit = nothing is submitted.
# Keep it under VALIDATE_TIMEOUT (900 s) and off the compute nodes.
set -euo pipefail

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
