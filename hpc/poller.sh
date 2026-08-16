#!/usr/bin/env bash
# XCALibre pipeline poller. Runs on the Ada login node.
#   poller.sh --once     one sweep (cron: * * * * *)
#   poller.sh --loop     sweep forever (tmux fallback when cron is unavailable)
#   poller.sh --status   print state
# Submits and reports only; never runs solver work on the login node.

# Cron's plain `bash -c` doesn't source /etc/profile.d, so neither Lmod's
# `module` function nor the `slurm` module (sbatch/squeue/sacct) are loaded
# here, unlike in an interactive shell. Bootstrap before -u, since Lmod's
# init scripts aren't nounset-clean.
if ! declare -F module >/dev/null 2>&1; then
  for f in /etc/profile.d/z00_lmod.sh /etc/profile.d/00-modulepath.sh; do
    [[ -r "$f" ]] && source "$f"
  done
fi
command -v sbatch >/dev/null 2>&1 || { declare -F module >/dev/null 2>&1 && module load slurm 2>/dev/null; }

set -uo pipefail

MODE="${1:---once}"
POLL_SECONDS="${POLL_SECONDS:-60}"
CONF="${CFD_CONF:-$HOME/cfd-pipeline/pipeline.conf}"
[[ -r "$CONF" ]] || { echo "missing conf: $CONF" >&2; exit 1; }
# shellcheck disable=SC1090
source "$CONF"

WHOAMI="${USER:-$(id -un)}"
STATE_DIR="${STATE_DIR:-$HOME/.cfd-pipeline}"
mkdir -p "$STATE_DIR/jobs" "$STATE_DIR/logs" "$STATE_DIR/publish"
PUBLISH_TMP_BRANCH="__publish_staging"   # scratch ref for the orphan publish
LOCK="$STATE_DIR/poller.lock"
HEARTBEAT="$STATE_DIR/heartbeat"
LOG="$STATE_DIR/logs/poller.log"
SUBMITS="$STATE_DIR/submits.log"

log()    { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >>"$LOG"; }
say()    { log "$*"; [[ -t 2 ]] && printf '%s\n' "$*" >&2; return 0; }
notify() { [[ -n "${NTFY_URL:-}" ]] && curl -fsS -m 10 -d "$1" "$NTFY_URL" >/dev/null 2>&1; return 0; }

# ------------------------------------------------------------------ git sync
ensure_clone() {
  [[ -d "$REPO_DIR/.git" ]] && return 0
  say "cloning $REPO_URL"
  git clone --quiet --branch "$QUEUE_BRANCH" "$REPO_URL" "$REPO_DIR"
}

sync_repo() {
  local remote local_sha
  remote=$(git -C "$REPO_DIR" ls-remote origin "refs/heads/$QUEUE_BRANCH" 2>/dev/null | awk '{print $1}')
  [[ -n "$remote" ]] || { say "WARN cannot reach origin"; return 1; }
  local_sha=$(git -C "$REPO_DIR" rev-parse HEAD 2>/dev/null || echo none)
  [[ "$remote" == "$local_sha" ]] && return 0
  say "new head $remote"
  git -C "$REPO_DIR" fetch --quiet --prune origin "$QUEUE_BRANCH" || return 1
  git -C "$REPO_DIR" reset --hard --quiet "origin/$QUEUE_BRANCH" || return 1
}

# ------------------------------------------------------------------ guards
active_jobs()      { squeue -h -u "$WHOAMI" -o '%i' 2>/dev/null | wc -l | tr -d ' '; }
submits_last_hour() {
  local now cutoff n=0 t
  now=$(date +%s); cutoff=$((now - 3600))
  [[ -f "$SUBMITS" ]] || { echo 0; return; }
  while read -r t; do [[ "$t" =~ ^[0-9]+$ ]] && ((t > cutoff)) && ((n++)); done <"$SUBMITS"
  echo "$n"
}

guards_ok() { # name -> 0 if we may submit now
  local act rate
  act=$(active_jobs); rate=$(submits_last_hour)
  if (( act >= MAX_ACTIVE_JOBS ));      then say "HOLD $1: $act active jobs";      return 1; fi
  if (( rate >= MAX_SUBMITS_PER_HOUR )); then say "HOLD $1: rate limit ($rate/h)"; return 1; fi
  return 0
}

validate_repo() {
  [[ -n "${VALIDATE_CMD:-}" ]] || return 0
  [[ -x "$REPO_DIR/${VALIDATE_CMD#./}" ]] || { say "WARN no $VALIDATE_CMD, skipping validation"; return 0; }
  ( cd "$REPO_DIR" && CFD_MODULES="$DEFAULT_MODULES" \
    timeout "${VALIDATE_TIMEOUT:-900}" bash -c "$VALIDATE_CMD" ) \
    >"$STATE_DIR/logs/validate.log" 2>&1
}

write_marker() { # name state jobid detail
  cat >"$STATE_DIR/jobs/$1.json" <<EOF
{"name":"$1","state":"$2","jobid":"$3","commit":"${COMMIT:-}","rundir":"${RUNDIR:-}",
 "updated":"$(date -u +%FT%TZ)","detail":"$4"}
EOF
}

snapshot() { # name -> sets RUNDIR to a frozen copy of the current commit
  RUNDIR="$SCRATCH_ROOT/${1}-$(date +%Y%m%d-%H%M%S)"
  mkdir -p "$RUNDIR/src" || return 1
  git -C "$REPO_DIR" archive HEAD | tar -x -C "$RUNDIR/src" || return 1
  # Manifest.toml is gitignored (deliberately, so it's not part of the
  # tracked tree), so `git archive` never includes it. Without it, the job
  # runs offline against no resolved manifest at all. Carry over the one
  # validate_repo() just instantiated/precompiled against in $REPO_DIR, so
  # what was validated is what actually runs.
  [[ -f "$REPO_DIR/Manifest.toml" ]] && cp "$REPO_DIR/Manifest.toml" "$RUNDIR/src/Manifest.toml"
}

# ------------------------------------------------------ mode: spec files
parse_spec() {
  local f=$1 line k v
  SOLVER_CMD="" CASE_DIR="."   # NAME is set by the caller
  NODES="$DEFAULT_NODES" NTASKS="$DEFAULT_NTASKS" CPUS_PER_TASK="$DEFAULT_CPUS_PER_TASK"
  TIME="$DEFAULT_TIME" PARTITION="$DEFAULT_PARTITION" ACCOUNT="$DEFAULT_ACCOUNT"
  MODULES="$DEFAULT_MODULES" GPUS="" MEM="" SBATCH_EXTRA=""
  while IFS= read -r line || [[ -n "$line" ]]; do
    line=${line%$'\r'}
    [[ "$line" =~ ^[[:space:]]*(#|$) ]] && continue
    k=${line%%=*}; v=${line#*=}
    k=${k//[[:space:]]/}
    v=${v#"${v%%[![:space:]]*}"}; v=${v%"${v##*[![:space:]]}"}
    case "$k" in
      SOLVER_CMD|CASE_DIR|NODES|NTASKS|CPUS_PER_TASK|TIME|PARTITION|ACCOUNT|MODULES|GPUS|MEM|SBATCH_EXTRA)
        printf -v "$k" '%s' "$v" ;;
      *) say "WARN $f: ignoring unknown key '$k'" ;;
    esac
  done <"$f"
  [[ -n "$SOLVER_CMD" ]]
}

render_template() {
  local tpl=$1 out=$2 line
  while IFS= read -r line || [[ -n "$line" ]]; do
    line=${line//@NAME@/$NAME};         line=${line//@NODES@/$NODES}
    line=${line//@NTASKS@/$NTASKS};     line=${line//@CPUS_PER_TASK@/$CPUS_PER_TASK}
    line=${line//@TIME@/$TIME};         line=${line//@PARTITION@/$PARTITION}
    line=${line//@MODULES@/$MODULES}
    line=${line//@SOLVER_CMD@/$SOLVER_CMD}
    line=${line//@CASE_DIR@/$CASE_DIR}; line=${line//@COMMIT@/$COMMIT}
    line=${line//@DEPOT@/${JULIA_DEPOT:-$HOME/.julia}}
    # RUNDIR last: SOLVER_CMD/CASE_DIR values (from the .run spec) may
    # themselves contain @RUNDIR@ and need it expanded too.
    line=${line//@RUNDIR@/$RUNDIR}
    if [[ "$line" == *@ACCOUNT_LINE@* ]]; then [[ -n "$ACCOUNT" ]] && line="#SBATCH --account=$ACCOUNT" || continue; fi
    if [[ "$line" == *@GPU_LINE@*   ]]; then [[ -n "$GPUS" ]] && line="#SBATCH --gpus=$GPUS" || continue; fi
    if [[ "$line" == *@MEM_LINE@*   ]]; then [[ -n "$MEM"  ]] && line="#SBATCH --mem=$MEM"   || continue; fi
    if [[ "$line" == *@EXTRA_LINE@* ]]; then [[ -n "$SBATCH_EXTRA" ]] && line="#SBATCH $SBATCH_EXTRA" || continue; fi
    printf '%s\n' "$line"
  done <"$tpl" >"$out"
}

submit_spec() {
  local name=$1 spec=$2
  COMMIT=$(git -C "$REPO_DIR" rev-parse --short HEAD)
  NAME="$name"
  parse_spec "$spec" || { say "REJECT $name: no SOLVER_CMD"; write_marker "$name" rejected "" "missing SOLVER_CMD"; return; }
  guards_ok "$name" || return
  if ! validate_repo; then
    say "REJECT $name: validation failed"; write_marker "$name" rejected "" "validation failed"
    cp "$STATE_DIR/logs/validate.log" "$STATE_DIR/publish/${name}.validate.log" 2>/dev/null
    notify "xcalibre: $name REJECTED (validation)"; return
  fi
  snapshot "$name" || { say "ERROR snapshot failed for $name"; return; }
  render_template "$SBATCH_TEMPLATE" "$RUNDIR/job.sbatch"

  if (( DRY_RUN )); then say "DRYRUN $name -> $RUNDIR/job.sbatch"; write_marker "$name" dryrun "" "DRY_RUN=1"; return; fi
  local jobid
  jobid=$(sbatch --parsable "$RUNDIR/job.sbatch" 2>"$STATE_DIR/logs/sbatch.err")
  if [[ -z "$jobid" ]]; then
    say "ERROR sbatch failed for $name: $(tail -3 "$STATE_DIR/logs/sbatch.err" | tr '\n' ' ')"
    write_marker "$name" failed "" "sbatch rejected"; notify "xcalibre: $name sbatch FAILED"; return
  fi
  date +%s >>"$SUBMITS"
  write_marker "$name" submitted "$jobid" "commit $COMMIT"
  say "SUBMIT $name job=$jobid commit=$COMMIT dir=$RUNDIR"
  notify "xcalibre: $name submitted (job $jobid, $COMMIT)"
}

process_specs() {
  shopt -s nullglob
  local spec name
  for spec in "$REPO_DIR/$RUNS_DIR"/*.run; do
    name=$(basename "$spec" .run)
    [[ -e "$STATE_DIR/jobs/$name.json" ]] && continue
    submit_spec "$name" "$spec"
  done
}

# ------------------------------------------------------ mode: repo hook
# Your original idea: a script in the repo decides everything and calls sbatch.
# The poller still owns the guards, the snapshot, and the reporting.
process_hook() {
  local sha key msg report line jname jid
  sha=$(git -C "$REPO_DIR" rev-parse HEAD)
  COMMIT="${sha:0:8}"
  key="c-$COMMIT"
  [[ -e "$STATE_DIR/jobs/$key.json" ]] && return
  msg=$(git -C "$REPO_DIR" log -1 --pretty='%s%n%b')
  if [[ -n "${RUN_TOKEN:-}" && "$msg" != *"$RUN_TOKEN"* ]]; then return; fi
  [[ -x "$REPO_DIR/${HOOK_CMD#./}" ]] || { say "WARN no $HOOK_CMD in repo"; return; }
  guards_ok "$key" || return
  if ! validate_repo; then
    say "REJECT $key: validation failed"; write_marker "$key" rejected "" "validation failed"
    cp "$STATE_DIR/logs/validate.log" "$STATE_DIR/publish/${key}.validate.log" 2>/dev/null
    notify "xcalibre: $COMMIT REJECTED (validation)"; return
  fi
  snapshot "$key" || { say "ERROR snapshot failed"; return; }

  report="$RUNDIR/.jobids"; : >"$report"
  say "HOOK running $HOOK_CMD for $COMMIT"
  ( cd "$RUNDIR/src" && \
    CFD_RUNDIR="$RUNDIR" CFD_COMMIT="$COMMIT" CFD_JOB_REPORT="$report" \
    CFD_ACCOUNT="$DEFAULT_ACCOUNT" CFD_PARTITION="$DEFAULT_PARTITION" \
    CFD_MODULES="$DEFAULT_MODULES" CFD_DEPOT="${JULIA_DEPOT:-$HOME/.julia}" \
    CFD_DRY_RUN="$DRY_RUN" \
    timeout "${HOOK_TIMEOUT:-300}" bash "./${HOOK_CMD#./}" ) \
    >"$RUNDIR/hook.log" 2>&1
  local rc=$?
  tail -n 100 "$RUNDIR/hook.log" >"$STATE_DIR/publish/${key}.hook.log" 2>/dev/null
  write_marker "$key" "hook-rc-$rc" "" "$(head -c 120 <<<"$msg" | tr '\n' ' ')"
  (( rc != 0 )) && { say "HOOK failed rc=$rc"; notify "xcalibre: hook failed for $COMMIT"; return; }

  # hook reports submissions as "name jobid" lines
  while read -r jname jid; do
    [[ "$jid" =~ ^[0-9]+$ ]] || continue
    date +%s >>"$SUBMITS"
    NAME="$jname" write_marker "$COMMIT-$jname" submitted "$jid" "via hook"
    say "SUBMIT $jname job=$jid commit=$COMMIT"
    notify "xcalibre: $jname submitted (job $jid, $COMMIT)"
  done <"$report"
}

# ------------------------------------------------------------------ tracking
job_state() {
  local id=$1 s=""
  command -v sacct  >/dev/null && s=$(sacct -j "$id" -n -X -P -o State 2>/dev/null | head -1)
  [[ -z "$s" ]] && command -v squeue >/dev/null && s=$(squeue -h -j "$id" -o '%T' 2>/dev/null | head -1)
  echo "${s%% *}"
}

track_jobs() {
  shopt -s nullglob
  local f name id state prev rundir
  for f in "$STATE_DIR"/jobs/*.json; do
    prev=$(sed -n 's/.*"state":"\([^"]*\)".*/\1/p' "$f" | head -1)
    case "$prev" in submitted|PENDING|RUNNING|CONFIGURING|COMPLETING) ;; *) continue ;; esac
    name=$(basename "$f" .json)
    id=$(sed -n 's/.*"jobid":"\([^"]*\)".*/\1/p' "$f" | head -1)
    rundir=$(sed -n 's/.*"rundir":"\([^"]*\)".*/\1/p' "$f" | head -1)
    COMMIT=$(sed -n 's/.*"commit":"\([^"]*\)".*/\1/p' "$f" | head -1)
    [[ -n "$id" ]] || continue
    state=$(job_state "$id"); [[ -z "$state" ]] && state=UNKNOWN
    case "$state" in
      PENDING|RUNNING|CONFIGURING|COMPLETING)
        [[ "$state" != "$prev" ]] && { RUNDIR=$rundir; write_marker "$name" "$state" "$id" ""; } ;;
      COMPLETED|FAILED|CANCELLED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|PREEMPTED|DEADLINE|BOOT_FAIL)
        RUNDIR=$rundir; write_marker "$name" "$state" "$id" "terminal"
        collect_artifacts "$name" "$rundir"
        say "DONE $name job=$id state=$state"
        notify "xcalibre: $name $state (job $id)" ;;
    esac
  done
}

collect_artifacts() {
  local name=$1 rundir=$2 f
  [[ -d "$rundir" ]] || return 0
  for f in "$rundir"/slurm-*.out "$rundir"/*.log; do
    [[ -f "$f" ]] || continue
    tail -n "${LOG_TAIL_LINES:-300}" "$f" >"$STATE_DIR/publish/${name}.$(basename "$f")"
  done
  if [[ -d "$rundir/summary" ]]; then
    mkdir -p "$STATE_DIR/publish/$name"
    cp -r "$rundir/summary/." "$STATE_DIR/publish/$name/" 2>/dev/null || true
  fi
}

# ------------------------------------------------------------------ publish
ensure_status_repo() {
  [[ -d "$STATUS_DIR/.git" ]] && return 0
  if ! git clone --quiet --branch "$STATUS_BRANCH" --single-branch "$REPO_URL" "$STATUS_DIR" 2>/dev/null; then
    git init -q "$STATUS_DIR"
    git -C "$STATUS_DIR" remote add origin "$REPO_URL"
    git -C "$STATUS_DIR" checkout -q --orphan "$STATUS_BRANCH"
  fi
  git -C "$STATUS_DIR" config user.name  "ada-poller"
  git -C "$STATUS_DIR" config user.email "ada-poller@localhost"
}

build_index() {
  local f name state id upd hb age
  hb=$(cat "$HEARTBEAT" 2>/dev/null || echo 0); age=$(( $(date +%s) - hb ))
  {
    echo "# XCALibre pipeline status"
    echo
    echo "- updated: $(date -u +%FT%TZ)"
    echo "- host: $(hostname -s)   heartbeat age: ${age}s"
    echo "- slurm jobs for $WHOAMI: $(active_jobs)   submits last hour: $(submits_last_hour)"
    echo "- trigger mode: ${TRIGGER_MODE:-spec}   dry run: ${DRY_RUN:-0}"
    [[ -f "$STATE_DIR/PAUSED" ]] && echo "- **PAUSED** since $(cat "$STATE_DIR/PAUSED" 2>/dev/null)"
    echo
    echo "| run | state | jobid | updated |"
    echo "|---|---|---|---|"
    for f in "$STATE_DIR"/jobs/*.json; do
      [[ -f "$f" ]] || continue
      name=$(basename "$f" .json)
      state=$(sed -n 's/.*"state":"\([^"]*\)".*/\1/p' "$f" | head -1)
      id=$(sed -n 's/.*"jobid":"\([^"]*\)".*/\1/p' "$f" | head -1)
      upd=$(sed -n 's/.*"updated":"\([^"]*\)".*/\1/p' "$f" | head -1)
      echo "| $name | $state | $id | $upd |"
    done
  } >"$STATUS_DIR/INDEX.md"
}

publish() {
  ensure_status_repo
  rm -rf "$STATUS_DIR/jobs" "$STATUS_DIR/results"
  mkdir -p "$STATUS_DIR/jobs" "$STATUS_DIR/results"
  cp "$STATE_DIR"/jobs/*.json "$STATUS_DIR/jobs/" 2>/dev/null || true
  ( cd "$STATE_DIR/publish" 2>/dev/null &&
    find . -type f -size -"$(( ${MAX_PUBLISH_BYTES:-10485760} / 1024 ))"k \
      -exec cp --parents {} "$STATUS_DIR/results/" \; ) 2>/dev/null || true
  tail -n 200 "$LOG" >"$STATUS_DIR/poller.log"
  build_index
  git -C "$STATUS_DIR" add -A
  git -C "$STATUS_DIR" diff --cached --quiet && return 0
  # Publish as a single root commit. Interactive scenes run to several MB
  # each; appending them would grow every clone of this branch without
  # bound. A force-push only moves the ref -- the old commits stay
  # reachable from it -- so the history has to be dropped deliberately.
  git -C "$STATUS_DIR" branch -q -D "$PUBLISH_TMP_BRANCH" 2>/dev/null || true
  git -C "$STATUS_DIR" checkout -q --orphan "$PUBLISH_TMP_BRANCH" \
    || { say "WARN orphan checkout failed, skipping publish"; return 0; }
  git -C "$STATUS_DIR" add -A
  git -C "$STATUS_DIR" commit -q -m "status $(date -u +%FT%TZ)"
  git -C "$STATUS_DIR" branch -q -M "$STATUS_BRANCH"
  git -C "$STATUS_DIR" push -q --force origin "$STATUS_BRANCH" 2>>"$LOG" || say "WARN status push failed"
  git -C "$STATUS_DIR" gc --auto --quiet 2>/dev/null || true
}

# ------------------------------------------------------------------ sweep
sweep() {
  date +%s >"$HEARTBEAT"
  if [[ -f "$STATE_DIR/PAUSED" ]]; then
    say "paused (touch-file $STATE_DIR/PAUSED present), skipping sweep"
    ensure_clone || return 1   # keep status branch fresh even while paused
    publish
    return 0
  fi
  ensure_clone || return 1
  sync_repo
  case "${TRIGGER_MODE:-spec}" in
    spec) process_specs ;;
    hook) process_hook  ;;
    both) process_specs; process_hook ;;
    *)    say "ERROR unknown TRIGGER_MODE=${TRIGGER_MODE}" ;;
  esac
  track_jobs
  publish
}

case "$MODE" in
  --status)
    ensure_status_repo 2>/dev/null; build_index 2>/dev/null
    cat "$STATUS_DIR/INDEX.md" 2>/dev/null || echo "no status yet"; exit 0 ;;
  --once)
    exec 9>"$LOCK"; flock -n 9 || { echo "sweep already running"; exit 0; }; sweep ;;
  --loop)
    exec 9>"$LOCK"; flock -n 9 || { echo "poller already running"; exit 0; }
    trap 'say "poller stopping"; exit 0' INT TERM
    say "poller loop started (pid $$)"
    while true; do sweep; sleep "$POLL_SECONDS"; done ;;
  *) echo "usage: $0 [--once|--loop|--status]" >&2; exit 2 ;;
esac
