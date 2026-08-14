#!/usr/bin/env bash
# Run once on Ada:  bash ~/cfd-pipeline/hpc/install_ada.sh
# Decides cron vs tmux empirically instead of guessing.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POLLER="$HERE/poller.sh"
CONF="$HOME/cfd-pipeline/pipeline.conf"
STATE="$HOME/.cfd-pipeline"
MARK="# xcalibre-pipeline"

echo "== preflight =="
for c in git sbatch squeue flock; do
  command -v "$c" >/dev/null && echo "  ok   $c" || echo "  MISSING $c"
done
[[ -r "$CONF" ]] && echo "  ok   pipeline.conf" || { echo "  MISSING $CONF"; exit 1; }
chmod +x "$HERE"/*.sh 2>/dev/null
# shellcheck disable=SC1090
source "$CONF"
git ls-remote "$REPO_URL" HEAD >/dev/null 2>&1 \
  && echo "  ok   github reachable" || { echo "  FAIL cannot reach $REPO_URL"; exit 1; }
mkdir -p "$STATE/jobs" "$STATE/logs" "$STATE/publish" "$SCRATCH_ROOT"

echo "== first sweep =="
"$POLLER" --once || true
[[ -f "$STATE/heartbeat" ]] || { echo "  FAIL no heartbeat, check $STATE/logs/poller.log"; exit 1; }
echo "  ok   heartbeat written"

echo "== choosing a supervisor =="
hb_age() { echo $(( $(date +%s) - $(cat "$STATE/heartbeat" 2>/dev/null || echo 0) )); }

use_tmux() {
  command -v tmux >/dev/null || { echo "  FAIL no cron and no tmux"; exit 1; }
  tmux has-session -t xcapoll 2>/dev/null && tmux kill-session -t xcapoll
  tmux new-session -d -s xcapoll "exec $POLLER --loop"
  sleep 3
  tmux has-session -t xcapoll 2>/dev/null \
    && echo "  using TMUX session 'xcapoll' (will NOT survive a login-node reboot)" \
    || { echo "  FAIL tmux session did not start"; exit 1; }
}

if command -v crontab >/dev/null && \
   { crontab -l 2>/dev/null | grep -qF "$MARK" || \
     { crontab -l 2>/dev/null; echo "* * * * * $POLLER --once >/dev/null 2>&1 $MARK"; } | crontab - ; }; then
  echo "  cron entry installed, waiting 80 s to see if crond actually runs it..."
  : >"$STATE/heartbeat"; echo 0 >"$STATE/heartbeat"
  sleep 80
  if (( $(hb_age) < 75 )); then
    echo "  using CRON (self-heals across reboots)"
  else
    echo "  crond did not fire; removing cron entry and falling back"
    crontab -l 2>/dev/null | grep -vF "$MARK" | crontab - 2>/dev/null
    use_tmux
  fi
else
  echo "  crontab unavailable"
  use_tmux
fi

echo
echo "== state =="
"$POLLER" --status
echo
echo "Watch:      tail -f $STATE/logs/poller.log"
echo "Status:     $POLLER --status"
echo "Stop cron:  crontab -l | grep -v '$MARK' | crontab -"
echo "Stop tmux:  tmux kill-session -t xcapoll"
