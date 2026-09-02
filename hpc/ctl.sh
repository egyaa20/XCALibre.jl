#!/usr/bin/env bash
# Pause/resume/status the pipeline. Cron keeps ticking every minute either
# way — a touch-file just makes each tick a no-op. This avoids editing
# crontab from a script, which is fiddlier to get exactly right.
#
# Usage on Ada:      ~/cfd-pipeline/hpc/ctl.sh {pause|resume|status|kill}
# Usage from laptop: ssh ada "~/cfd-pipeline/hpc/ctl.sh pause"
set -uo pipefail

STATE_DIR="${STATE_DIR:-$HOME/.cfd-pipeline}"
FLAG="$STATE_DIR/PAUSED"
mkdir -p "$STATE_DIR"

case "${1:-}" in
  pause)
    date -u +%FT%TZ >"$FLAG"
    echo "paused. cron will keep ticking but every sweep is now a no-op."
    echo "already-running jobs on SLURM are NOT affected — this only stops"
    echo "new submissions and repo pulls."
    ;;
  resume)
    rm -f "$FLAG"
    echo "resumed. next cron tick (within 60s) will sweep normally."
    ;;
  status)
    if [[ -f "$FLAG" ]]; then
      echo "PAUSED since $(cat "$FLAG")"
    else
      echo "ACTIVE"
    fi
    hb="$STATE_DIR/heartbeat"
    if [[ -f "$hb" ]]; then
      age=$(( $(date +%s) - $(cat "$hb") ))
      echo "last heartbeat: ${age}s ago"
    else
      echo "no heartbeat yet — poller may never have run"
    fi
    crontab -l 2>/dev/null | grep -q xcalibre-pipeline \
      && echo "cron entry: present" || echo "cron entry: ABSENT"
    tmux has-session -t xcapoll 2>/dev/null \
      && echo "tmux session: present" || echo "tmux session: absent"
    ;;
  kill)
    crontab -l 2>/dev/null | grep -v xcalibre-pipeline | crontab - 2>/dev/null
    tmux kill-session -t xcapoll 2>/dev/null
    echo "cron entry removed, tmux session killed (if either existed)."
    echo "the pipeline will not resume on its own. re-run install_ada.sh to restart it."
    ;;
  *)
    echo "usage: $0 {pause|resume|status|kill}" >&2
    exit 2
    ;;
esac
