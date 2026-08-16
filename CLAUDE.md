# Repo protocol — read before acting

This repo is driven remotely from a phone, and pushes to `queue` launch SLURM
jobs on Ada within 60 seconds. A careless push costs real allocation hours.

## Branches
- `main`  — reviewed work. Never force-push.
- `dev`   — normal development. Push freely.
- `queue` — **watched by Ada.** Merging here submits jobs.
- `pipeline-status` — written by Ada only, force-pushed. Read it, never write it.

## Launching a run (TRIGGER_MODE=spec)
1. Work on `dev`; make sure `./validate.sh` passes.
2. Add `runs/<descriptive-name>.run`. The filename is the dedup key — a name
   that already ran is ignored forever. Never reuse one; encode the varied
   parameter: `tgv-re1600-weno5-nx256.run`.
3. Merge `dev` into `queue` and push. That is the trigger.
4. Max two new `.run` files per push unless asked otherwise.

## Launching a run (TRIGGER_MODE=hook)
Edit `pipeline/on_update.sh`, then push to `queue` with `[run]` in the commit
message. Without that token nothing fires.

## Checking a run
```
git fetch origin pipeline-status
git show origin/pipeline-status:INDEX.md
```
Artifacts under `results/<run>/`. Raw fields are not in git — they are on Ada
scratch at the `rundir` recorded in `jobs/<run>.json`.

If `heartbeat age` in INDEX.md exceeds a few minutes, the poller is dead and
nothing is being submitted. Say so rather than assuming the queue is just idle.

## Hard rules
- Never `git push --force` except to `pipeline-status` (which is Ada's job).
- Never commit files > 1 MB to `main`/`dev`/`queue`: no meshes, restarts, or
  field output. Ada may publish artifacts up to 10 MB to `pipeline-status`
  (interactive scenes are the reason); that branch is rewritten as a single
  orphan commit per publish, so those blobs never accumulate in history.
- Never edit `hpc/`, `pipeline.conf`, or `validate.sh` unless explicitly asked.
  Those enforce the guards that stop a bad instruction from burning the
  allocation. Weakening them is never an implied part of "make the run work".
- "Run X" means add a spec and push, not execute a solver locally.
- If a run failed, read the log in `pipeline-status` and propose a fix. Do not
  resubmit an identical configuration.

## Scientific conventions
- New schemes/models need a verification case (MMS or an established benchmark)
  before they reach `queue`.
- Cite the commit SHA on any figure or note. The pipeline freezes the tree at
  submit time, so results map to exactly one commit.
