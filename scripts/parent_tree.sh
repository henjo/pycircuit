#!/usr/bin/env bash
# THE PARENT TREE for `benchmarks/step_machinery.py --check` (2026-10-04,
# testing for development, stage 3): a git worktree at HEAD -- before a
# commit, HEAD is the parent of the working tree -- or at the revision given.
#
#   scripts/parent_tree.sh [REV]          # default HEAD
#
# Where: $PYCIRCUIT_PARENT_TREE, else ~/.cache/pycircuit/wt_parent (the
# check's default).  Stale worktree records are pruned first.
set -eu
cd "$(dirname "$0")/.."
WT=${PYCIRCUIT_PARENT_TREE:-$HOME/.cache/pycircuit/wt_parent}
REV=$(git rev-parse "${1:-HEAD}")
git worktree prune
if [ -e "$WT/.git" ]; then
    git -C "$WT" checkout -q --detach "$REV"
else
    mkdir -p "$(dirname "$WT")"
    git worktree add -q --detach "$WT" "$REV"
fi
echo "parent tree $WT at $(git -C "$WT" log --oneline -1 | cut -c1-72)"
