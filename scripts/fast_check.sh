#!/usr/bin/env bash
# THE SWITCHES-OFF CHECK (2026-10-04, testing for development, stage 1).
#
# The fast tier, recorded twice -- every fast path on, then all seventeen off --
# and compared bit for bit.  Each fast path's contract is the same bytes as
# the Python path it replaces (the evaluate core, the pass batches, the limit
# walk, the limiter kernel, the stamp plan, the fused passes, the zero-source
# skip, the predictor's multistep fast path, the Newton solve in C, the
# error test in C, PSP's limiter in C, the converged point without its
# unread G, the fused device kernels, the stamped readiness, the source
# pass as a plan, the predictor's fit in C, the core's passes for the
# stage methods); this checks every one of them on every circuit of the
# fast tier at once.
# ~6 min.
#
#   scripts/fast_check.sh [OUTDIR] [-- extra pytest arguments]
#
# A test that asserts a fast path SERVED fails in the off run by design: the
# comparison covers the tests that passed on both sides and names the rest.
# Exit status: the comparison's (0 = every compared call identical).  The
# fast-path counts (family `paths`) are not recorded here: with the paths off
# they differ by design.
set -u
cd "$(dirname "$0")/.." || exit 2
OUT=${1:-$(mktemp -d -t fast_check.XXXXXX)}
[ $# -gt 0 ] && shift
[ "${1:-}" = "--" ] && shift
EXTRA=("$@")
PY=.venv/bin/python
OFF=(PYCIRCUIT_TRAN_CORE=0 PYCIRCUIT_HDL_BATCH=0 PYCIRCUIT_HDL_LIMIT_WALK=0
     PYCIRCUIT_HDL_CLIMIT=0 PYCIRCUIT_STAMP_PLAN=0 PYCIRCUIT_HDL_FUSE=0
     PYCIRCUIT_HDL_ZERO_U=0 PYCIRCUIT_PRED_FAST=0
     PYCIRCUIT_NEWTON_C=0 PYCIRCUIT_LTE_C=0 PYCIRCUIT_PSP_LIMIT_C=0
     PYCIRCUIT_SKIP_UNREAD_J=0 PYCIRCUIT_HDL_CFUSE=0 PYCIRCUIT_WATCH=0
     PYCIRCUIT_SOURCE_PLAN=0 PYCIRCUIT_CORE_PASSES=0
     PYCIRCUIT_RADAU_C=0)

run() {   # run LABEL [VAR=value ...]
    local label=$1; shift
    mkdir -p "$OUT/$label"
    env "$@" PYTHONPATH=benchmarks/tranrec TRANREC_OUT="$OUT/$label" TRANREC_FAMILIES=transient,pss,pac \
        PYCIRCUIT_LEAKS_REPORT="$OUT/$label/leaks" \
        "$PY" -m pytest pycircuit -q -p no:cacheprovider -p tran_recorder --tier fast \
        "${EXTRA[@]}" > "$OUT/$label.log" 2>&1
    echo "fast_check: $label: $(grep -E '[0-9]+ (passed|failed)' "$OUT/$label.log" | tail -1)"
}

echo "fast_check: recordings in $OUT"
run on
run off "${OFF[@]}"
"$PY" benchmarks/tranrec/compare.py "$OUT/on" "$OUT/off" --passed-only
status=$?
## the off run's failures, by file: tests asserting that a fast path served
## (or that a fault planted in one is seen) -- anything else is a finding
echo "fast_check: off-run failures by file (expected: tests of the fast paths themselves):"
grep -E '^FAILED ' "$OUT/off.log" | sed -e 's/^FAILED //' -e 's/::.*//' | sort | uniq -c
## a failure with every fast path ON is the suite failing: not a pass
if grep -qE '^FAILED |^ERROR ' "$OUT/on.log"; then
    echo "fast_check: the run with the fast paths ON has failures (see $OUT/on.log)"
    status=1
fi
[ $status -eq 0 ] && echo "fast_check: every fast path gave the bytes of its Python path" \
                  || echo "fast_check: FAILED -- differences or failures above"
exit $status
