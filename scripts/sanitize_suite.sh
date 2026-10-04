#!/usr/bin/env bash
# THE SANITIZED SUITE (2026-10-04, testing for development, stage 4).
#
# The suite with every C object the package builds compiled under
# AddressSanitizer and UndefinedBehaviorSanitizer (`PYCIRCUIT_C_SANITIZE=1`,
# `_hdl_cbackend.SANITIZE`): an out-of-bounds access or undefined behaviour
# in a kernel, the pass driver, the walk or the evaluate core stops its
# worker and is reported with the file and line.
#
#   scripts/sanitize_suite.sh [--deep] [OUTDIR] [-- extra pytest arguments]
#
# By default the FAST TIER (`--tier fast`, which reaches every C path) with
# Python's own allocator for its small objects: the arrays the C reads and
# writes are numpy's, allocated through malloc, so ASan guards them.  Before
# any commit that touches C.  `--deep`: the whole suite, and every Python
# allocation through ASan's allocator too (`PYTHONMALLOC=malloc`: cffi's
# small buffers get redzones as well) -- once per speed round; it takes
# hours (measured 2026-10-04: ~27 of the slowest tests in 17 min at -n 4).
#
# Exit status: nonzero when a sanitizer reported, a test failed, or a class
# that builds C ran numpy instead (a sanitized object that would not load).
set -u
cd "$(dirname "$0")/.." || exit 2
DEEP=0
if [ "${1:-}" = "--deep" ]; then DEEP=1; shift; fi
OUT=${1:-$(mktemp -d -t sanitize.XXXXXX)}
[ $# -gt 0 ] && shift
[ "${1:-}" = "--" ] && shift
PY=.venv/bin/python
ASAN=$(gcc -print-file-name=libasan.so)
case "$ASAN" in /*) ;; *) echo "sanitize_suite: no gcc AddressSanitizer runtime"; exit 2;; esac
## (libstdc++ beside it: Python is not C++, and without libstdc++ loaded
## when ASan starts, its `__cxa_throw` interceptor has nothing to call --
## jaxlib throws at import and ASan aborts on its own CHECK)
STDCXX=$(gcc -print-file-name=libstdc++.so)
mkdir -p "$OUT/reports" "$OUT/state"
echo "sanitize_suite: output in $OUT"
## (JAX on the CPU only: its CUDA plugin's cuDNN aborts under the
## sanitizer runtime, and no JAX path runs the package's C)
## the tests a sanitized build changes by design, deselected by name:
## the objdump count of `pow` calls (instrumented code calls differently)
## and the wall-time guards (ASan's allocator and checks are not the
## product's speed)
DESELECT=(
    --deselect 'pycircuit/circuit/tests/test_hdl_cbackend.py::test_const_merges_the_repeated_calls_and_keeps_the_bytes'
    --deselect 'pycircuit/circuit/tests/test_perf_guards.py'
    ## (the two JAX tests that held on the GPU backend only -- here JAX runs
    ## on the CPU -- were made backend-independent on 2026-10-04: the traced
    ## blocks to eight ulps, the cold start an expected failure off the GPU)
)
if [ $DEEP -eq 1 ]; then
    SCOPE=(--tier all); MALLOC=(PYTHONMALLOC=malloc)
else
    SCOPE=(--tier fast); MALLOC=()
fi
echo "sanitize_suite: ${SCOPE[*]}${MALLOC:+, ${MALLOC[*]}}"
env "${MALLOC[@]}" LD_PRELOAD="$ASAN $STDCXX" \
    ASAN_OPTIONS="detect_leaks=0:halt_on_error=1:log_path=$OUT/reports/asan" \
    UBSAN_OPTIONS="print_stacktrace=1:halt_on_error=1:log_path=$OUT/reports/ubsan" \
    PYCIRCUIT_C_SANITIZE=1 PYCIRCUIT_TEST_TIMINGS=0 \
    JAX_PLATFORMS=cpu \
    PYCIRCUIT_STATE_DUMP="$OUT/state" PYCIRCUIT_LEAKS_REPORT="$OUT/leaks" \
    "$PY" -m pytest pycircuit -q -s -p no:cacheprovider -n 4 --max-worker-restart=0 \
    "${SCOPE[@]}" "${DESELECT[@]}" "$@" > "$OUT/suite.log" 2>&1
status=$?
echo "sanitize_suite: $(grep -E '[0-9]+ (passed|failed)' "$OUT/suite.log" | tail -1)"
## (UBSan, combined with ASan, writes its report to stderr whatever its
## log_path says, and a test's stderr is pytest's capture file, lost with
## the aborted worker: the run is uncaptured, `-s`, so the report reaches
## the log; ASan's own reports go to the files above as well)
if grep -qE 'runtime error:|ERROR: AddressSanitizer' "$OUT/suite.log"; then
    echo "sanitize_suite: SANITIZER REPORTS IN THE LOG:"
    grep -E -A6 'runtime error:|ERROR: AddressSanitizer' "$OUT/suite.log" | grep -E 'runtime error|ERROR|SUMMARY| #[0-3] ' | head -20
    status=1
fi
if ls "$OUT"/reports/* >/dev/null 2>&1; then
    echo "sanitize_suite: SANITIZER REPORTS:"
    for f in "$OUT"/reports/*; do
        echo "--- $f"; grep -E 'ERROR|runtime error|SUMMARY| #[0-9] ' "$f" | head -20
    done
    status=1
fi
## every class that resolves to C did (a sanitized object that failed to
## load would leave it on numpy, and its tests passing on numpy)
"$PY" - "$OUT/state" <<'EOF' || status=1
import json, os, sys
d = sys.argv[1]
bad = set()
for f in sorted(os.listdir(d)):
    for cls, st in json.load(open(os.path.join(d, f)))['state'].items():
        for s in (str(st.get('status', '')), str(st.get('limit_status', ''))):
            if 'compile failed' in s or 'unloadable' in s:
                bad.add((cls, s))
for cls, s in sorted(bad):
    print(f'sanitize_suite: {cls} ran numpy: {s}')
sys.exit(1 if bad else 0)
EOF
[ $status -eq 0 ] && echo "sanitize_suite: clean" || echo "sanitize_suite: FAILED (above; $OUT/suite.log)"
exit $status
