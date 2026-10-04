#!/usr/bin/env bash
# VALGRIND WITHOUT ROOT, for instruction counts (2026-10-04, testing for
# development, stage 6): the distribution's package unpacked into a
# directory of your own -- nothing installed -- and the helper that switches
# cachegrind's counting on and off around a call (`CACHEGRIND_START/
# STOP_INSTRUMENTATION`, valgrind >= 3.22).  `benchmarks/step_machinery.py
# --count` uses it when the kernel's counters are closed
# (`perf_event_paranoid` > 2).
#
#   scripts/get_valgrind.sh [DIR]      # default ~/.local/opt/valgrind
#
# Needs apt-get (Debian/Ubuntu: `apt-get download` runs unprivileged) and gcc.
set -eu
DIR=${1:-$HOME/.local/opt/valgrind}
mkdir -p "$DIR"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
(cd "$TMP" && apt-get download valgrind >/dev/null)
dpkg -x "$TMP"/valgrind_*.deb "$DIR"
cat > "$TMP/count.c" <<'EOF'
#include <valgrind/cachegrind.h>
void pyc_count_start(void) { CACHEGRIND_START_INSTRUMENTATION; }
void pyc_count_stop(void) { CACHEGRIND_STOP_INSTRUMENTATION; }
EOF
gcc -O2 -shared -fPIC -I "$DIR/usr/include" -o "$DIR/pycircuit-count.so" "$TMP/count.c"
"$DIR/usr/bin/valgrind.bin" --version >/dev/null 2>&1 || \
    env VALGRIND_LIB="$DIR/usr/libexec/valgrind" "$DIR/usr/bin/valgrind.bin" --version
echo "valgrind $(env VALGRIND_LIB="$DIR/usr/libexec/valgrind" "$DIR/usr/bin/valgrind.bin" --version) in $DIR"
