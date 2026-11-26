#!/usr/bin/env bash
# CHOLMOD WITHOUT OPENMP (2026-11-26, doc/hdl_roadmap_260824.md 362): the
# system's SuiteSparse version of CHOLMOD built from its source with OpenMP
# off, into pycircuit's cache, where `_cholklu` opens it before the system's
# (THE LIBRARY in its note).  The system's CHOLMOD opens an OpenMP region for
# each update of a supernode, whatever the thread count; without OpenMP a
# factorization of ibmpg6t's matrix takes 0.40x the time, the same bytes.
# Only CHOLMOD is built: it is linked against the system's AMD, CAMD, COLAMD,
# CCOLAMD, SuiteSparse_config and BLAS/LAPACK.  Rerun it after the system's
# SuiteSparse is upgraded (a version other than the header's is not opened:
# `_cholklu` warns and takes the system's).
#
#   scripts/build_cholmod.sh [DIR]   # default ${XDG_CACHE_HOME:-~/.cache}/pycircuit/cholmod
#
# Needs curl, cmake (>= 3.22), a C compiler and the system's SuiteSparse
# headers and libraries (Debian/Ubuntu: libsuitesparse-dev); fetches the
# SuiteSparse release's source from GitHub (~95 MB).
set -eu
DIR=${1:-${XDG_CACHE_HOME:-$HOME/.cache}/pycircuit/cholmod}
H=
for d in /usr/include/suitesparse /usr/local/include/suitesparse /usr/include; do
    [ -f "$d/SuiteSparse_config.h" ] && { H=$d; break; }
done
[ -n "$H" ] || { echo "no SuiteSparse_config.h: install the system's SuiteSparse headers (libsuitesparse-dev)" >&2; exit 1; }
v() { sed -n "s/^#define SUITESPARSE_$1_VERSION *\([0-9]*\).*/\1/p" "$H/SuiteSparse_config.h"; }
V=$(v MAIN).$(v SUB).$(v SUBSUB)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
curl -fsSL -o "$TMP/ss.tar.gz" "https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/refs/tags/v$V.tar.gz"
echo "SuiteSparse $V: sha256 $(sha256sum "$TMP/ss.tar.gz" | cut -d' ' -f1)"
## (CHOLMOD's cmake reads SuiteSparse_config's modules beside it)
tar -C "$TMP" -xzf "$TMP/ss.tar.gz" "SuiteSparse-$V/CHOLMOD" "SuiteSparse-$V/SuiteSparse_config"
cmake -S "$TMP/SuiteSparse-$V/CHOLMOD" -B "$TMP/build" -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$DIR" -DCMAKE_INSTALL_LIBDIR=lib \
    -DCHOLMOD_USE_OPENMP=OFF -DSUITESPARSE_USE_OPENMP=OFF \
    -DSUITESPARSE_USE_CUDA=OFF -DCHOLMOD_USE_CUDA=OFF \
    -DBUILD_STATIC_LIBS=OFF -DSUITESPARSE_DEMOS=OFF -DBLA_VENDOR=Generic >"$TMP/cmake.log" 2>&1 \
    || { cat "$TMP/cmake.log" >&2; exit 1; }
cmake --build "$TMP/build" -j "$(nproc)" >"$TMP/build.log" 2>&1 || { tail -40 "$TMP/build.log" >&2; exit 1; }
cmake --install "$TMP/build" >/dev/null
if nm -D --undefined-only "$DIR/lib/libcholmod.so.5" | grep -q GOMP; then
    echo "$DIR/lib/libcholmod.so.5 calls OpenMP: not the build wanted" >&2
    exit 1
fi
echo "CHOLMOD of SuiteSparse $V without OpenMP: $DIR/lib/libcholmod.so.5"
