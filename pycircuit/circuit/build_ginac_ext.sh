#!/bin/bash
# Build the optional GiNaC linear-solve extension (_ginac_ext).
#
# Requires: GiNaC (discoverable via `pkg-config ginac`), pybind11, python dev
# headers, and g++.  If GiNaC is only in a user prefix, point pkg-config at it:
#     export PKG_CONFIG_PATH=$HOME/.local/ginac/lib/pkgconfig
#
# Usage:
#     PYTHON=/path/to/venv/bin/python ./build_ginac_ext.sh
# (defaults to `python3`).  The resulting _ginac_ext.so is gitignored; without
# it, GinacToolkit / the ginac_toolkit singleton are simply unavailable.
set -e
cd "$(dirname "$0")"

PYTHON=${PYTHON:-python3}
PYINC=$("$PYTHON" -c "import sysconfig; print(sysconfig.get_path('include'))")
PBINC=$("$PYTHON" -c "import pybind11; print(pybind11.get_include())")

if ! pkg-config --exists ginac; then
    echo "error: GiNaC not found by pkg-config (set PKG_CONFIG_PATH?)." >&2
    exit 1
fi

g++ -O2 -std=c++17 -shared -fPIC \
    -I"$PYINC" -I"$PBINC" \
    _ginac_ext.cpp -o _ginac_ext.so \
    $(pkg-config --cflags --libs ginac)

echo "built $(pwd)/_ginac_ext.so"
