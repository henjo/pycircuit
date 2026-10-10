
## 2026-10-10 — 94. Stage R1d, step 1: one interface text per unit, the shared types once (machine code identical)

Each runtime unit `X.c` now has a header `csrc/X.h` -- its struct (cut from the `.c` byte for
byte) and its `hdl_fn` declaration (cut from the module's cdef).  The module's cdef IS that
header's text (`_cbuild.csrc('X.h')`) and `X.c` includes it, so cffi lays the struct out from
the very text the compiler compiles, and the compiler checks the definition against the
declaration cffi reads.  Before, each struct was written twice (C and cdef) and R1a's test
compared them; that comparison would now compare a text with itself, so the test asserts the
construction instead (the cdef is the header, the unit includes it once, no `.c` declares a
type, no header includes another).

The function types the units passed each other by hand-written copies -- `core_fn_t` and
`walk_fn_t` (3 copies each), `dgesv_t` (4), `dgetrf_t` and `dgemv_t` (2 each), and `dgetrs_t`
(ONE name for TWO signatures: 32-bit integers in newton.c and sens.c, 64-bit in radau_tc.c) --
are declared once in `csrc/pcx_abi.h`, under names that say what they are: `pass_fn_t` (core.c's
`hdl_fn_t` and pass.c's `hdl_batch_fn_t`, the same contract), `fused_fn_t`, `lim_fn_t`,
`dgesv64_t` / `dgetrs64_t` / `dgemv64_t` (numpy's ILP64 OpenBLAS), `dgetrf32_t` / `dgetrs32_t`
(SciPy's LP64), `klu_refactor_t` / `klu_solve_t`.  `_cbuild.csrc` expands `#include "x.h"` in
place, so a unit's key and its pinned pair cover every line compiled.

Proof: every unit's MACHINE CODE identical -- the `.text` and `.rodata` bytes and the exported
symbols of all ten units, compiled as the store compiles them, before and after (an A/A of the
instrument identical first); each new cdef declares what the old one did (normalised text and
cffi's layout); no unit gains a `-Wall -Wextra` warning.  The 12 pinned pairs that digest a
unit's text moved on the twin only (re-recorded; no reference moved).  `limit.c` and
`psp_limit.c` untouched (no model limiter rebuilt).  Full suite 5028 passed; the sanitizer's
fast tier clean (4722 passed).
