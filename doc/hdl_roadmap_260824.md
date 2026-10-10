
## 2026-10-10 — 96. Stage R1d, step 3: the runtime units built without the HDL kernels' prelude

Every object was built from `hdl._KERNEL_C` + its source: the prelude of the printed kernels
(`<math.h>`, the libm functions re-declared `__attribute__((const))`, `_sel`, `_npmax`,
`quartic_root` ...).  The runtime units used none of it knowingly, but one depended on it:
`sens.c` calls `isfinite` and never included `<math.h>` -- the prelude did.  Now
`_cbuild.load_unit(csrc, cdef)` builds a unit as the complete C file it is (`prelude=False`
through `source_key`, `_build`, `load_kernel`; the key hashes an empty prelude part, so no HDL
kernel's key moved), `sens.c` includes `<math.h>`, and the ten unit drivers no longer import the
HDL backend to build (`_cbuild.load_unit`, `_cbuild.CompileError`).

Proof: nine units' machine code identical (the objects the store built, hashed).  `radau_tc`'s
differs, by one thing: the prelude declared `hypot` `const`, and it is not -- it raises the
floating-point flags this unit tests.  Its four `hypot` calls sit between a `feclearexcept` and
the `fetestexcept` that judges them, each result used before that test, so data dependencies
already fix their place: the same four call sites either way, 26 instructions fewer without the
attribute; the values and the flags cannot move.  The 102 Radau tests (their C-against-Python
bit-identity among them) pass.  test_c_abi gains: every unit compiles ON ITS OWN under `-Wall
-Wextra -Werror` (alive: `sens.c` without its include is refused, implicit `isfinite`).  The two
tests that build a planted variant of a unit (core.c, pass.c) build it as the unit is built.
One pinned pair moved on its twin (sens.c's include).  Full suite 5050 passed; the sanitizer's
fast tier clean (4758 passed).
