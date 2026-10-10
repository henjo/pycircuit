
## 2026-10-10 — 97. Stage R1d, step 4: `#line` -- a report names the repository file and line (R1 done)

`_cbuild.csrc` marks each header it expands: `#line 1 "pycircuit/circuit/csrc/x.h"` before it,
`#line N "pycircuit/circuit/csrc/unit.c"` after it (the path from the source tree's root, the
same wherever the package lives -- it is part of the key).  A unit is compiled with nothing before
it, so its own lines need no other marker; a file including no header (limit.c, psp_limit.c) is
its text exactly -- no model limiter's key moved.  A failed build or a sanitizer report (the
sanitized builds carry `-g`) named a temporary `.c` that `_build` had already deleted; it now
names the file and line in the tree.

The plan's other half -- the `.c` kept beside its `.so` -- is not built, by reason: a unit's key
is the hash of its exact text, so an object served for the current source was compiled from
exactly the files its markers name; a kept copy would add a file per object (4566 in an 11 GB
store) for nothing a report needs.  The HDL kernels (generated text, no file in the tree) are where
a kept `.c` would help; left until a report in one needs it.

Proof: the machine code of all ten units identical (no `-g`, `#line` emits nothing).
test_c_abi gains: an `#error` planted before the last line of newton.c, of newton.h and of
pcx_abi.h is reported at that file and line (alive: with the markers stripped the report names
`unit.c:193`; a first plant at a file's middle sat in a block comment, where `#error` is inert).
12 pinned pairs moved on the twin only.  Full suite 5053 passed; the sanitizer's fast tier clean
(4762 passed).

R1 (stages R1a-R1d) is done: the C units are files, one interface text each, shared types and
codes named once, complete without the HDL prelude, and reported at their own lines -- every step
proved by unchanged machine code (radau_tc's one difference: the prelude's false `const` on
`hypot`).
