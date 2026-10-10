
## 2026-10-10 — 95. Stage R1d, step 2: the codes named once, in the headers (machine code identical)

The numbers the C units return and are called with were bare in the C and mirrored by hand in
Python: each unit's results (`BAIL` tables in five modules, `status != 1`, `status == 8`), the
core's passes (`want`: 1 G, 2 C, 4 i; the fused mask's q bit 8; `_PASS_BITS`, `WANT`), its
formulas (`FORMULA`: 0 gear2 .. 3 theta, 4 the passes alone), the error test's reference modes
(`_RELREF`) and CBLAS's `101, 111`.  A real trap among them: "the limiter walk stopped" is 8 in
newton.c and radau.c and 10 in radau_tc.c -- the Python tested each by its own number.

Now each unit's header names its results (`enum newton_status { NEWTON_CONVERGED = 1,
NEWTON_BAIL_X = 2, ... }`, a line saying what each means; lte.h also `enum lte_relref`), and
pcx_abi.h names the core's (`enum pcx_pass`, `enum pcx_formula`, `enum pcx_cblas`).  The C says
the names; Python reads its tables from the header text (`_cbuild.constants`, `_cbuild.bails`):
`BAIL = _cbuild.bails(_CODES, 'NEWTON_BAIL_')`, `status != _CODES['NEWTON_CONVERGED']`,
`FORMULA` / `WANT` / `_PASS_BITS` from pcx_abi.h -- the strings the Python reports (`'walkstop'`,
`'getrf'` ...) are the names' suffixes, unchanged.

Proof: the machine code of all ten units identical again; every Python table equal to its old
literal (the eight compared as data); cffi reads each header's enums as `constants` does.
test_c_abi gains: the two readers agree on every header; pcx_pass's bits are `1 << k` of the
pass order Python builds the fused masks from (`_tran_core.PASSES`, `_hdl_cbackend.FUSE_PASSES`);
no `.c` says a code by its number (the check finds 80 in the previous commit's C, 0 now); an
enumerator spelling no value refuses.  11 pinned pairs moved on the twin only (pass.c and walk.c
among them: they include pcx_abi.h, which gained the enums -- the watcher's, which does not,
stayed).  Full suite 5040 passed; the sanitizer's fast tier clean (4776 passed).
