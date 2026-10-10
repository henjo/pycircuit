
## 2026-10-10 — 105. Stage 1: the sources and their breakpoints in C (`csrc/sources.c`, `_tran_sources.Table`)

Sized first: on spice_4049 the source pass is 48.5 k a call, served by the direct pass
(`_stamp_plan.source_direct`) -- its validity check 22.2 k (what no watcher sees: the elements'
dicts and types, the classes' `u` and `function`, `ParameterDict`'s lookup), the evaluation
16.8 k.  So the C sources save little per PYTHON step; their place is the C loop (stage 4), where
they run every step and their validity once a loop call.  A prerequisite, as the plan says.

`csrc/sources.c` (+ `sources.h`): SEVEN waveforms -- the plan said five; DC (the base
`TimeFunction`, `f` the int 0), Sin, Pulse, PWL, Exp, AM, SFFM -- in the Python's own operations
and order: Python's float `%` and `//` (CPython's `float_rem` and `_float_div_mod`) for the
pulse's period, the `where` chain's last true branch with each branch's arithmetic as written,
`_sin_scalar` and its guards, Exp's three branches, the AM / SFFM products as Python groups
them; the direct pass's per-source `v = a + f(t)` and its scatter in element order from +0.0
(a voltage source's `-v` into its branch row, a current source's `+v`, `-v`); each
`next_event` (the pulse's candidates from the period index, strictly after `t`).  DECLINE,
NEVER DIFFER: anything the Python makes by another route bails -- a non-finite intermediate,
a value past the fast add's 1e300 guard, a sine past `_sin_scalar`'s guards, a PWL time that
finds no segment (the Python's None raises).

`_tran_sources.table_for(sd)`: the table of a direct pass, its sources in the pass's own order
and rows, built only where every value comes the fast float way (the waveform exactly one of
the seven, the numeric toolkit's `exp` / `sin` / `where` / `pi`, every parameter a float or a
small int, the DC value and the waveform not both ints -- the Python adds those as ints, an
exponential's taus nonzero -- the Python raises -- numpy's exp / sin the C library's on this CPU
by the HDL backend's probe, `libm_check`).  Per call: the toolkit's functions, the sine's
switch, numpy's underflow state.  `Table.holds()`: every input as read (each source's
`function`, its class, its attributes, a PWL's points, the DC value).

Tests (`test_tran_sources.py`): every waveform V and I on a 3001-point grid and its corners
(+-1 ulp): `u` byte for byte and the next event, none declined; a hypothesis twin per waveform
(signed zeros, subnormals, 1e299 .. 1e301, infinities, NaN, ints among the values; times
negative and past the period) -- wherever the C answers it is the Python's (185 vectors
compared, 53 bails, 63 declines at build in the gate profile); each decline; a table whose input
changed does not hold; a parameter a ulp off is caught.  The unit in `test_c_abi` (complete,
`-Werror`, its header its cdef); a pinned pair "the sources and their C" (the waveforms' `f` /
`next_event`, `_sin_scalar`, `_SourceDirect.run` against the C and its table).

NOT wired into the Python step's source pass: the benchmark decks hold one DC source each
(spice_4049, spice_ring), whose pass is plumbing, not waveforms -- wired only if it counts FEWER
against an A/A, which a single DC source will not show.  Stage 4 calls the table.
Full suite 5084 passed; the sanitizer's fast tier clean (4816 passed).
