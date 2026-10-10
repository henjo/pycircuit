
## 2026-10-10 — 98. R1's line held, and the one cost it hid: radau_tc's residual loops reloaded after every `hypot`

R1's line (written before measuring): bit-identical everywhere, `step_machinery --check` no SLOWER,
against the tree before R1 (02235bcc).  Held: all 18 cases "no change", bytes, statistics and
fast-path counts SAME.  But the child leaned slower in most rounds (+0.1..+1.1 %, some intervals
clear of zero: ladder_gear +1.1 % [+0.5, +1.7]) -- under the check's band (+3 % flags every case,
+1 % none), so counted rather than shrugged off: `--count` (instructions, perf) flagged ONE case,
mos1_radau's step +1.00 % (2.759 M -> 2.787 M), with ladder_radau +0.33 % (spread 0.01 %); every
other case within its tree's own spread.

The cause is step 3's one machine-code change.  The prelude declared `hypot` `const`; without that
lie, the compiler must assume each `hypot` call may write memory, so after every call in the
residual check's loops it reloaded the struct's array pointers (`s->Ax`, `s->Ai`, `s->y1`,
`s->ab`, the bound `s->Ap[j + 1]`).  Fixed by saying it honestly: the arrays in local pointers and
the inner bound read once (the arithmetic and its order untouched).  Recounted against the same
parent: mos1_radau's step -0.10 %, ladder_radau -0.07 %, pss_radau -0.03 % -- no MORE anywhere.
Lesson: a time verdict of "no change" with a one-sided lean is worth one `--count`; it named the
case and the cause in a minute.
Two pinned pairs moved on the twin (radau_tc.c) and were re-recorded; the 103 Radau tests pass.  Full
suite 5053 passed; the sanitizer's fast tier clean (4746 passed).
