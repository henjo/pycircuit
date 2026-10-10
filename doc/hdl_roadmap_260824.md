
## 2026-10-10 — 100. Stage R2b measured and REFUTED: one run context per call does not pay

R2b's premise: the C Newton's ~28 per-call checks become ONE run context, made once a run and
re-validated by ONE stamp test per call.  Sized first (`insitu.py lines spice_4049
_tran_newton_c:solve`, 408.5 k a call): the C 113.7 k, the branch screen 87.8 k (stage 3's), the
source pass 49.2 k (stage 1's), the run-level re-checks ~95 k.  R2's line was RE-SIZED before R2b
was built (the plan, 2026-10-10): the R2-addressable ~158 k down >= 35 %.  A throwaway ceiling
prototype (every run-level check skipped, unsafely) put the most R2b could save at 46.5 k a call.

Built exactly: a `_Run` context after a call whose every check passed, at the counter's present
value; its test -- the counter, plus ONE pass over what no dict watcher sees (the transient's own
attributes the checks read, the machinery's class and module attributes, the switches,
`ParameterDict`'s epoch, the linear solver's choice, values in unwatched dicts, each element's
`__dict__` object, each pass's code, the circuit's size).  First run: built but held 0 of 400
calls -- the step's source memo (`_RUN_STATE` 'step') is a new dict every step, and the context
compared it by identity where the check asks only that one exists (fixed: held 397 of 400).
Then COUNTED against its parent (`step_machinery --count`): spice_4049 +1.24 %, mos1 and stage
+1.0 % a step -- MORE.  Bytes and statistics the same.

Why: the stamps already skip what a dict watcher can see; what is left on every call is exactly
what NO watcher sees, item by item -- and the fused test must compare every one of those too (~70
identities through mapped tuples), which costs what today's inline checks cost, plus the
context's own reads.  The work a context could truly skip (the limiter closure built to test `is
not None`, the walk's lookup, the core record's call) is a few k.  Reverted, by the plan's rule (a
stage showing MORE is not kept, even for clarity).  The design is not lost: validated ONCE PER
C-LOOP CALL -- many steps a call -- it is the C loop's entry gate (stage 4), where it pays; the
diff is kept with the plan's working notes.
