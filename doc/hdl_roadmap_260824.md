
## 2026-10-10 — 102. CORRECTION to 100 and 101: `--count` has a 1-1.7 % offset between identical trees -- R2b was neutral, not MORE

An A/A (`step_machinery --count` with the parent tree at HEAD and this tree identical to it,
checked file by file): spice_4049 +1.66 % / +1.69 %, mos1 +1.20 % / +1.22 % (its step +0.89 /
+0.88 %), stage +1.17 % / +1.19 % (its step +1.22 / +1.12 %), mos1_adaptive +1.41 / +1.33 % --
the child MORE, the same within 0.1 % run to run: a SYSTEMATIC offset, not noise.  Roles
swapped (this tree the parent, the worktree the child): -0.2 .. -0.3 % -- so it follows the
role and the tree both (the address layout each process gets with randomisation off: object
hashes, allocation sizes), and moves when either tree's contents change.

So a count difference below ~2 % between two trees is not attributable without an A/A taken in
the same configuration.  Consequences:

* 100 (R2b) said "counted MORE: +1.24 %, +1.0 % a step".  Against this offset that is NEUTRAL
  (~0, within the band).  R2b stays reverted -- no gain was measured, and it adds a context to
  keep exact -- but it was not shown to cost.  Its design still goes to stage 4's entry.
* 101 (R2c batch 1, -3.49 % spice_4049, -2.6 / -2.7 % a step) was counted in the same default
  configuration, whose A/A offset is positive: its true effect is at least that (about -5 % if
  the offset then was today's).  It stands.
* Every later count of the plan names its A/A, taken in the same configuration beside it.
