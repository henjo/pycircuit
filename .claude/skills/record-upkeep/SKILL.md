---
name: record-upkeep
description: Keep a project's records in step with its code — which documents, memory files, docstrings and tests go stale after a change, and in what order to refresh them. Use after any commit that changes behaviour, adds or removes a capability, or overturns a recorded claim; and before reporting a piece of work finished.
---

# Refreshing the record after a change

Writing durable documentation is `documenting-moving-code`. This is the
narrower, duller, more often skipped job: **after the code changes,
which records are now lying?**

Measured on one campaign: the plan document was falsified by work the
plan itself asked for **within a day**, twice; twelve line references
went stale, were re-verified, and were all stale again the next day; a
capability table understated the DSL by six rows for weeks; a memory
file said a branch was "complete and pushed" ten commits later; and a
model docstring explained that it was written a certain way because of
a bug that had just been fixed.

None of that was carelessness. Records rot silently because **nothing
fails when they do** — which is exactly why it needs a checklist rather
than good intentions.

## The trigger

Run this when a change:

- adds, removes or fixes a **capability** (not an internal refactor);
- **overturns a claim** written down anywhere — including one written
  this week, including one you wrote;
- changes a **measured number** that appears in prose;
- makes a **known gap** narrower, wider or gone.

And always before saying a piece of work is finished.

## The sweep, in dependency order

Do it in this order, because each step can invalidate the next.

**1. The claim you just falsified.** Find it *first*, before writing
anything new — you probably knew it existed. Correct it in place with a
date, keep the wrong version briefly, and say how it was wrong. The
wrong version is what stops the next reader re-deriving it.

**2. Specifications and capability maps.** These are the records people
plan work from, so a stale one costs someone a week. Check the
"implemented" list and the "not implemented" list — the second rots
faster and in the more expensive direction.

**3. The plan or roadmap.** Status tables, sequencing, "not started"
rows, and any suite/benchmark number quoted in prose. Add what the work
found: an item now done, a new open item with its measurement, a claim
the work disproved.

**4. Cross-references and anchors.** Grep the docs for references to
what you touched. Prefer **symbol names over line numbers** — line
numbers are a measurement of something nobody is holding still. Verify
each named symbol still exists; a rename is invisible until someone
greps.

**5. Docstrings that explain a workaround.** When the thing being
worked around is fixed, the docstring becomes actively misleading —
it now teaches a constraint that no longer exists. These are the
hardest to find because they live nowhere near the fix; grep for the
symptom, not the fix.

**6. Tests named after a claim.** A test called
`..._is_not_supported` or `..._still_hides_...` is an assertion about
the world with an expiry date. When the claim flips, the test's *name*
is wrong even if it passes. Some are written deliberately to expire —
honour that and invert them rather than deleting them.

**7. Persistent memory.** Status entries ("complete", "pushed", "N
tests") and the one-line index that summarises them. These are read at
the start of every future session with no code nearby to contradict
them, so a stale one is believed longer than any other kind.

**8. Skills.** Only if the work taught something transferable. Prefer
extending an existing skill; a new one earns its place when the domain
is genuinely distinct and there is more than one hard-won finding in
it.

## What to write

**The measurement, not the adjective.** "225 Jacobian evaluations
against 25" survives; "convergence was poor" cannot be re-checked and
cannot expire.

**The wrong turn, with its mechanism.** In any campaign the successes
are a few techniques applied patiently; the errors are where the
transferable content is. A correction that only states the new truth
throws away the more useful half.

**Numbers that let a fix be measured.** When recording something you
are NOT fixing, write down what it costs, so whoever picks it up can
tell whether they improved it.

## Traps

⚠ **Do not update a test just to make it pass after a deliberate
change.** Work out what property it was protecting, keep that, and drop
only the incidental. If you cannot tell the difference, you cannot tell
a regression from an improvement either — see `validation-design`.

⚠ **Suspect the reason before the rule.** A note saying "do X because
Y" has two lifetimes and Y is checked by nobody. When a rule looks
obsolete, check whether only its explanation rotted; deleting a correct
rule because its stated reason expired is the expensive mistake.

⚠ **A record you did not check is not a record you verified.** Saying
"docs updated" after touching one file, when four are stale, is worse
than saying nothing, because it stops the next person looking.

## Closing the loop

State plainly what you refreshed and what you deliberately left alone.
"The roadmap and memory are current; `hdl.md` §7 still describes the
old surface and I did not touch it" is a useful sentence. "Everything
is up to date" is a claim you almost certainly cannot support.
