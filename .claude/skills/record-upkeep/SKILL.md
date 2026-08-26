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

⚠ **A PERFORMANCE FIGURE NO TEST ASSERTS WILL DRIFT, and nothing will
tell you.** The triggers above all assume some change *targeted* the
number. The worse case is the number no change targeted: it decays a
little per commit, every suite stays green, and the prose keeps
claiming the original. `hdl.rst` advertised the DSL's overhead as
`1.14x` for weeks; re-measured on the same machine it was **1.25x**.
The benchmark that produced it was still in the tree and still passing
— it asserted the two waveforms AGREED, never that the ratio held.

So: when you touch a code path, **re-run the published figures that
measure it**, even when your change was not about performance. And
prefer a figure with a test around it — assert a bound the number must
stay under, so the drift fails a run instead of aging in prose.
**Published figures rot in the flattering direction**, because that is
the direction nobody re-checks.

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

## Do not edit a document you have handed to a running agent

An orchestration rule, learned by breaking it.

A model-library agent was told "the tree is clean at `<sha>`" and asked
to update the roadmap as part of its work. While it ran, three more
commits landed on that same roadmap from the orchestrator. The agent
read the file, read it again later, found it 180 lines longer with a
status row it had previously seen saying something else — **and
concluded it had misread its own earlier observation.** It had not.

Nothing was lost (the edits were in different sections and merged
cleanly), but the agent spent real effort reasoning about an anomaly
that was manufactured, and it had no way to tell a moving file from its
own error. That is the expensive part: **an agent cannot distinguish
"the ground moved" from "I was wrong", so it will usually assume the
second.**

- Treat an assigned document like an assigned source file. If you told
  an agent to update `PLAN.md`, `PLAN.md` is theirs until they report.
- If you must edit it, **tell them** — a one-line message costs nothing
  and converts a mystery into a fact.
- State the baseline as a fact that may expire ("clean at `<sha>` as of
  now; I may commit docs while you run") rather than as a standing
  truth.

The same applies to attribution. Sweeping your own edits into an
agent's commit with `git add -A` produces a commit whose message
credits a batch for changes it did not make — including, in that case,
a file the agent had been explicitly forbidden to touch. Stage
deliberately, or say in the message which parts were yours.


## Scope a table edit to the table -- and check the file did not shrink

A status-table row was patched with a regex that replaced text "up to
the next line that starts with a number". The next such line was in a
section 1300 lines below. The commit deleted a third of the plan
document, including the section a running agent was working from, and
nothing failed: no test reads a document. The agent noticed only
because its reference section had vanished.

- Slice the document to the table (or the section) FIRST, edit inside
  the slice, then splice it back. Never search past the region you
  mean to change.
- Before and after editing a record, compare something cheap that a
  deletion would move: the line count, or the list of `## ` headings.
  Print both. A record edit that shrinks the file by more than the
  edit is a deletion, whatever the script intended.
- Keep the previous commit's copy around until the diff has been read.
