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
`1.14x` for weeks; re-measured on the same machine it was **1.22x**
(median of four runs — a single run said 1.25x, and quoting that would
have overstated the drift). The benchmark that produced it was still in
the tree and still passing — it asserted the two waveforms AGREED, never
that the ratio held.

So: when you touch a code path, **re-run the published figures that
measure it**, even when your change was not about performance. And
prefer a figure with a test around it — assert a bound the number must
stay under, so the drift fails a run instead of aging in prose.
**Published figures rot in the flattering direction**, because that is
the direction nobody re-checks.

⚠ **A SPEEDUP FIGURE EXPIRES WITHOUT ANYONE TOUCHING IT.** "This cache
is worth 1.16x" is a fact about the cache *and everything else in that
call*, so it moves when any of the rest does — here 1.16x became 1.24x
in an hour with no edit to the cache. It is the only staleness where
**re-running the same benchmark on the same code gives a different
answer**, so it defeats every check based on "did anyone touch this?".
Write "re-run this" beside such a figure rather than a number to quote.
(`measure-before-building` has the mechanism.)

*How* to re-take a figure so the new one is trustworthy — interleaving,
pairing the variants in one process on a busy machine, and why a
profiler's per-call times are not the ones to publish — is
`measure-before-building`. This skill is about noticing that the
published one has rotted.

⚠ **A RATIO THAT DID NOT MOVE IS NOT A NULL RESULT.** When the published
figure compares two things that share a code path, a change to the
shared part speeds up *both* and the ratio sits still. Reading that as
"nothing happened" buries a real win — here a 10% simulator-wide gain
showed as an unchanged 1.15x. Check the absolute numbers on both sides
before concluding a change did nothing.

⚠ **A SPLIT TEST-RUN RECIPE GOES STALE THE MOMENT YOU ADD A TEST FILE.**
If the suite is run in chunks derived from a file count, adding a file
during the work silently drops it out of every chunk — it is then
untested while the totals still look right. This has now happened twice
in this tree, the second time *within a single session*, to a file
written that same hour. Derive the ranges again after adding any test
file, and **verify coverage rather than trusting the arithmetic**: diff
the union of the chunks against the full file list, and reconcile the
outcome count against `--collect-only` (they differ by exactly the
module-level `importorskip` skips, which collection does not count).

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

## Audit the NAVIGATION SURFACES, not the body

Three stale records were found in one session here, and all three were
the same kind of thing: a status table, a resume block, and a suite
count. The reasoning underneath every one of them was sound. **Records
rot at the top**, in the parts written to be read first and skimmed
thereafter.

The worst case measured: a plan document's Stage 4 status table listed
ten items as "not started". Every one had landed, been superseded, or
been refused -- **by commits dated the same day the table was written**.
It was written mid-campaign, the campaign continued past it, and the
table was never revisited. It then stood for four weeks as the first
thing any reader of that document saw.

⚠ **A stale "not started" is far more expensive than a stale "done".**
A stale "done" costs someone a re-check. A stale "not started" costs
them the work itself. Two rows there were of the second kind: one item
read "not started" when its own entry measurements had already refuted
it -- inviting someone to build a thing that had been disproven -- and
another read "needs items 4a-4e" for a parameter that no longer existed
and whose name now raises `TypeError`.

**So when a campaign closes an item, edit the table in the same commit**,
not at the end of the session. And when picking up any document, date
its navigation surfaces before trusting them: `git blame` the status
table and compare against `git log` since. If the table is older than
the last commit to the subsystem it describes, audit it before planning
from it.

⚠ **This is a review cost, not just untidiness.** Where the bottleneck
is somebody else reading the branch, a document that overstates open
work spends their attention on items that do not exist.

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

## "Nothing references it" is not "nothing is lost"

Clearing nine committed scratch files out of a repo root, I ran the
obvious check first: grep every name across the tree, confirm zero
references, delete. That check is necessary and it is not sufficient,
and on its own it would have been wrong twice.

- **It nearly deleted a documented example.** `example10.py` sat among
  the scratch files and looked exactly like them. Two `.rst` files
  include it; removing it would have broken the docs build. The grep
  found it only because the search covered `.rst` and `.md`, not just
  `.py`. **Search the documentation sources, not only the code.**
- **It could not see value.** Two of the nine turned out to be the
  original prototypes of code that later shipped — one of a coupled
  `(dx, dh)` Newton/step-size solve, one of a delay-line ring buffer for
  a traced JAX loop. Nothing was lost by deleting them, but that was
  established by **reading all nine and tracing each idea into the
  current tree**, which took ten minutes and is the only thing that
  could have established it.

The two questions are different and only one of them is about value:

| question | answered by | tells you |
|---|---|---|
| will deleting break something? | grep for references, across code **and** docs | whether it is safe |
| is anything lost? | read it, then find where its idea lives now | whether it is worth keeping |

Report both, separately. "Unreferenced, and its idea shipped as X" is a
sentence someone can approve. "Nothing references it" invites a yes to a
question they did not know they were being asked.

⚠ **A file that is broken can still be broken for the wrong reason.**
Five of those nine set a global (`circuit.default_toolkit`) at *module
scope*, and pytest imports every collected module into one process. One
file therefore failed under collection with a completely different error
than it produced when run alone — a sibling had replaced the global
first. Before drawing a conclusion from a failure in a batch, reproduce
it in isolation.

## A commit message is a record, and the shell can eat half of it

`git commit -m "... `pcnr_vec = dict(...)` ..."` in a double-quoted
shell string runs **command substitution** on everything between the
backticks. The shell prints a syntax-error warning, `git` still commits,
and the message reads:

> *"generate_code then REBUILDS it as  instead of passing it through"*

The sentence survives, grammatical and confident, with its subject
deleted. Nothing fails. The commit is already pushed by the time anyone
reads it.

- **Write the message to a file and use `git commit -F`**, with a tool
  that does no shell interpolation (a Python heredoc, an editor). This
  is the fix, not "remember to escape backticks" -- escaping is one more
  thing to get right per commit, and the failure is silent.
- After committing anything with code in the message, **read it back**:
  `git log -1 --format='%B'`. Check the phrases you meant to write are
  there. Absence leaves no marker to grep for -- the substitution
  removes the backticks too, so you cannot search for what is missing.
- The one detectable trace is a **mid-sentence double space** where the
  content used to be: `grep -nE "[a-z,)]  +[a-z]"` over recent messages
  finds it. Worth running once after a session of code-heavy commits.

⚠ **Fixing it means `--amend` and a force push.** That is fine seconds
after the push and rude later; check `git rev-parse HEAD` against
`origin/<branch>` first, use `--force-with-lease`, and confirm the tree
is unchanged (`git rev-parse HEAD^{tree}`) so the amend is provably
message-only. Ask before force-pushing a branch someone else is reading
-- the message being wrong is cheaper than a reviewer's history moving
under them.
