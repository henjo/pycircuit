---
name: documenting-moving-code
description: Write and maintain documentation, docstrings, test names and project records for a codebase that is still changing. Use when writing docs or worked examples, when recording a known gap or limitation, when naming a test after a claim, and whenever you are about to trust a note you found in the repo.
---

# Documentation for code that is still moving

**A note that records the state of the code is trusted like a
measurement, and it is not one.** It was true when written, nothing
tells you when it stops being true, and the reader who follows it has no
way to know. That is the whole of this skill.

Everything below is measured from one campaign, where four such notes
were found stale in a single session and one of them cost an entire
misdirected investigation.

## A recorded gap expires, and nothing announces it

The four, and what each did:

- a test named **`..._terms_are_absent`**, whose body only ever asserted
  that the CARD's coefficients were nonzero — true either way. The terms
  had been implemented weeks earlier. It was written as a help, "so the
  next person measuring a residual knows where to look", and it worked
  exactly as designed on a residual it had nothing to do with;
- a **class docstring** listing nine terms as "still absent" that were
  all built, and quoting an accuracy figure four orders of magnitude out
  of date;
- the project's own **research record**, current to seven commits back,
  with a headline table two days stale — and the memory file pointed at
  it as the place the plan lived;
- a **published number** taken from before a reference regeneration,
  which would have shown a 0.02% disagreement in a document whose
  headline claim was 1e-6.

The pattern is the same every time: the note is *load-bearing for
someone else's decision* and has no owner once written.

**The habit that catches it: after finishing a piece of work, re-read
what the repo says about the state.** Not the code — the prose. Test
names, docstrings, the record, the README's numbers. Ask of each: was
this written before or after what I just did?

See `validation-design`, "Recorded decisions have a shelf life", for the
same failure in a deliberate trade-off rather than in a note.

## Name things by what they do now

A test name is an assertion, and it decays like any other. So does a
variable name, a dict key, and a section heading.

`..._terms_are_absent` kept passing after the terms went in, because its
body checked something else. The name was the only thing carrying the
claim, and names are not executed.

Two defences, both cheap:

- **make the body assert what the name says.** If the name claims a term
  is modelled, assert that it is — `inspect.getsource` on the function is
  a legitimate test when the alternative is an unchecked adjective;
- **prefer names that stay true.** "What this measures" survives; "what
  is missing" does not.

A key that lies is the same defect one level down: a dict entry named
`Gmob` holding `Gmob_dL` reads correctly at every call site and is wrong
at exactly one.

## Run every example — do not write it

Two code samples written for a documentation page, both plausible, both
wrong: a sweep API called with a signature it does not have, and a
result accessor that does not exist. Neither would have survived being
run once.

This is not a lapse in care; it is what happens when prose and code are
written in the same pass. **Transcribe examples from a run**, in that
direction. If it is too awkward to run, it is too awkward for the reader
to run.

## Re-measure every number at publish time

Numbers in prose come from memory or from an earlier run, and both go
stale silently. The reference in this campaign was regenerated twice
mid-session; a value quoted from before that would have contradicted the
document's own headline.

**Prefer documentation that regenerates.** Sphinx's `exec-rst`, doctests,
a table built by a script at build time — a number that is recomputed
cannot rot, and one that is typed will. Where a live block needs
something a reader may not have (a PDK, a licence, a GPU), give it a
recorded fallback and say which one is showing.

## When you cannot build it, verify structurally — and know the gap

Without the doc toolchain installed, hand-verification is still worth
doing: check heading underlines, directive indentation, and *run the
executable blocks standalone*. That caught a short underline and two
dead APIs.

It could not catch what the build caught the moment it was available: a
docstring in an unrelated module whose displayed formula was parsed as
markup, because an indented block after a paragraph is a blockquote and
`|x - y|` reads as a substitution reference. **Autodoc pulls docstrings
into the build, so a mathematical absolute value in a comment is a
documentation error.**

Install the toolchain if you can. Structural checking is the right thing
to do without it and is not a substitute for it.

## What earns a place in a record

Write down the measurement, not the adjective. "Up to 1.73 in
saturation, 0.2% in the linear region" is a regression guard that a
later reader can check; "we do not model X well" is a claim with no
expiry and no test.

And record the wrong turns. In a campaign of any length the successes
are mostly the same few techniques applied patiently; the errors are
where the transferable content is, and they are the entries a reader
returns to.
