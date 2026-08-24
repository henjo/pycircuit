---
name: reference-differential
description: Find the root cause when an implementation disagrees numerically with a reference implementation (vendor model, golden data, a second implementation) and line-by-line source comparison keeps passing. Use for "our model is N% off and every formula matches", for compact-model work against a PDK, and whenever a residual survives an audit.
---

# Finding a numerical discrepancy that source reading will not find

## The principle

**A textual comparison and a numerical one are different tests, and when
they disagree the numerical one is right.**

Reading your code against the reference tests *your reading*. It cannot
find an error that lives in how a term is used, in a value that matches
but is applied to the wrong quantity, or in a formula that is
transcribed perfectly and still produces a different number. When ten
quantities have been read against the vendor and matched and the
residual is still there, stop reading. Measure.

## The procedure

Work in this order. Each step is cheap relative to the one after it, and
each has found something here.

### 1. Audit the comparison itself, before the model

The most expensive mistake is measuring the wrong thing carefully.

- **Is the reference what you think it is?** A deck that instantiates
  `X1 ... some_device` may be instantiating a **subcircuit**, not a bare
  model: wrappers insert series sources, compute instance parameters,
  and pass geometry/stress arguments the model then uses. Read the
  subcircuit.
- **Is the measured quantity what you think it is?** Reference data is
  usually *terminal* current, which includes every path to that
  terminal. Check the other terminals: if `i_g` and `i_b` are 1e-14
  against an `i_d` of 1e-8, the drain current is channel current and you
  may proceed.
- **Is the window clean?** Reference curves have floors — leakage,
  numerical noise, a term you do not model. Measuring below one means
  measuring *its* floor, not the physics. Find the floor and stay above
  it.
- **Are the two sides at the same operating conditions?** Temperature,
  bias convention, geometry, instance parameters.

### 2. Bound each candidate — do not fit

Turn each suspect term off, one at a time, and see how far the number
moves. One run eliminates most of the field.

**Leverage is not error.** A term large enough to explain the residual
is a *candidate*; it is not the cause. Conflating the two is the single
most common failure in this kind of work, and the most seductive,
because the arithmetic looks like a proof.

Read the eliminations as hard as the hits:

- a term whose removal changes nothing is out;
- a term that also acts in a regime where you already agree is out
  (it would have broken that agreement too);
- a term whose removal breaks something currently exact is *right*.

### 3. Choose an observable where the candidates differ qualitatively

Fitting a shape in one variable rarely separates two mechanisms —
**check collinearity before fitting**. Over a realistic range, candidate
basis functions are often 99%+ collinear, and then no amount of
precision separates them and the best-fitting one is not thereby the
true one.

Look instead for an observable on which they differ in *kind*:

- gain terms scale a curve; electrostatic terms change its slope;
- some terms act only in one regime (weak vs strong inversion, linear
  vs saturation);
- some touch charge and not current, or vice versa.

Convert the discrepancy into the unit the physics lives in. A current
*ratio* in an exponential regime is a *voltage*:
`dV = ln(ours/ref) / (d ln I / dx)`. That says both how large the error
is in millivolts and — from whether it is flat — whether it is an offset
at all.

### 4. The differential against the reference

**This is the step that works when reading has failed.**

Modify the *shared input* — the model card, the parameter file, the
config — so a suspect term is disabled, and let **both** implementations
read the modified input. Then compare.

- If the disagreement vanishes, the discrepancy is in that term.
  Causal, not correlational, and established without reading a line.
- If it persists, the term is eliminated no matter how well its shape
  matched.

Then get the factor: hold the reference at its real value and scan
**your** parameter until the curves coincide. The value that lands is
the reference's effective contribution as a fraction of yours.

### 5. Read the reference's own internals

If the differential says a term is wrong and every ingredient still
matches, the mechanism is invisible to reading and you need the
reference to tell you its own numbers: add printf/`$strobe` of the
intermediates to a copy of the reference source, rebuild, run one point.

## Record the factor; do not apply it

A scale factor that lands the number without a mechanism **hides the bug
rather than fixing it**. Pin it as a test with the measurement that
found it, and say plainly that the mechanism is open.

## Traps, each of which cost real time here

- **A residual is not an error.** Root error is `residual / |dF/dx|`. If
  the equation is flat where you are measuring, an impressively small
  residual can sit on a large displacement. Validating a solver by its
  residual says it solved the equation, not that the root is right.
- **A tidy, interpretable answer is not thereby a correct one.** A
  contaminated window here produced a clean, flat 0.5% that read exactly
  like a thermal-voltage error. The clean window gave 0.01%.
- **Two measurements that cannot both be true are a finding.** Do not
  reconcile them by inventing a mechanism that makes the contradiction
  sound like physics. Follow the contradiction; it is pointing at the
  error.
- **A parameter that matches is not a term that is right.** Matching the
  reference's exported parameter proves the scaling, not the use.
- **A zero coefficient does not test a scaling**, and a bias no sweep
  visits tests nothing. Two cards is not twice the validation — it is
  the difference between exercising a term and assuming it.
- **A key that lies is worse than a long name.** A dict entry named
  `Gmob` holding `Gmob_dL` reads correctly at every call site and is
  wrong at exactly one.
- **Do not assert a failure in a test.** "This breaks above 1e28" fails
  the day someone stops it breaking. Scan for the boundary and assert it
  has not got worse.

## Harness notes for this repo

- Reference data: `pycircuit/circuit/tests/data/psp103_ihp_iv.json`,
  generated by `benchmarks/psp_reference.py`. Regeneration must be
  **additive** — check every existing number reproduces bit-identically.
- The vendor exports its scaled parameters as `lp_*` operating-point
  outputs. Comparing term by term against those is the first tool to
  reach for; `benchmarks/psp_reference.py:LP_PARAMS` is the list.
- The vendor also exports quantities it does not name directly: two
  outputs plus one shared quantity can determine a third
  (`sig` + `cigid` ⇒ `CGeff`).
- Reading our own compiled intermediates: take `_hdl_info['funcs'][...]`,
  split `fn._src`, and re-`exec` it with a dict capture after every
  assignment. That gives every chain variable at a bias.
- Differential setup: `cp -r` the PDK to scratch and edit the one card
  file. **`shutil.copytree` produces a copy ngspice cannot use.**
- ngspice **segfaults (rc = −11)** on some intermediate parameter
  values. Scan your own model against the reference's recorded curve
  rather than re-running the vendor per point.
- PSP subthreshold: the reference carries a ~2e-12 A junction-leakage
  floor this core does not model. Use the window **1e-9 to 1e-6 A**.

## Worked example

The n-channel drain current was 5.5% high in weak inversion, 0.1%
correct above threshold, with an exact subthreshold slope, on one
channel type only.

Ten quantities were read against the vendor and matched. The surface
potential was verified exact to 1e-8 against a 40-digit root. Every
parameter matched `lp_*`. Nothing absent could account for it, checked
against the card's own switches. The residual stayed.

Bounding eliminated six candidates in one run and left `ALP2`, whose
term is weak-inversion-specific by construction — the right shape. That
was *not* sufficient: no ingredient could be off by the required factor,
and saying "it is the ALP2 term" on shape alone was wrong once already.

The differential settled it in two runs. Zeroing `ALP2` in a copy of the
card, read by **both** models, moved the agreement from 1.0549 to
0.9976. Scanning ours against the vendor's fixed curve put the crossing
at 0.286 — our term is 3.5× too strong.

Mechanism still open. Recorded as a factor, pinned by a test that needs
no vendor tooling, with the next step named: strobe `dL1`/`FdL` out of a
rebuilt `.osdi`.
