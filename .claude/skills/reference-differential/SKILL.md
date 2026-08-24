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
- **Is the reference CONVERGED?** If it came out of a solver it carries
  that solver's tolerances. ngspice's defaults left 9.6e-4 of error in
  the golden curves here — the largest single error in the comparison
  once the model was good enough to see it, and it presents as a KINK at
  one bias rather than as an offset, which invites a hunt for a
  mechanism. It cost two false leads: a "threshold shift switching on"
  in the transfer sweeps and a "saturation-knee defect" in the output
  sweeps, the latter sending me after `THESATB`/`THESATG`. Regenerate at
  a tightened tolerance and diff before concluding anything about
  physics. See `validation-design`, "The reference inherits its
  generator's tolerances".
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

**This is the step that ends it.** If the differential says a term is
wrong and every ingredient still matches, the mechanism is invisible to
reading, and you need the reference to tell you its own numbers: add
debug outputs to a COPY of the reference source, rebuild, run one point.

Two bugs here were found this way after ten term-by-term readings found
nothing, and each took a single bias point. Strobe a whole CHAIN of
intermediates at once, not one suspect: the ratios then localise the
error themselves. `alpha`, `Gmob`, `FdL` at 1.000 with `dps` at 0.958
says the error is upstream of `dps`; `zsat` at 0.929 being exactly
`dps²` says it is not a separate error. One run replaced a week of
inference.

**Check first that you can rebuild the reference at all**, with an
unmodified source, before spending time on the patch. Here that check
was skipped and cost an hour: `openvaf-r 20260616` produces OSDI that
**ngspice-47 cannot load — it segfaults, deterministically, on a
pristine rebuild**. The shipped `.osdi` was built with something else.
Bisecting the patch was wasted work; the patch was never the problem.

⚠ **`rc=0` does not prove the model loaded.** ngspice exits 0 on a
MISSING `osdi` file, so a successful exit code is not evidence. Assert
on the output, not the status.

If the reference's own host cannot load your rebuild, use the compiler's
sibling simulator (VACASK, for openvaf-r) rather than fighting the
mismatch.

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
- **⚠ A QUANTITY THE REFERENCE COMPUTES TWICE.** Both final bugs here
  had this shape, and it is invisible to term-by-term reading *because
  every term that is present is transcribed correctly*:
  - one NAME standing for two different quantities — PSP's `Vdsx` is
    `Vds²/(√(Vds²+0.01)+0.1)`, NOT the smoothed `|Vds|` this element
    also needed; using one for both made a log four times too large;
  - one TERM present in one place and missing in the other — the series
    resistance `GR` was in the midpoint mobility and absent from the
    source-end one, 26% of that quantity.

  - one QUANTITY computed once and USED at several sites — PSP's
    conditioned `Vsbstar` feeds both the quasi-Fermi level and the GATE
    DRIVE (`Vgb1 = Vgs + Vsbstar - VFB`); having it at the first site
    only was 331 µV of threshold at `Vsb = 1`, and the largest
    remaining discrepancy in the whole comparison.

  When the reference evaluates the same physics at several points
  (source end, midpoint, drain end), check each SITE, not just each
  formula. And share the helper rather than writing the expression
  twice — two copies is how the sites drift apart.

  **The reading habit this asks for is different from the usual one.**
  Grep the reference for every USE of a quantity's name, not for its
  definition. You verify definitions by reading them; you can only
  verify uses by ENUMERATING them, because the one you are missing is
  by construction the one you are not looking at.
- **Match the reference's own conditions exactly.** The last 0.04% here
  was a **0.15 K** difference between the temperature the parameters
  were scaled at and the one the model was evaluated at. Not a model
  error at all; a harness one. Check temperature, bias convention and
  units before concluding anything about physics — and note that this
  one had been *documented and bounded* as an acceptable trade for
  months. It was acceptable when the model was 1.6% off and stopped
  being so at 3e-5. See `validation-design`, "Recorded decisions have a
  shelf life".
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
- **Reading PSP's own internals** — ngspice-47 cannot load an
  `openvaf-r` build, but **VACASK** can (same author). Convert a COPY of
  the PDK with
  `PDK_ROOT=... PDK=ihp-sg13g2 PYTHONPATH=~/local/vacask/lib/vacask/python
  python3 -m sg13g2tovc`; it rebuilds the PDK **and compiles the
  Verilog-A**, so a patched PSP103 runs. Declare extra `OPP(...)`
  outputs, assign them where the quantity is computed, `save
  p('m1:nsg13_lv_nmos', dbg_x)`, read the `.raw` with
  `rawfile.rawread`. The converter chokes on the local LDE edits in
  `sg13g2_moslv_mod.lib` — use `git show HEAD:` versions in the copy.
- Comparisons must set **`epar.T = 273.15 + 27`**: the parameters are
  scaled at the card's `TR`, and leaving `epar` at its 300.0 K default
  is a 0.15 K mismatch worth 0.04% of the noise density.

## Worked example

The n-channel drain current was 5.5% high in weak inversion, 0.1%
correct above threshold, with an exact subthreshold slope, on one
channel type only.

Ten quantities were read against the vendor and matched. The surface
potential was verified exact to 1e-8 against a 40-digit root. Every
parameter matched the vendor's exported `lp_*`. Nothing absent could
account for it, checked against the card's own switches. The residual
stayed.

Bounding eliminated six candidates in one run and left the term with
the right shape — but no ingredient of it could be off by the required
factor, and saying "it is that term" on shape alone was wrong once.

The differential settled WHERE: zeroing that term in a copy of the card,
read by BOTH models, moved the agreement from 1.0549 to 0.9976.

Reading the reference's own internals settled WHY, in one bias point:
`s2` was 0.14424 against PSP's 0.03597, every other ingredient exact.
The cause was `Vdsx` — which is not the smoothed `|Vds|` its name
suggests. A second strobe found the last residual the same way: the
series resistance was missing from the source-end mobility.

**All twelve sweeps then matched the vendor at median 1.000**, and with
the evaluation temperature matched too, `ids` and the noise density come
to 3e-5 — the reference's own printed precision.

The lesson the whole exercise sharpens: **reading tests your reading.**
Ten careful readings found nothing that two strobes of the reference's
own numbers found immediately.
