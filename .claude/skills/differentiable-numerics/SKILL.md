---
name: differentiable-numerics
description: Write and debug numeric code that is compiled through symbolic differentiation — compact models, DSL-generated elements, anything where a Jacobian is produced from an expression. Use for NaN or non-finite Jacobians, overflow/invalid-value warnings from generated code, "the value is right but the derivative is wrong", and before adding a new math primitive.
---

# Numerics for code that is symbolically differentiated

## The two things that make this different

**1. The derivative is a different program from the value.** It is
generated, it is bigger, and it divides by things the value never
divides by. A model can be finite everywhere and have a NaN Jacobian at
one bias — which is the *worse* failure, because nothing looks wrong and
Newton is silently poisoned.

**2. Both arms of every conditional are evaluated.** A compiled
`Piecewise` becomes `where`/`select`, which picks *afterwards*. So the
discarded arm still overflows, still raises floating-point flags, and
still contributes its NaN to any expression that combines the arms.
Each arm's INPUT must be clamped to that arm's own domain.

## Rules

### Clamping fixes the value and breaks the derivative

`sqrt(max(x, 0))` differentiates to `d(max)/dx / (2*sqrt(max(x,0)))` —
`0/0` for every `x <= 0`. Smooth instead of clip. Reach for
`hypsmooth` first; `safe_sqrt`, `safe_ln`, `safe_div`, `safe_abs` are
built on it.

### `Abs` is a quiet trap

Sympy does not know a model symbol is real, so it differentiates
`Abs(u)` as `(re(u)*re'(u) + im(u)*im'(u)) / Abs(u)` — `0/0` at zero.
Use `safe_abs`. Raising an `Abs` to a power compounds it: the printed
derivative carries `sign(u)/(u*Abs(u))`.

### `Max`/`Min` go on ATOMS, never on expressions

Differentiation expands the argument into every branch condition. When
the argument contains a smoothing function, sympy rationalises it there
into a form whose denominator cancels to exactly zero. The value is
finite; the Jacobian divides by zero. Bind the argument with `var()`
first.

### Floor the base of a power

`rat**ax` differentiates to `ax*rat**ax/rat`. Floor `rat` at a small
positive value whenever the base can vanish, and always when the
exponent is below 1 (the derivative diverges).

### A "safe_" helper spends exponent range -- do not buy what you do not need

`safe_ln` regularises as `log(0.5*(sqrt(u^2 + 4e-60) + u))`. It SQUARES
its argument, so it overflows for `u > 1e154` where plain `log` would
not have blinked. `safe_div` squares its `eps` and `hypsmooth` squares
its; every one of them halves the exponent range it can carry.

That is the right trade when the guard can fire. It is a pure loss when
it cannot. Here the arguments were `1 + expl(...)`, at least 1 by
construction, and using `safe_ln` on them cost a finiteness bound of
1e36 and returned `-inf` from a quantity whose every ingredient was
finite.

**Before reaching for a `safe_` helper, ask what it is guarding against
and whether that case is reachable.** If the argument is provably
positive, the plain function is both correct and wider.

### A PRODUCT of two bounded exponentials is not bounded

`expl` is finite for every double, so `expl(a)` and `expl(b)` are both
safe -- and `expl(a)*expl(b)` is not. Past its seam `expl` continues
polynomially rather than saturating, so each factor reaches ~1e186 at
extreme bias and the product overflows.

When the code only ever needs the LOGARITHM of such a product, carry
the exponents instead and never form it: `log(1 + e^z)` is bounded by
`|z| + 0.7`, where `log(1 + e^a * e^b)` is bounded by nothing. Written
branch-free as `max(z,0) + log(1 + exp(-|z|))`, the exponential only
ever sees a non-positive argument.

### A term switched OFF by a zero parameter is still evaluated

If the model compiles from symbolic parameters, `if param != 0` is TRUE
at build time whatever the value, so the block is in the expression
regardless. Passing zero only multiplies the result by zero at the end
-- and `0 * inf = nan`.

Two consequences:

- **a block must be finite even when its parameter is zero.** "It is off
  on this card" is not a defence;
- **zeroing a parameter is not a valid way to test whether a block is
  the cause of a numerical failure.** It cost a wrong conclusion here:
  the finiteness bound fell from 1e36 to 1e26, switching the parameter
  off did not restore it, and the block was the cause anyway. Compare
  against a build WITHOUT the code -- `git stash` -- not against a build
  with the code and a zero coefficient.

### Budget the exponent range twice: value and derivative

A regulariser that squares its denominator halves the usable range. Its
*derivative* squares that again and quarters it — so the Jacobian dies
two decades of exponent before the residual does. State both ranges in
the docstring; they are different numbers.

### Underflow is not automatically the safe direction

Dropping a term from a *sum* can flip a sign. `-2*b*inv**2 + inv`
underflows the first term for large `b` and returns `+1/b^2` where the
answer is `-1/b^2`: finite, plausible, wrong. Prefer groupings whose
factors are bounded — `inv*(1 - 2*b*b*inv)` contains
`b^2/(b^2+eps^2)`, which is in `[0,1]` for every real `b`.

### Hold every intermediate — the let-chain IS the arithmetic

Bare sub-expressions get substituted into one another before printing,
so the compiled form evaluates products the written form never creates
(`nu**4` against `a*nu*tau`). The same algebra typed into a REPL can
return a finite number where the compiled chain returns NaN. `var()` is
not an optimisation here; it fixes the order of operations.

### Write a derivative in terms of itself when the natural form loses range

If the generated derivative overflows where the value does not, make the
quantity a primitive with a hand-written `fdiff` expressed via `self`.
That turns an overflowing `denominator**2` into an underflowing
`reciprocal**2`.

### A new primitive needs four things

1. `fdiff` — or sympy emits an unevaluated `Derivative` that kills
   `lambdify`;
2. `_imp_` — sympy's hook for "this has a numeric implementation", so a
   **plain** `lambdify` works with no modules map. Public helpers are
   lambdified by callers outside your module; requiring a namespace
   entry makes your change breaking;
3. a `_print__name` method on any strict printer (`_ChainPrinter` is a
   `NumPyPrinter` and raises on unknown functions);
4. a namespace entry **where the namespace is built**, not at each call
   site that passes a modules map. Relying on every future call site is
   not a design.

## How to find these

- **Finite differences against the analytic Jacobian**, swept across the
  domain — not at one point. Kinks hide.
- **Warnings as errors** (`np.errstate(all='raise')`, `simplefilter`)
  over a wide bias sweep. A green suite carrying overflow warnings is
  telling you something; trace them.
- **Sweep to absurd bias** and scan for the first non-finite, separately
  for value and Jacobian, on both signs. Scan for the boundary; never
  assert that something breaks.
- **Use CARD parameters, not defaults.** Defaults leave optional blocks
  off and every quantity smaller. Three separate bugs here were
  invisible with defaults.
- **Dump the compiled intermediates.** Take `fn._src`, re-`exec` it with
  a dict capture after every assignment, and find the first non-finite
  variable. That converts "somewhere in 2000 lines" into one name.

## One expression cannot branch on the solution

A reference model may ORDER its terminals with a real `if` that swaps
`Vgs`/`Vsb`/`Vds` and records a sign. A compiled expression cannot: the
branch condition depends on the unknowns, so the sign has to be carried
arithmetically instead — typically through a smooth `|Vds|`.

That substitution is correct and necessary, and it creates a specific
hazard: **the smooth `|Vds|` is now a quantity the reference does not
have**, and it is easy to reuse it wherever the reference writes
something with a similar name. Here the reference's own `Vdsx` is
`Vds²/(√(Vds²+0.01)+0.1)` — a softened drain bias, not an absolute value
— and using one quantity for both roles cost 5.5% of the current.

Whenever you replace a branch with arithmetic, write down which of the
reference's quantities your new one stands for, and check every OTHER
use of that name separately.

## A smoothing constant is meaningless without its scale

`hypsmooth(x, e)`, `safe_div(a, b, eps)`, any regulariser — the constant
is compared against something. Find out what, and write the constant
RELATIVE to it.

A channel-shortening clamp here used `hypsmooth(x, 1e-3)`, and the
smoothed quantity was immediately divided by a scaled parameter `VP`
and put through a logarithm. So the answer depended on `e/VP` — and
`VP` was **0.322 on one channel type and 7.38e-6 on the other, a factor
of 44000**. The same constant was 0.3% of one scale and 135x the other,
and it was the entire remaining error on the second: 2.5e-4 against
5.3e-6.

Two things follow:

- **an absolute constant encodes an assumption about a parameter you do
  not control.** It was defensible on the card it was tuned against and
  indefensible on the next one. Write `1e-3 * vp`, not `1e-3`;
- **pick the LOOSEST value that has converged, not the tightest that
  works.** A wider smoothing is a rounder corner and a better-behaved
  Jacobian, so scan until the answer stops moving and then step back to
  the knee. Here `1e-2*VP` was still costing accuracy and `1e-4*VP`
  bought 1% over `1e-3*VP`, so `1e-3*VP` is the pick.

Scan the constant over decades and tabulate. If the answer moves, the
constant is doing physics rather than regularisation, and its value is a
model parameter you have not admitted to having.

## Symmetric in exact arithmetic is not symmetric in floating point

An exact structural property -- antisymmetry, conservation, a
cancellation the model relies on -- can be broken by a CANCELLATION even
when the algebra is perfect.

Here a correction was written as `vsb - vsbst`, a difference of two
quantities differing by 3e-4. At `Vsb = 1` that keeps ten digits; at
`Vsb = 1e40` it keeps none, and the exact source/drain antisymmetry
broke by 2% out there. The algebra was symmetric; the arithmetic was
not.

**Form a small quantity WHERE IT IS SMALL.** The conditioning above is
`mina(vlow, 0, aphi) - phix1`, bounded by the smoothing constants and
of order 1e-4 at every bias. Computing it there and letting the large
quantity be derived FROM it (`vsbst = vsb - vsbcnd`) inverts the
dependency and the cancellation never happens.

The general form: when `big_a - big_b` is known to be small, there is
almost always an expression for the small thing that does not mention
the big ones. Find it. And prefer holding the small quantity as the
primitive, deriving the big one, rather than the reverse.

Test it at ABSURD magnitude, not merely at large. The bug here was
invisible at 1.5 V and obvious at 1e7.

## Warnings in a green suite are unfinished work

24 `invalid value encountered in scalar divide` warnings sat in this
suite for a long time, raised by two ordinary two-transistor circuits.
They were a `0/0` being evaluated on the way to real solutions. Tracing
them was what turned an "absurd bias only" curiosity into a real fix.
