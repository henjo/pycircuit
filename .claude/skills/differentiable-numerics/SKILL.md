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
