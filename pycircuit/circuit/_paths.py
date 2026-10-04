"""THE FAST-PATH COUNTERS (2026-10-04; testing for development, stage 2).

Every fast path of the simulator -- the evaluate core (`_tran_core`), the
pass batches (`_hdl_batch`), the stamp plan, the limiter walk and kernel,
the fused passes, the zero-source skip, the memos -- steps aside WITHOUT
CHANGING AN ANSWER: that is its contract, and every bit-identity check
holds when one is quietly switched off.  The run is only slower.  So each
decline is counted here with its reason, and so is each served call at the
coarse sites (a core call, a batch pass, a walk, a memo hit); the gate
records each test's counts (`benchmarks/tranrec`, family `paths`) and
compares them like the outputs, so a fast path that stopped serving is
named, with the test and the reason.

Keys read `<path>.<site>:<outcome>` (`core.fj:shadow_cir`,
`batch.G:served`, `walk:unbound`).  A key beginning `once:` is a fact
decided once per class, plan or process (a backend resolution, a core's
build): it lands on whichever test a worker runs first, so a comparison
reports those and never fails on them.

`COUNTS` is a container (the leak detector snapshots scalars only), and a
decline reads `return _paths.no(key)`: counted, and None.  A served count
on a hot path binds `COUNTS` to a module global and increments it in
place.  The cost, counted in executed bytecodes (deterministic): +0.6 % on
the 20-MosLevel1 step, +1.1 % on the PSP stage's -- about a dozen
increments a step.  Nothing here changes what any path computes.
"""
import collections
import types

#: every count since the process started; a test's are a difference.  A
#: `defaultdict(int)`, a C type: an increment takes the dict's own subscript
#: (~36 ns; a `Counter`, a Python subclass, took ~73)
COUNTS = collections.defaultdict(int)


def no(key):
    """Count `key` (a decline); a decline returns this call's None."""
    COUNTS[key] += 1


def genuine(obj, qual, mod):
    """`obj` is what `mod`'s source defines as `qual` -- not a caller's
    stand-in (a lambda, a wrapper: another name or module, or `__wrapped__`).
    A fast path that stands in for Python code reads that code's pieces
    once, and only when every one is genuine: one read under a caller's
    patch would hold the patch (the C error test's tests found it)."""
    return (getattr(obj, '__module__', None) == mod and getattr(obj, '__qualname__', None) == qual
            and (isinstance(obj, type) or (type(obj) is types.FunctionType
                                          and '__wrapped__' not in obj.__dict__)))


def snapshot():
    return dict(COUNTS)


def since(before):
    """The counts added since `snapshot()` returned `before` (zeros
    dropped)."""
    out = {}
    for k, v in COUNTS.items():
        d = v - before.get(k, 0)
        if d:
            out[k] = d
    return out
