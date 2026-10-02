"""EVALUATION SESSIONS: an analysis says which of a circuit's `i`, `q`, `G`,
`C` it is about to request at one state, so that an element able to compute
several of them in one pass -- a compiled HDL model's fused functions,
`_hdl_cse.fused` -- does, and hands the others out from that pass
(2026-10-02; the fused-evaluation plan's F3).

The analyses keep calling the four methods as before (the circuit passes,
the Jacobian counters the tests wrap, PCNR's instance shadows all see the
same calls); a session only lets an element answer the second, third and
fourth from the first one's work.  It must name EXACTLY what the site will
request at that state: an element computes the whole named set on the first
request, and a set larger than the site's needs computes derivative chains
nobody reads (and can raise floating-point warnings the separate calls never
did).  Outside a session nothing changes.

`evaluating(*which)` opens a session; `evaluating(session=s)` re-enters one,
for a site whose requests at one state are split across two scopes (the
multistep step's branch screen reads `C` at the converged point, and its
caller then reads `q` and `G` there).
"""
import contextlib
import contextvars

_CURRENT = contextvars.ContextVar('pycircuit_evaluation_session',
                                  default=None)

#: the method names a session can carry
METHODS = frozenset(('i', 'q', 'G', 'C'))


class Session:
    """The methods an analysis will request at one state.  Its identity
    scopes the elements' memos (`_hdl_cse.take`)."""

    __slots__ = ('which',)

    def __init__(self, which):
        which = frozenset(which)
        if not which <= METHODS:
            raise ValueError(f'evaluation session: unknown method(s) '
                             f'{sorted(which - METHODS)}')
        self.which = which


def current():
    """The session in force, or None."""
    return _CURRENT.get()


@contextlib.contextmanager
def evaluating(*which, session=None):
    """Within the block, the circuit's `i`/`q`/`G`/`C` requests at one
    state are exactly `which` (or `session`'s, re-entered).  Yields the
    session."""
    s = Session(which) if session is None else session
    token = _CURRENT.set(s)
    try:
        yield s
    finally:
        _CURRENT.reset(token)
