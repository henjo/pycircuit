"""Negative-result regression test for the naive CNA ``Nullor`` hypothesis.

Tlelo-Cuautle & Sarmiento-Reyes, *A Pure Nodal-Analysis Method Suitable for
Analog Circuits Using Nullors*, JART 2003 (``doc/cna_conclusions.md``,
``doc/cna_plan.md``). Their Fig. 8/9 worked example is an op-amp inverting
amplifier with transfer function ``Vo/Vi = -Rf/Ri`` (their eq. 3).

Stage 1 of the CNA revival plan hypothesized that a single 4-terminal
element stamping an all-zero ``G`` and declaring no ``branches`` -- the
mainline ``Nullor`` minus its branch -- would reach the paper's "Compacted
PNA" matrix size for free, since pycircuit's ``IS`` is already NA-compatible
and never needs the paper's voltage-source-to-nullor transform that inflates
their own node count. Measured here: **refuted**. The naive element is not
a correct nullor at all -- for a nullor with one terminal of each port
grounded (as in this circuit), it produces a genuinely singular reduced
system (``det() == 0``, exactly, in ``sympy``), not merely a less-compact
one. Real compaction needs the paper's row/column merge-and-delete rules
applied at circuit-assembly time (node aliasing/elimination before the
matrix is built), which is a materially bigger feature than an element
class -- see ``doc/cna_plan.md``'s recorded Stage 1 result for the reasoning
and the decision point it raises.

The candidate element lives only in this test (not in ``elements.py`` --
it's provably wrong, so it doesn't belong in the public API) so this
regression stays cheap to run if anyone reconsiders reviving CNA later.
"""
import sympy
from sympy import symbols, simplify

from pycircuit.circuit.circuit import Circuit, defaultepar
from pycircuit.circuit import SubCircuit, R, IS, Nullor, gnd
from pycircuit.circuit.analysis_ss import AC
from pycircuit.circuit.toolkit import symbolic

Ri, Rf, i_in, s = symbols('Ri Rf i_in s', real=True, positive=True)


class _NaiveZeroGNullor(Circuit):
    """The refuted Stage-1 hypothesis: no branch, all-zero G, same terminals
    as the mainline Nullor. Kept local to this test -- see module docstring.
    """
    terminals = ('inp', 'inn', 'outp', 'outn')

    def G(self, x, epar=defaultepar):
        return self.toolkit.zeros((self.n, self.n))


def _build(nullor_cls):
    """Fig. 8's inverting amplifier: Iin -> node 1 -(Ri)-> node 2 -(Rf)-> node
    3, nullor across (node2, gnd)-(node3, gnd). node1 is the input voltage
    node (undriven except through Ri and the current source), node3 is Vout.
    """
    cir = SubCircuit(toolkit=symbolic)
    cir['Iin'] = IS(1, gnd, iac=i_in)
    cir['Ri'] = R(1, 2, r=Ri)
    cir['Rf'] = R(2, 3, r=Rf)
    cir['Nullor'] = nullor_cls(2, gnd, 3, gnd)
    return cir


def test_mainline_nullor_matches_published_transfer_function():
    """Sanity check on the test circuit itself: reproduces eq. 3."""
    cir = _build(Nullor)
    res = AC(cir, toolkit=symbolic).solve(s, complexfreq=True)
    tf = simplify(res.v(3, gnd) / res.v(1, gnd))
    assert simplify(tf - (-Rf / Ri)) == 0, (
        'test circuit does not reproduce the published Vo/Vi = -Rf/Ri; '
        'got %s' % tf)


def test_naive_zero_g_nullor_is_singular_not_merely_less_compact():
    """The refuted hypothesis, pinned: this produces a singular system.

    If this test ever fails (i.e. the naive element stops being singular),
    something about the surrounding G-assembly machinery changed in a way
    that's directly relevant to reconsidering the CNA revival -- see
    doc/cna_plan.md's Stage 1 result before touching this.
    """
    cir = _build(_NaiveZeroGNullor)
    cir.update_iparv()

    import numpy as np
    n = cir.n
    G = sympy.Matrix(cir.G(np.zeros(n, dtype=object)))

    gnd_idx = cir.nodes.index(cir.get_node('gnd'))
    keep = [i for i in range(n) if i != gnd_idx]
    reduced = G[keep, keep]

    assert reduced.det() == 0, (
        'expected the naive zero-G Nullor to give a singular reduced '
        'system for this topology; det() = %s' % reduced.det())


def test_naive_zero_g_nullor_has_no_branch():
    """Confirms the one thing the hypothesis got right: no extra unknown.

    Not useful on its own (the element is wrong), but documents precisely
    what "compact" would have meant if the element had been correct.
    """
    cir_mna = _build(Nullor)
    cir_cna = _build(_NaiveZeroGNullor)

    assert len(cir_mna.nodes) == len(cir_cna.nodes)
    assert len(cir_mna.branches) == 1
    assert len(cir_cna.branches) == 0
