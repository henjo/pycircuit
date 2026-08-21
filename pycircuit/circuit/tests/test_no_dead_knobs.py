"""The dead-knob instrument (doc/transient_review_260820.md, cross-cutting test).

Two parts, each honest about what it catches:

1. a per-function unused-argument scan over the transient-related modules --
   the half that found F18 (and `ywr_error_ratio`'s dead `x_last`) on first
   contact;
2. a parameter-reachability walk for Transient and JAXTransient, scoped across
   the module AND the base-class module, with the string-variable forwarding
   loop handled by treating literal name-tuples as references.

What this deliberately does NOT catch: options that are read but dishonoured
(F2, F5, F9, F10's shapes).  Those are guarded by their behavioral tests and
by nothing cheaper.

Every allowlist entry carries the reason it is allowed.  An entry whose reason
has expired (a scheduled fix landed) should be removed with that fix.
"""

import ast
import os

HERE = os.path.dirname(__file__)
CIRCUIT = os.path.normpath(os.path.join(HERE, '..'))

SCAN_MODULES = ['transient.py', 'jaxtransient.py', 'integrator.py',
                'stepcontroller.py', 'nrsolver.py', '_lte_kernels.py']

## (module, function, argument) -> why it is allowed to ignore the argument.
UNUSED_ARG_ALLOWLIST = {
    ## Deleted by F7 (Phase 2); remove this entry with that fix.
    ('transient.py', 'fang_timestep', 'refnode'):
        'dead; scheduled for deletion by F7',
    ## Helper signatures kept uniform with their sibling that does use ctrl.
    ('transient.py', '_band_centre', 'ctrl'): 'signature uniformity',
    ('transient.py', '_lte_in_band', 'ctrl'): 'signature uniformity',
    ## The Integrator ABC fixes one signature for all methods; lower-order
    ## integrators legitimately ignore history they do not difference.  Each
    ## case is documented at the implementation.
    ('integrator.py', 'check_order_drop', 'h_curr'): 'ABC signature (Euler)',
    ('integrator.py', 'check_order_drop', 'h_last'): 'ABC signature',
    ('integrator.py', 'check_order_drop', 'is_first_step'): 'ABC signature',
    ('integrator.py', 'compute_derivatives', 'iq_last'): 'ABC signature',
    ('integrator.py', 'compute_derivatives', 'h_last'): 'ABC signature',
    ('integrator.py', 'compute_derivatives', 'is_first_step'): 'ABC signature',
    ('integrator.py', 'compute_derivatives', 'toolkit'): 'ABC signature',
    ('integrator.py', 'companion_dh', 'h_last'): 'ABC signature',
    ('integrator.py', 'compute_lte', 'h_last2'):
        'ABC signature; Euler differences one past point (documented there)',
    ## StepController.evaluate_step is one interface for three estimator
    ## families; each ignores the inputs the other families need.
    ('stepcontroller.py', 'evaluate_step', 'x_hist'):
        'charge-based controllers ignore the solution history',
    ('stepcontroller.py', 'evaluate_step', 'h_clamped'):
        'PIController has no lower band to suppress (F10 revisits)',
    ('stepcontroller.py', 'evaluate_step', 'q_curr'):
        'SolutionLTEController is charge-free by design',
    ('stepcontroller.py', 'evaluate_step', 'q_last_hist'):
        'SolutionLTEController is charge-free by design',
    ('stepcontroller.py', 'evaluate_step', 'iq_last_hist'):
        'SolutionLTEController is charge-free by design',
    ('stepcontroller.py', 'evaluate_step', 'J'):
        'SolutionLTEController needs no J^-1 mapping',
    ('stepcontroller.py', 'evaluate_step', 'irefnode'):
        'SolutionLTEController takes the argmax over the full vector',
    ## NonLinearSolver interface; SchurCoupledNewton and JAXNewtonSolver do
    ## not dispatch on these (JAXNewtonSolver is F15 delete-or-port).
    ('nrsolver.py', 'solve_system', 'scaler'): 'interface; see F15',
    ('nrsolver.py', 'solve_system', 'row_names'): 'interface; see F15',
    ('nrsolver.py', 'solve_system', 'linsolver'): 'interface; see F15',
    ('nrsolver.py', '_structural_singularity', 'toolkit'):
        'signature uniformity with its sibling',
    ('nrsolver.py', '_worst_row_report', 'toolkit'):
        'signature uniformity with its sibling',
}

## Parameters declared but not referenced in (module, base module) -- each with
## the reason.
PARAMETER_ALLOWLIST = {
    ('JAXTransient', 'vabstol'):
        'F6(a) moved the Newton floor to iabstol; F6(b) re-threads vabstol '
        'per-row (Phase 3) -- remove this entry with that fix',
    ('JAXTransient', 'analysis'):
        "inherited; the traced loop hardcodes analysis='tran' -- unify when "
        'the loop learns other analysis names',
}


def _trivial_body(fn):
    body = [n for n in fn.body
            if not (isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant))]
    if all(isinstance(n, ast.Pass) for n in body):
        return True
    return len(body) == 1 and isinstance(body[0], ast.Raise)


def test_no_function_accepts_an_argument_it_never_reads():
    offenders = []
    for mod in SCAN_MODULES:
        path = os.path.join(CIRCUIT, mod)
        tree = ast.parse(open(path).read())
        for fn in ast.walk(tree):
            if not isinstance(fn, ast.FunctionDef) or _trivial_body(fn):
                continue
            args = [a.arg for a in fn.args.args + fn.args.kwonlyargs
                    if a.arg not in ('self', 'cls')
                    and not a.arg.startswith('_')]
            loaded = {n.id for n in ast.walk(fn) if isinstance(n, ast.Name)}
            for a in args:
                if a not in loaded and (mod, fn.name, a) not in UNUSED_ARG_ALLOWLIST:
                    offenders.append('%s:%d %s(%s)' % (mod, fn.lineno, fn.name, a))
    assert not offenders, (
        'arguments accepted and never read (the F8/F18 defect class):\n  '
        + '\n  '.join(offenders))


def _par_references(paths):
    names, literal_tuples = set(), []
    for path in paths:
        tree = ast.parse(open(path).read())
        for node in ast.walk(tree):
            if (isinstance(node, ast.Attribute)
                    and isinstance(node.value, ast.Attribute)
                    and node.value.attr == 'par'):
                names.add(node.attr)
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == 'getattr' and node.args):
                a0 = node.args[0]
                if (isinstance(a0, ast.Attribute) and a0.attr == 'par'
                        and len(node.args) > 1
                        and isinstance(node.args[1], ast.Constant)):
                    names.add(node.args[1].value)
            ## the string-variable forwarding loop (_solve_operating_point):
            ## a literal tuple of parameter names IS a set of references, and
            ## treating it so also guards the tuple itself against typos.
            if (isinstance(node, ast.Tuple) and len(node.elts) >= 4
                    and all(isinstance(e, ast.Constant)
                            and isinstance(e.value, str) for e in node.elts)):
                literal_tuples.append([e.value for e in node.elts])
    return names, literal_tuples


def _assert_parameters_reachable(cls, module_files):
    declared = {p.name for p in cls.parameters}
    paths = [os.path.join(CIRCUIT, f) for f in module_files]
    refs, tuples = _par_references(paths)
    for t in tuples:
        ## Only a tuple consisting ENTIRELY of declared parameter names is a
        ## forwarding tuple -- a partial overlap is some other literal tuple
        ## (TransientStatistics.__slots__ shares 'max_step' with the
        ## parameters, for instance) and must not be misread as forwarding.
        ## A forwarding tuple with a typo therefore does not count as
        ## references, and its genuine names get flagged dead -- which still
        ## surfaces the typo, just under the other assertion.
        if set(t) <= declared:
            refs |= set(t)
    allowed = {name for (c, name) in PARAMETER_ALLOWLIST if c == cls.__name__}
    dead = declared - refs - allowed
    assert not dead, (
        '%s declares parameters nothing reads: %s' % (cls.__name__, sorted(dead)))


def test_every_transient_parameter_is_read():
    from pycircuit.circuit.transient import Transient
    _assert_parameters_reachable(Transient, ['transient.py', 'analysis.py'])


def test_every_jaxtransient_parameter_is_read():
    from pycircuit.circuit.jaxtransient import JAXTransient
    _assert_parameters_reachable(JAXTransient, ['jaxtransient.py', 'analysis.py'])
