"""Shooting tests: shooting modal.  Split out of test_analysis_shooting.py on
2026-09-27 (in its original order); shared helpers are in
`_shooting_fixtures.py`, the HDL elements in `_shooting_elements.py`.
"""
from pycircuit.circuit import *
from pycircuit.circuit.shooting import (PAC, algebraic_conditioning,
                                        topological_index)
import warnings
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter as _HdlParameter, white_noise)
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest
import functools as _functools
from pycircuit.circuit.tests._shooting_fixtures import (_a9_vdp,
    _coloured_vdp,
    _orbit_modulated_vdp)


def test_the_orbital_covariance_resolves_onto_the_floquet_modes():
    """A9 step 2: `K_orb` resolved onto the Floquet directions.

    The two routes we already own meet here — `oscillator_covariance`
    gets `K_orb` from a bordered Kronecker solve, `floquet_modes` gets the
    eigen-directions from the monodromy — and Traversa & Bonani's eq (22)
    sums over exactly these mode pairs.

    ⚠⚠ **THE RECONSTRUCTION CANNOT BE EXACT, AND THAT IS STRUCTURAL, NOT A
    TOLERANCE.** A DAE monodromy has annihilated (null) directions, which
    `floquet_modes` drops; on this fixture that leaves **2 modes against a
    4-wide covariance**, so `U cw U†` is rank ≤ 2 and `K_orb` is not.
    Measured residual 1.9e-3 relative — the part of `K_orb` living in the
    slaved algebraic directions. ⚠ Asserting machine precision here would
    be asserting that a rank-2 object equals a rank-4 one; the honest
    gates are the ones below.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 400, x0=np.array([2.0, 0.0]),
                  maxiterations=300)
    assert pss.converged

    pac = PAC(cir)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cw, _wi = pac.orbital_mode_weights(pss)
        modes, K = _wi['modes'], _wi['K']
    U = np.column_stack([m['u0'] for m in modes])
    V = np.column_stack([m['v0'] for m in modes])

    ## 1. the basis is biorthonormal -- everything else rests on it
    bio = float(np.max(np.abs(V.conj().T @ U - np.eye(len(modes)))))
    assert bio < 1e-10, \
        'the Floquet basis is not biorthonormal (max|V^H U - I| = %.3e), ' \
        'so the projection weights are not what they claim to be' % bio

    ## 2. ⚠ THE PHASE MODE CARRIES ESSENTIALLY NO ORBITAL WEIGHT. This is
    ## what `oscillator_covariance`'s split MEANS -- the along-orbit
    ## growth `n d uu^T` has been removed, so what remains should not sit
    ## on the phase direction. Measured 1.56e-19 against 1.26e-05.
    assert abs(cw[0, 0]) < 1e-8 * abs(cw[1, 1]), \
        'the phase mode carries orbital weight %.3e against the amplitude ' \
        'mode\'s %.3e. `oscillator_covariance` is supposed to have taken ' \
        'the along-orbit growth out, so a large value here means the ' \
        'split leaked' % (abs(cw[0, 0]), abs(cw[1, 1]))

    ## 3. and the AMPLITUDE mode accounts for the covariance
    rec = U @ cw @ U.conj().T
    rel = float(np.linalg.norm(rec - K)) / float(np.linalg.norm(K))
    ## ⚠⚠ AND THIS BOUND IS A PROPERTY OF WHERE THIS FIXTURE PUTS ITS NOISE,
    ## NOT OF THE METHOD.  The basis omits the ANNIHILATED modes, and on the
    ## same circuit noised in a FAST branch instead of at the oscillator node
    ## they carry 99.96% of `K_orb` -- see
    ## `test_the_orbital_mode_basis_is_complete_only_for_noise_in_the_slow_subspace`.
    ## Read this as "the sibling fixture injects into the slow subspace", not
    ## as "the non-null modes account for the covariance".
    assert rel < 1e-2, \
        'the retained modes capture only %.3f of K_orb; if this has grown, ' \
        'the covariance has significant support outside the non-null ' \
        'Floquet directions and a modal orbital spectrum would be ' \
        'incomplete' % (1.0 - rel)
    assert abs(abs(cw[1, 1]) / np.linalg.norm(K) - 1.0) < 5e-2, \
        'the amplitude mode no longer accounts for the orbital covariance ' \
        '(weight %.3e against ||K_orb|| %.3e)' \
        % (abs(cw[1, 1]), np.linalg.norm(K))


def test_orbital_correlation_is_gated_three_ways():
    """A9 step 3: eq (22)'s `C_lhj`, with eq (23) as the gate.

    Three routes to `R_yy(0)`, the stationary transverse covariance:
      A. `orbital_correlation` — the modal Fourier sum, eq (22).
      B. the DEFINITION — a 1-D Lyapunov integral along the single orbital
         mode, no Fourier machinery, built here from the same modes.
      C. `oscillator_covariance`'s samples, cycle-averaged with the
         along-orbit growth removed — shares nothing with A or B.

    ⚠ A ≈ B to 1e-3 says the transcription of eq (22) is right. A ≈ C says
    the modes and `CY/2` are right. It was C failing by 1.75× — while A
    and B agreed — that isolated a scale defect to the shared input and
    found `floquet_modes` mis-normalising the state block.

    ❌ The 3 % shape residual against C is OPEN and bounded here at 5 %,
    not tuned away; the clean comparison needs the plain-path Lyapunov
    solve, which is still gear-only.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 400, x0=np.array([2.0, 0.0]),
                  maxiterations=300)
    assert pss.converged
    pac = PAC(cir)
    m = cir.n - 1

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        R4, _ = pac.orbital_correlation(pss, maxharmonics=4)
        R, Cc = pac.orbital_correlation(pss, maxharmonics=32)
        Kf, info = pac.oscillator_covariance(pss, samples=True)
        modes = pss.floquet_modes(pss)
    assert np.linalg.norm(R - R.T) < 1e-12 * np.linalg.norm(R), 'R not symmetric'
    assert np.linalg.norm(R - R4) < 1e-6 * np.linalg.norm(R), \
        'the harmonic sum has not converged by H=4 on van der Pol (%.3e)' \
        % (np.linalg.norm(R - R4) / np.linalg.norm(R))

    ## B: the definition. Single real orbital mode.
    orb = [k for k, md in enumerate(modes) if abs(abs(md['lam']) - 1.0) > 1e-6]
    assert len(orb) == 1
    md = modes[orb[0]]
    assert abs(np.imag(md['mu'])) < 1e-12
    mu2 = float(np.real(md['mu']))
    P2 = np.real(md['p'][:, :-1]); Q2 = np.real(md['q'][:, :-1])
    Nn = P2.shape[1]; Tp = float(pss.period); h = Tp / Nn
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        CY2 = 0.5 * np.real(np.asarray(pac._cy_reduced(pss, 0.0)))
    g = np.array([Q2[:, k] @ CY2 @ Q2[:, k] for k in range(Nn)])
    nper = max(int(np.ceil(-40.0 / (2 * mu2 * Tp))), 1)
    taus = np.arange(0, nper * Nn) * h
    wts = np.exp(2 * mu2 * taus) * h
    sig2 = np.array([float(np.sum(wts * g[(k - np.arange(0, nper * Nn)) % Nn]))
                     for k in range(Nn)])
    Rdef = np.mean(np.stack([sig2[k] * np.outer(P2[:, k], P2[:, k])
                             for k in range(Nn)]), axis=0)
    relAB = np.linalg.norm(R - Rdef) / np.linalg.norm(Rdef)
    assert relAB < 2e-3, \
        'eq (22) sum and the definition integral disagree by %.3e; they ' \
        'share only the modes, so this is the transcription' % relAB

    ## C: cycle-mean TRANSVERSE Lyapunov covariance -- the OBLIQUE projection
    ## Pi K Pi^T with Pi = I - u v^T/(v^T u), which is Demir's v1^T y = 0.
    ## ⚠ Subtracting only the secular growth is NOT the transverse part: it
    ## leaves the phase direction's bounded within-period variance and read
    ## 2-6 % against this sum, falling as 1/Q. That was the reference being
    ## the wrong object, and it cost an afternoon.
    Ps = [np.asarray(P, float)[:m, :m] for P in info['samples']]
    G = [np.asarray(gg, float)[:m, :m] for gg in info['growth_samples']]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        v0, pinfo = pss.ppv()
    ## `samples[j]` is node j: prepending `v0` shifted every node by one
    ## (the 'first-order Lyapunov route' of 2026-09-20 was this shift)
    vs = [np.asarray(sv, float)[:m] for sv in pinfo['samples']]
    proj = []
    for j in range(min(len(Ps), len(vs))):
        w, Uj = np.linalg.eigh(G[j])
        uj = Uj[:, np.argmax(w)] * np.sqrt(max(float(w.max()), 0.0))
        den = float(vs[j] @ uj)
        if abs(den) < 1e-300:
            continue
        Pi = np.eye(m) - np.outer(uj, vs[j]) / den
        proj.append(Pi @ Ps[j] @ Pi.T)
    Pm = np.mean(np.stack(proj), axis=0)
    ratio = np.linalg.norm(R) / np.linalg.norm(Pm)
    relAC = np.linalg.norm(R - Pm) / np.linalg.norm(Pm)
    assert abs(ratio - 1.0) < 2e-3, \
        'magnitude against the projected Lyapunov cycle-mean is %.6f' % ratio
    assert relAC < 3e-3, \
        'the eq (22) sum disagrees with the obliquely-projected Lyapunov ' \
        'covariance by %.3e (expected ~1e-4 quadrature). If this has grown ' \
        'to a few percent, the projection has been dropped and the phase ' \
        'direction\'s bounded variance is back in the reference' % relAC

def test_the_orbital_residual_was_the_reference_not_the_sum():
    """A9's 2-3 % residual, CLOSED by the plain-path wiring, in two steps.

    First the wiring falsified the pair-artefact story: on euler-plain,
    `n = m`, no pair, the residual against the growth-subtracted Lyapunov
    cycle-mean is still 2.4 %. Then the correct reference removed it: the
    transverse covariance is the OBLIQUE projection `Pi K Pi^T`, Demir's
    `v1^T y = 0`, and against THAT the eq (22) sum agrees to ~6e-4.

    ⚠ THE SIGNATURE THAT NAMED IT: the growth-subtracted residual falls as
    1/Q_lambda (5.4 / 2.4 / 1.2 / 0.6 % at Q = 4 / 8 / 16 / 32) -- orbital
    variance ~ Q against a CONSTANT phase-direction bounded part. Neither
    "physics" (which would grow with Q) nor "numerical" (flat).
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
    pss = PSS(cir, method='euler', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, timestep=T / 1600, x0=np.array([2.0, 0.0]),
                  maxiterations=400)
    assert pss.converged
    pss.monodromy = 'native'   # this test measures the one-step method's OWN monodromy
    assert pss.factored_period().width == cir.n - 1, 'expected n = m'
    pac = PAC(cir)
    m = cir.n - 1
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        Rm, _ = pac.orbital_correlation(pss)
        Kf, info = pac.oscillator_covariance(pss, samples=True)
        v0, pinfo = pss.ppv()
    Ps = [np.asarray(P, float)[:m, :m] for P in info['samples']]
    G = [np.asarray(g, float)[:m, :m] for g in info['growth_samples']]
    ts = np.asarray(info['times'], float)[:len(Ps)]
    Tp = float(pss.period)
    ## the WRONG reference, kept as the documented signature
    Pg = np.mean(np.stack([Ps[j] - (ts[j] / Tp) * G[j]
                           for j in range(len(Ps))]), axis=0)
    rel_g = float(np.linalg.norm(Rm - Pg)) / float(np.linalg.norm(Pg))
    ## the RIGHT reference
    ## `samples[j]` is node j: prepending `v0` shifted every node by one
    ## (the 'first-order Lyapunov route' of 2026-09-20 was this shift)
    vs = [np.asarray(sv, float)[:m] for sv in pinfo['samples']]
    proj = []
    for j in range(min(len(Ps), len(vs))):
        w, Uj = np.linalg.eigh(G[j])
        uj = Uj[:, np.argmax(w)] * np.sqrt(max(float(w.max()), 0.0))
        den = float(vs[j] @ uj)
        if abs(den) < 1e-300:
            continue
        Pi = np.eye(m) - np.outer(uj, vs[j]) / den
        proj.append(Pi @ Ps[j] @ Pi.T)
    Pp = np.mean(np.stack(proj), axis=0)
    rel_p = float(np.linalg.norm(Rm - Pp)) / float(np.linalg.norm(Pp))
    assert rel_p < 3e-3, \
        'against the obliquely-projected reference the sum is off by %.3e; ' \
        'expected ~6e-4' % rel_p
    assert 1e-2 < rel_g < 5e-2, \
        'the growth-subtracted reference reads %.3e off; it is supposed to ' \
        'be 2.4%% here -- the phase direction\'s bounded variance. If it has ' \
        'vanished, oscillator_covariance\'s split changed; if it has grown, ' \
        'so did that variance' % rel_g
    assert rel_p < rel_g / 10.0, \
        'projecting did not remove most of the residual (%.3e -> %.3e)' \
        % (rel_g, rel_p)


def test_the_modal_spectrum_reads_a_coloured_source_per_input_sideband():
    """`modal_spectrum` with a STATIONARY coloured source (2026-09-25).  Input
    sideband `m` carries the source at ``w - m w0``; it now reads `CY`
    there, as `pnoise` does, instead of one `CY` for every sideband (which
    the refusal guarded).  The phase widths take `c` from the white part
    alone.  On `_coloured_vdp` (400 points, H = 8, 16 sidebands) at +3 /
    +10 / -10 f_amp, `total / (pnoise/2) - 1`:

        white      +5.2e-4  +5.7e-4  +4.4e-4    (the modal floor, as before)
        coloured   +5.2e-4  +5.7e-4  +4.4e-4
        filtered   +5.4e-4  +6.5e-4  +3.8e-4

    and coloured/filtered totals +1.3e-4 (the `pnoise` pair's own
    agreement).  ⚠ The PHASE / ORBITAL split differs between the two
    realisations (phase 3.95e-9 against 3.32e-9 at 10 f_amp): the filtered
    circuit's filter state enters its PPV, so its phase is a different
    coordinate.  The totals are physical; the split is not.  One `CY` for
    every sideband (the model the refusal protected) reads 2.3 .. 5.1x.
    Below the phase model's validity (`phase_psd`'s corner / power bound)
    it refuses."""
    import warnings as _w
    tot = {}
    for kind in ('coloured', 'filtered'):
        _c, pss, pac, ov = _coloured_vdp(kind)
        f0 = 1.0 / float(pss.period)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            _v, info = pss.ppv()
            f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
            offs = np.array([3.0, 10.0, -10.0]) * f_amp
            ms = pac.modal_spectrum(pss, offs, ov, maxharmonics=8, maxsidebands=16)
            pn = np.array([float(np.real(pac.pnoise(pss, f0 + o, ov,
                                                    maxsidebands=16,
                                                    sweeptype='absolute')[0]))
                           for o in offs])
        ratio = ms['total'] / pn
        assert np.max(np.abs(ratio - 1.0)) < 1e-3, (kind, ratio)
        tot[kind] = ms['total']
        if kind == 'coloured':
            with pytest.raises(ValueError, match='linearised skirt'):
                pac.modal_spectrum(pss, np.array([1e-12]), ov, maxharmonics=8)
    assert np.max(np.abs(tot['coloured'] / tot['filtered'] - 1.0)) < 5e-4


def test_the_modal_spectrum_takes_a_coloured_source_that_follows_the_orbit():
    """`modal_spectrum` and `phase_psd` with a COLOURED source whose level
    follows the orbit (2026-09-26).  Each component is a unit process
    through its own columns `G(t)`: the rows are the harmonics of the
    PRODUCT `q_l^T G` (the convolution `R_p = sum_k T_{p-k} G_k` over
    every harmonic on the grid), band `p` weighted by the colour at
    `|w - p w0|`; `c(f)` likewise with the PPV.  `G` is the element's
    SIGNED amplitude where it states one, else the root of its PSD; a
    per-band colour takes its root per band frequency.

    Gated against the stationary realisation (`_orbit_modulated_vdp`),
    gear, at 3 / 10 / -10 f_amp:
      * a signed flicker `k V_v flicker_noise(1)`: every part <= 6.1e-13,
        `phase_psd` 1.5e-12 (2.2e-13 on an asymmetric orbit);
      * a Lorentzian at `(k V(v, b))^2` (per band, a level that keeps its
        sign): parts <= 8.4e-14, `phase_psd` 3.4e-14 -- and no sign warning;
      * the corner MOVING with V: against `pnoise(cyclostationary=True)`,
        the same quasi-static model by another fold, on radau 1.2e-13 (with
        no white part the phase widths are 0 and the two sums coincide).
    The PSD-specified flicker (`(k V_v)^2 / f`, V_v changes sign) is warned
    on, and reads 0.6 / -0.9 of the parts -- the `|m|` process.
    Poisons: the modulation frozen to its mean 1.00; the convolution's sign
    flipped 0.82."""
    import warnings as _w

    def run(kind, method='gear'):
        _c, pss, pac, ov = _orbit_modulated_vdp(kind, method=method)
        f0 = 1.0 / float(pss.period)
        with _w.catch_warnings(record=True) as caught:
            _w.simplefilter('always')
            _v, info = pss.ppv()
            f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
            offs = np.array([3.0, 10.0, -10.0]) * f_amp
            ms = pac.modal_spectrum(pss, offs, ov, maxharmonics=8, maxsidebands=16)
            ph = pac.phase_psd(pss, np.abs(offs[:2]))
        blind = any('touches zero along the orbit' in str(x.message)
                    for x in caught)
        return ms, ph, blind, (pss, pac, ov, offs, f0)

    for ref, kind in (('flicker_ref', 'flicker'), ('lorentz_ref', 'lorentz')):
        A, B = run(ref), run(kind)
        for k in ('phase', 'orbital', 'correlation', 'total'):
            err = np.max(np.abs(B[0][k] / A[0][k] - 1.0))
            assert err < 1e-9, (kind, k, B[0][k] / A[0][k] - 1.0)
        ## `phase_psd` is frequency-aware by default (2026-09-26): each
        ## side is a bordered GMRES solve on a different circuit, 2.4e-7
        ## on the Lorentzian pair (the DC fold agrees to 3e-14)
        assert np.max(np.abs(B[1] / A[1] - 1.0)) < 1e-6, (kind, B[1] / A[1])
        assert not B[2], kind
    assert run('flicker_psd')[2], 'a sign-blind root went unwarned'
    ms, _ph, _b, (pss, pac, ov, offs, f0) = run('lorentz_moving', method='radau')
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pn = float(np.real(pac.pnoise(pss, f0 + offs[1], ov, maxsidebands=16,
                                      cyclostationary=True, sweeptype='absolute')[0]))
    assert abs(ms['total'][1] / pn - 1.0) < 1e-9, (ms['total'][1], pn)


def test_the_floquet_modes_carry_a_source_on_an_algebraic_node():
    """⚠⚠ Until 2026-09-26 a one-step method's `floquet_modes` gave the
    adjoint `q` EXACTLY 0 on an algebraic node: the replay carries `C^T q`,
    `pinv(C^T)` recovers `q` from it, and on a node with no capacitance its
    minimum-norm choice is 0.  `modal_spectrum` (and every consumer of the
    modes) then read a source there as absent, silently -- exactly 0 for
    a source into `n` on radau, trbdf2, trap and the GLM (their twin);
    -4.9e-3 against pnoise for a source across `n` and the tank.  Gear's
    transposed solve carries the DAE adjoint and was right.

    An algebraic state's column of the adjoint equation holds no
    derivative, whatever the mode's exponent, so its entry is SLAVED -- the
    PPV's constraint fill (`_algebraic_adjoint_fill`), now at each sample
    of each mode.  Measured (van der Pol, the source into `n` times V_v):
      * `q_n = -k V_v q_v` at every sample of every mode to 2.8e-16 on
        radau, trbdf2, trap, glm3 and a non-uniform radau grid -- gear's
        own transposed solve satisfies it to 5.6e-16 (4.2e-16 on its
        non-uniform continuous adjoint), so the fill's sign is gear's;
      * radau's modal parts against the modulated element (no algebraic
        node): 1.2e-14;
      * a source ACROSS `n` and the tank against pnoise: -1.5e-5 / -1.6e-6
        at 3 / 10 f_amp (was -4.9e-3); on the asymmetric orbit -5.9e-5 /
        -1.35e-4 at +-10 f_amp.
    Poisons: the fill dropped (1.00; -1.5e-3 / -8.3e-3 across); its sign
    flipped -- INVISIBLE on the symmetric orbit (`(1 - kV) q` is `(1 + kV)
    q` shifted by T/2, the same power spectrum), +5.1e-3 / -5.2e-3 on the
    asymmetric one, which is why the across case runs there."""
    import warnings as _w
    c, pss, pac, ov = _orbit_modulated_vdp('white_ref', method='radau')
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        modes = pss.floquet_modes(pss)
    names = [str(x) for x in c.nodes]
    irn = pss.irefnode
    kn = names.index('n') - (names.index('n') > irn)
    kv = names.index('v') - (names.index('v') > irn)
    V = np.asarray(pss.waveform[1], dtype=float)[names.index('v')]
    for md in modes:
        q = np.asarray(md['q'])
        N = q.shape[1] - 1
        err = np.max(np.abs(q[kn, :N] + 0.05 * V[:N] * q[kv, :N]))
        assert err <= 1e-12 * np.max(np.abs(q[kv])), (md['lam'], err)
    res = {}
    for kind in ('white_ref', 'white'):
        _c, pss, pac, ov = _orbit_modulated_vdp(kind, method='radau')
        f0 = 1.0 / float(pss.period)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            _v, info = pss.ppv()
            f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
            offs = np.array([0.3, 3.0, 10.0, -10.0]) * f_amp
            res[kind] = pac.modal_spectrum(pss, offs, ov, maxharmonics=8, maxsidebands=16)
    for k in ('phase', 'orbital', 'correlation', 'total'):
        err = np.max(np.abs(res['white_ref'][k] / res['white'][k] - 1.0))
        assert err < 1e-9, (k, res['white_ref'][k] / res['white'][k] - 1.0)
    _c, pss, pac, ov = _orbit_modulated_vdp('white_across', method='radau',
                                            a=0.3)
    f0 = 1.0 / float(pss.period)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        _v, info = pss.ppv()
        f_amp = -np.log(float(info['second_multiplier'])) * f0 / (2 * np.pi)
        offs = np.array([10.0, -10.0]) * f_amp
        ms = pac.modal_spectrum(pss, offs, ov, maxharmonics=8, maxsidebands=16)
        pn = np.array([float(np.real(pac.pnoise(pss, f0 + o, ov,
                                                maxsidebands=16,
                                                sweeptype='absolute')[0]))
                       for o in offs])
    ratio = ms['total'] / pn
    assert np.max(np.abs(ratio - 1.0)) < 5e-4, ratio


def test_oscillator_covariance_takes_a_coloured_source_in_its_transverse_part():
    """Coloured noise in `oscillator_covariance` (2026-09-25; D1, D2).  A
    coloured source's phase does not diffuse, so it has no bounded-plus-
    growth split; it enters the TRANSVERSE covariance alone -- the deflated
    forced responses, replayed from their bounded part and projected at
    every node (`_transverse_responses`), over the band, the grid adapting
    to the orbital lines (0.02 f0 wide at each harmonic the mode couples:
    `_orbital_lines`).  `K_orb`, `d`, `c_from_growth` stay the WHITE
    sources' (here none: all exactly 0), with a warning.

    Gate, two routes that share only the component model: the cycle-mean
    transverse variance at the output against ``int
    modal_spectrum['orbital'] df`` (one-sided) over the harmonics (broadening conserves
    line power; calibrated on the white fixture: +1.1e-3, and
    `orbital_correlation` -8.5e-4).  Measured on the Lorentzian van der Pol
    (gear): -3.1e-3 at 200 points, -7.8e-4 at 400 (gear's h^2).  The
    remaining refusals name their alternatives."""
    import warnings as _w
    from scipy.integrate import trapezoid
    _c, pss, pac, ov = _coloured_vdp('coloured', npts=200)
    f0 = 1.0 / float(pss.period)
    with pytest.raises(TypeError, match='COLOURED'):
        pac.oscillator_covariance(pss)
    with pytest.warns(RuntimeWarning, match="WHITE sources' alone"):
        K_orb, info = pac.oscillator_covariance(pss, samples=True,
                                                   colour_fmin=1e-6 * f0)
    d = info['d']
    assert d == 0.0 and not np.any(K_orb)
    tr = np.asarray(info['transverse_samples'])
    np.testing.assert_allclose(info['K_transverse'], tr[0], rtol=1e-12, atol=0)
    np.testing.assert_allclose(info['K_coloured'], tr[0], rtol=1e-12, atol=0)
    N = len(pss.factored_period().steps)
    cyc = float(np.mean(tr[:N, ov, ov]))
    tot = 0.0
    dl = np.geomspace(1e-6, 0.5, 150)
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        for j in range(1, 5):
            offs = (np.unique(np.concatenate((-(1.0 - dl[::-1]), -dl[::-1], dl)))
                    if j == 1 else np.concatenate((-dl[::-1], dl))) * f0
            ms = pac.modal_spectrum(pss, offs, ov, harmonic=j, maxharmonics=8, maxsidebands=16)
            tot += trapezoid(ms['orbital'], offs)
    assert abs(cyc / tot - 1.0) < 5e-3, cyc / tot - 1.0
    with pytest.raises(NotImplementedError, match='residue sum'):
        pac.orbital_correlation(pss)
    with pytest.raises(NotImplementedError, match='no diffusion constant'):
        pac.diffusion_constant(pss)


def test_the_white_only_covariance_routines_refuse_colour_with_the_reason():
    """⚠ A ROUTINE THAT FOLDS `CY` AT ONE FREQUENCY MUST SAY SO.

    `oscillator_covariance` (through `_lyapunov_pieces`) and
    `orbital_correlation` read `CY` at `2 pi / T` as if it held at every
    frequency.  (`oscillator_covariance` integrates the colour over a band
    since 2026-09-26 and, without one, asks for it -- a TypeError, a
    required argument, since 2026-09-29; `orbital_correlation` still
    refuses colour outright.)  (`covariance` did too until 2026-09-25; it now integrates
    the colour over a band and refuses only without one --
    `test_a_coloured_covariance_integrates_the_band_against_the_closed_form`.)  On a coloured source that returns a plausible
    number -- the shape A4d itself names -- so they refuse, and the
    refusal is the same test `_refuse_coloured` applies everywhere: `CY`
    at `w0` against `CY` at `10 w0`.
    """
    _c, pss, pac, _ov = _coloured_vdp('coloured', npts=240)
    for name, exc, call in (
            ('oscillator_covariance', TypeError,
             lambda: pac.oscillator_covariance(pss)),
            ('orbital_correlation', NotImplementedError,
             lambda: pac.orbital_correlation(pss)),
    ):
        with pytest.raises(exc, match='COLOURED'):
            call()


@pytest.mark.slow
def test_the_orbital_mode_basis_is_complete_only_for_noise_in_the_slow_subspace():
    """⚠⚠ A9's modal basis OMITS the annihilated modes, and they can carry ~all of it.

    `orbital_mode_weights` resolves `K_orb` onto the NON-NULL Floquet
    directions, so `Σ cw[k,k'] u_k u_k'^H` reproduces only the part of the
    covariance living on them.  The sibling test asserts that residual is
    `< 1e-2` and reads it as "the retained modes account for the covariance".

    **That holds because its fixture injects at the OSCILLATOR node.**  Move
    one current source and nothing else, on `_osc_with_ladder`'s circuit at
    `nslow = 4`::

        injected at            ||K_orb||    reconstruction residual
        the oscillator node    2.70e-05     1.80e-03   (0.18%)
        a SLOW ladder node     3.94e-01     3.56e-01   (36%)
        a FAST ladder node     6.87e+02     9.996e-01  (99.96%)
        a faster one           3.33e+03     9.999e-01  (99.99%)

    The annihilated modes are killed by the period map, so they reach the
    stationary covariance only through the `j = 0` term — but that term is not
    small when the noise is injected there, **and that is where device noise
    actually is**: every resistor in a bias or tuning network.  So a modal
    orbital spectrum built on this basis is complete only for noise entering
    the slow subspace, which is the minority case rather than the normal one.

    ⚠ AND IT MAKES THE RECONSTRUCTION RESIDUAL A DETECTOR, NOT A BOUND.  It
    catches a dropped NON-NULL mode well — which is what `orbital_mode_weights`
    claims for it — but it SATURATES at the floor the null modes carry, so it
    cannot certify a truncation below that floor however many modes are kept.

    ⚠ Independently reproduced by a peer session on a different oscillator with
    a different `K_orb` construction (69% there).  The mechanism transfers; the
    magnitude does not, and neither denominator counts the same modes.

    This test exists so the sibling's `rel < 1e-2` is never read as a property
    of the method.  It is a property of where that fixture puts its noise.
    """
    import warnings as _w
    from pycircuit.circuit.elements import IS as _IS
    circuit.default_toolkit = circuit.numeric

    def build(noise_node, nslow=4, nladder=14, Q=16.0, npts=200):
        tper = 2.0 * np.pi
        mu = 1.0 / (2.0 * np.pi * Q)
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0))
        prev = 'v'
        for j in range(nladder):
            nd = 'p%d' % j
            cir.add_node(nd)
            tau = (tper * 10.0 ** (-1.0 + 2.0 * j / max(nslow - 1, 1))
                   if j < nslow else tper * 1e-4)
            cir['r%d' % j] = R(prev, nd, r=1e3)
            cir['c%d' % j] = C(nd, gnd, c=tau / 1e3)
            prev = nd
        cir['n'] = _IS(noise_node, gnd, i=0.0, noisePSD=1e-6)
        T = 2.0 * np.pi / np.sqrt(max(1.0 - mu ** 2 / 4.0, 1e-9))
        pss = PSS(cir, method='gear', reltol=1e-11)
        x0 = np.zeros(cir.n - 1)
        x0[0] = 2.0
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=200)
        assert pss.converged, 'noise at %s did not converge' % noise_node
        return cir, pss

    def floor_of(node):
        cir, pss = build(node)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            cw, _wi = PAC(cir).orbital_mode_weights(pss)
            modes, K = _wi['modes'], _wi['K']
        cw = np.asarray(cw)
        U = np.column_stack([m['u0'] for m in modes])
        rec = U @ cw @ U.conj().T
        return (float(np.linalg.norm(rec - K)) / float(np.linalg.norm(K)),
                float(np.linalg.norm(K)), U.shape[1])

    osc, nK_osc, nmodes = floor_of('v')
    slow, _nK_s, _m = floor_of('p0')
    fast, nK_f, _m = floor_of('p6')

    ## (1) THE SIBLING'S CASE, reproduced -- if this is not small the contrast
    ## below has nothing to contrast against.
    assert osc < 1e-2, \
        'injecting at the oscillator node used to leave only 1.8e-03 outside ' \
        'the non-null basis and now leaves %.3e; the sibling test\'s reading ' \
        'rests on this' % osc

    ## (2) ⚠ AND IT IS THE INJECTION POINT THAT DECIDES IT.  Two orders between
    ## the same circuit noised in two places.
    assert fast > 0.9, \
        'noise in a FAST ladder branch should leave ~all of K_orb outside the ' \
        'non-null basis (0.9996 on record) and leaves %.4f. If this has ' \
        'fallen, the annihilated modes have stopped carrying the covariance ' \
        'and A9\'s basis is more complete than recorded -- re-measure before ' \
        'relying on it.' % fast
    assert slow > 10.0 * osc, \
        'even a SLOW ladder node should be far worse than the oscillator ' \
        'node (0.356 against 0.0018); got %.4f against %.4f' % (slow, osc)
    assert fast / osc > 100.0, \
        'the whole finding is the SPREAD across injection points: %.4f vs ' \
        '%.4f is only %.1fx' % (fast, osc, fast / osc)

    ## (3) VACUITY GUARD: the basis must actually be a truncation here, or a
    ## large residual would mean something else entirely.
    fp = build('v')[1].factored_period()
    assert nmodes < fp.width, \
        'the mode basis (%d) is not a truncation of the map (%d), so a ' \
        'reconstruction residual cannot be about omitted modes' \
        % (nmodes, fp.width)
    assert nK_f > nK_osc, \
        'injecting into a small fast capacitor should give a much LARGER ' \
        'covariance (6.9e+02 against 2.7e-05); got %.3e against %.3e -- if ' \
        'not, the source is not landing where this test thinks' \
        % (nK_f, nK_osc)


def _hostile_oscillator(npts=400):
    """THE fixture that is neither half-wave symmetric nor unit-reactance.

    ⚠⚠ WHY THIS EXISTS.  Every oscillator fixture in this file was van der
    Pol with `c = L = 1`, half-wave symmetric, starting at `[2, 0]` where the
    adjoint is axis-aligned.  Three coincidences, and on 2026-09-07 they hid
    two defects in `floquet_modes` from every gate in this file -- including
    A9's three-way gate, which is a good gate and caught a different adjoint
    defect the same week:

      * `q^T p = 1` where the DAE conserves `q^T C p` -- invisible when `C` is
        the identity up to sign;
      * the replayed adjoint is `C^T q`, used as `q` -- invisible when the
        seed is axis-aligned so `C^T q` is parallel to `q`.

    Together they put the orbital covariance 81x LOW on an asymmetric orbit
    against a Monte Carlo, while every test stayed green.

    `a = 0.30`, `c = 4`, `L = 1/4`: measured half-wave asymmetry 0.100,
    `|lam2| = 0.969`, adjoint separation `|cos(q_1, q_2)| = 0.70`.  Both
    properties, and the self-checks below refuse a fixture that has lost
    either, so it cannot be quietly tuned back to the blind one.
    """
    cir, pss = _a9_vdp(cval=4.0, lval=0.25, a=0.30)
    if npts != 400:
        import warnings as _w
        T = float(pss.period)
        pss = PSS(cir, method='gear', reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / npts, x0=np.array([2.0, 0.0]),
                      maxiterations=400)
        assert pss.converged
    ## self-checks: the fixture must actually BE hostile
    W = np.delete(np.asarray(pss.waveform[1], dtype=float), pss.irefnode,
                  axis=0)[0]
    half = len(W) // 2
    asym = float(np.max(np.abs(W[:half] + W[half:2 * half]))) / float(
        np.max(np.abs(W)))
    assert asym > 0.05, \
        'the "hostile" fixture has half-wave asymmetry %.3f; it must be ' \
        'asymmetric or it cannot see the C^T q defect' % asym
    x0r = np.delete(np.asarray(pss.waveform[1], dtype=float)[:, 0],
                    pss.irefnode)
    Cm = np.asarray(pss._C_at(x0r), dtype=float)
    assert np.max(np.abs(np.abs(np.diag(Cm)) - 1.0)) > 0.5, \
        'the "hostile" fixture has unit reactances (diag C = %r); it must ' \
        'not, or it cannot see the q^T p normalisation defect' % (
            np.diag(Cm).tolist(),)
    return cir, pss


def test_the_three_way_orbital_gate_holds_on_the_hostile_fixture():
    """A9's three-way gate, on the fixture it was blind without.

    Same three routes as `test_orbital_correlation_is_gated_three_ways` --
    eq (22), the definition integral, and the obliquely-projected Lyapunov
    cycle-mean -- on `_hostile_oscillator`.  On the symmetric unit-reactance
    fixture that gate passed at 0.3 % THROUGH two defects that put the
    answer 81x off elsewhere.  Here, with the fixes in, measured:

        npts   relAB (eq22 vs definition)   relAC (eq22 vs Lyapunov)
         400        9.1e-05                     1.6e-02
         800        4.6e-05                     8.1e-03

    `relAC` halves per doubling: the O(h) residual of the adjoint replay on
    an asymmetric orbit (`doc/shooting_history.md`, `_warn_if_orbit_is_asymmetric`).  The
    bound is 2x the 400-point measurement.  ⚠ Reverting the `C^-T` transform
    in `floquet_modes` takes `relAC` here to 2.98e-01 and this test RED
    (MEASURED, by doing exactly that) -- while the symmetric gate stays GREEN
    under the same mutation, blind -- which
    is the whole point: a fixture on which the defect is visible.
    """
    import warnings as _w
    ## ⚠ 2026-09-20: TWO grids, because the reference route is FIRST ORDER
    ## and this gate at one grid measured that, not eq (22).  See the ladder
    ## in the docstring: with the second-order adjoint the eq (22) route
    ## self-converges at >= 2nd order (2.7e-03 / 2.7e-04 / 1.5e-05 against
    ## N = 3200) while the Lyapunov route halves per doubling (7.8e-02 /
    ## 3.9e-02 / 1.7e-02); their distance at any one N is the reference's
    ## error.  The old first-order `q` read 1.6e-02 here by CANCELLATION
    ## (5.0e-02 / 2.5e-02 / 1.1e-02 self-convergence, and 2.0e-03 against the
    ## reference at 3200 -- closer than the exact adjoint, which cannot be
    ## right).  Pinned: the value at 400 (4.66e-02) and its halving from 200
    ## (9.41e-02), which is the reference's order.
    def gate(npts):
        cir, pss = _hostile_oscillator(npts=npts)
        pac = PAC(cir)
        m = cir.n - 1
        Tp = float(pss.period)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            CY2 = 0.5 * np.real(np.asarray(pac._cy_reduced(pss, 0.0)))
            R, _ = pac.orbital_correlation(pss, maxharmonics=8)
            modes = pss.floquet_modes(pss)
            v0, info = pss.ppv()
            Kf, ci = pac.oscillator_covariance(pss, samples=True)
            d = ci['d']
        return cir, pss, pac, m, Tp, CY2, R, modes, v0, info, Kf, d, ci
    cir, pss, pac, m, Tp, CY2, R, modes, v0, info, Kf, d, ci = gate(400)

    ## route B: the definition integral
    md = [x for x in modes if abs(abs(x['lam']) - 1.0) > 1e-6][0]
    mu2 = float(np.real(md['mu']))
    P2 = np.real(md['p'][:, :-1]); Q2 = np.real(md['q'][:, :-1])
    Nn = P2.shape[1]; h = Tp / Nn
    g = np.array([Q2[:, k] @ CY2 @ Q2[:, k] for k in range(Nn)])
    nper = max(int(np.ceil(-40.0 / (2 * mu2 * Tp))), 1)
    taus = np.arange(0, nper * Nn) * h
    wts = np.exp(2 * mu2 * taus) * h
    s2 = np.array([float(np.sum(wts * g[(k - np.arange(0, nper * Nn)) % Nn]))
                   for k in range(Nn)])
    Rdef = np.mean(np.stack([s2[k] * np.outer(P2[:, k], P2[:, k])
                             for k in range(Nn)]), axis=0)
    relAB = np.linalg.norm(R - Rdef) / np.linalg.norm(Rdef)
    assert relAB < 1e-3, 'eq (22) vs the definition integral: %.3e' % relAB

    ## route C: obliquely-projected Lyapunov cycle-mean
    Ps = [np.asarray(x, float)[:m, :m] for x in ci['samples']]
    G = [np.asarray(x, float)[:m, :m] for x in ci['growth_samples']]
    ## `samples[j]` is node j: prepending `v0` shifted every node by one
    ## (the 'first-order Lyapunov route' of 2026-09-20 was this shift)
    vs = [np.asarray(sv, float)[:m] for sv in info['samples']]
    proj = []
    for j in range(min(len(Ps), len(vs))):
        w, U = np.linalg.eigh(G[j])
        uj = U[:, np.argmax(w)] * np.sqrt(max(float(w.max()), 0.0))
        den = float(vs[j] @ uj)
        if abs(den) < 1e-300:
            continue
        Pi = np.eye(m) - np.outer(uj, vs[j]) / den
        proj.append(Pi @ Ps[j] @ Pi.T)
    Pm = np.mean(np.stack(proj), axis=0)
    relAC = np.linalg.norm(R - Pm) / np.linalg.norm(Pm)
    ## ⚠ 2026-09-20, second correction: the 4.66e-02 pinned here earlier was
    ## THIS TEST's own one-node shift of the phase-vector list ([v0] +
    ## samples), not the reference route.  Unshifted, the two routes agree
    ## at the 1e-3 level and the residual no longer halves with the grid.
    assert relAC < 2e-3, \
        'eq (22) disagrees with the Lyapunov reference by %.3e on the ' \
        'hostile fixture at 400 points (measured 7.9e-04 with the second-' \
        'order adjoint and the phase-vector list unshifted)' % relAC
    ## and the disagreement is the REFERENCE's: it halves with the grid
    _c2, _p2, _pac2, _m2, _T2, _CY2, R2, _mo2, v02, info2, _K2, _d2, ci2 = gate(200)
    Ps2 = [np.asarray(x, float)[:m, :m] for x in ci2['samples']]
    G2 = [np.asarray(x, float)[:m, :m] for x in ci2['growth_samples']]
    vs2 = [np.asarray(sv, float)[:m] for sv in info2['samples']]
    proj2 = []
    for j in range(min(len(Ps2), len(vs2))):
        w, U = np.linalg.eigh(G2[j])
        uj = U[:, np.argmax(w)] * np.sqrt(max(float(w.max()), 0.0))
        den = float(vs2[j] @ uj)
        if abs(den) < 1e-300:
            continue
        Pi = np.eye(m) - np.outer(uj, vs2[j]) / den
        proj2.append(Pi @ Ps2[j] @ Pi.T)
    Pm2 = np.mean(np.stack(proj2), axis=0)
    relAC200 = np.linalg.norm(R2 - Pm2) / np.linalg.norm(Pm2)
    ## SECOND order between the two routes once the list is unshifted:
    ## measured 3.07e-03 / 7.90e-04 / 2.12e-04 at 200 / 400 / 800 (3.9, 3.7)
    assert relAC200 < 6e-3, relAC200
    assert 3.0 < relAC200 / relAC < 4.8, (relAC200, relAC)


def test_every_period_harmonic_is_a_fourier_integral_on_a_non_uniform_grid():
    """⚠ THE SITES THE pnoise FIX LEFT ON AN INDEX DFT, measured then converted.

    Against a KNOWN answer -- a driven RC, first harmonic `H/(2j)` -- on a 3:1
    grid, before: `carrier_phasor` and `fpss` 7.5 % off and NOT converging
    (the same flat error as the folds had).  After, the trapezoid-weighted sum
    at the true times: second order, like the gear samples it reads.  The
    uniform grid keeps the index DFT bit for bit.

    And one thing the measurement found that is not a quadrature: the phase
    multiplier of an oscillator is 1 only while the DISCRETISATION is
    time-translation invariant.  A non-uniform grid breaks that at O(h^2)
    (1 - 5.1e-05 at N = 400), it fell outside `ppv`'s 1e-6 window, and `ppv`
    then reported the PHASE multiplier as `second_multiplier` -- silently,
    f_amp 600x too small.  `modal_spectrum` refuses such a grid and says why
    (converted anyway: with the mode forced it closes on pnoise to 1.0001 at
    N = 1600, 8-13 % off before).
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    w = 2 * np.pi / T
    Rv, Cv = 1e3, 0.3e-9
    truth = 1.0 / (1.0 + 1j * w * Rv * Cv) / 2j

    def fracs(n):
        f = 1.0 + 0.5 * np.sin(2 * np.pi * np.arange(n) / n)
        return f / f.sum()

    def run(n, nonuniform):
        c = SubCircuit()
        c.add_node('in')
        c.add_node('out')
        c['V'] = VSin('in', gnd, va=1.0, freq=1.0 / T, phase=0.0)
        c['R'] = R('in', 'out', r=Rv)
        c['C'] = C('out', gnd, c=Cv)
        pss = PSS(c, method='gear', reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res = pss.solve(period=T, timestep=T / n, maxiterations=60,
                            break_events=False,
                            grid=(fracs(n) if nonuniform else None))
        assert pss.converged
        cp = PAC(c, toolkit=circuit.numeric).carrier_phasor(pss, c.get_node_index('out'), 1)
        fp1 = complex(res['fpss'].v('out')[1])
        f1 = float(res['fpss'].sweep_values[1])
        return cp, fp1, f1

    errs = []
    for n in (200, 400):
        cp, fp1, f1 = run(n, True)
        assert abs(f1 * T - 1.0) < 1e-12
        ## fpss is RMS-folded: sqrt(2) times the coefficient
        assert abs(fp1 - np.sqrt(2) * cp) < 1e-12 * abs(cp)
        errs.append(abs(cp - truth) / abs(truth))
    assert errs[0] < 2e-3 and errs[1] < 5e-4, errs
    assert errs[0] / errs[1] > 3.0, errs
    ## the uniform grid is the index DFT, and as accurate as before
    cpu, fpu, _f = run(200, False)
    assert abs(cpu - truth) / abs(truth) < 2e-3
    assert abs(fpu - np.sqrt(2) * cpu) < 1e-12 * abs(cpu)

    ## the phase multiplier off the unit circle, and what ppv makes of it
    def vdp(nonuniform):
        mu = 1.0 / (2.0 * np.pi * 8.0)
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=4.0)
        cir['L'] = L('v', gnd, L=0.25)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        Tv = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
        pss = PSS(cir, method='gear', reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=Tv, timestep=Tv / 200, x0=np.array([2.0, 0.0]),
                      maxiterations=300, break_events=False,
                      grid=(fracs(200) if nonuniform else None))
            assert pss.converged
            lam2 = float(pss.ppv()[1]['second_multiplier'])
        return cir, pss, lam2

    _c, _p, lam_u = vdp(False)
    cir, pss, lam_n = vdp(True)
    phase = min(abs(abs(m['lam']) - 1.0) for m in pss.floquet_modes(pss))
    assert 1e-5 < phase < 1e-3, phase          # the premise: outside the window
    assert abs(lam_n / lam_u - 1.0) < 1e-3, (lam_n, lam_u)
    ## 2026-09-20: this used to pin a REFUSAL on the non-uniform gear grid;
    ## the phase mode is now found by its tangent alignment and the solve RUNS,
    ## warning with the measured departure (see `_phase_mode_split`)
    with _w.catch_warnings(record=True) as _rec:
        _w.simplefilter('always')
        PAC(cir, toolkit=circuit.numeric).modal_spectrum(
            pss, np.array([1e-3]), 0, maxharmonics=8)
    assert any('off the unit circle' in str(r.message) for r in _rec)


def test_floquet_modes_under_gear_are_second_order_on_a_uniform_grid_and_radau_is_exact_on_a_non_uniform_one():
    """⚠ GEAR'S ADJOINT MODES WERE FIRST ORDER; NOW SECOND, ON A UNIFORM GRID.

    Found as an 'unexplained 0.2 %' of modal_spectrum's parts on a 3:1 grid,
    it was gear's on any grid: `floquet_modes` reconstructed `q` from the
    discrete adjoint PAIR's first block `w1 = Jf^T t` through pinv(C^T) --
    `a0 * q(t + 2h/3)`, staggered by a fraction of a step -- so the invariant
    `q^T C p` drifted along the orbit and HALVED per doubling.  The per-step
    transposed solve `t` is the adjoint at the NEXT node; with `q_j = a0 *
    ts[j-1] exp(mu t_j)` the spread QUARTERS (measured 9.3e-04 / 2.3e-04 /
    5.6e-05 at N = 200 / 400 / 800, ratios 4.09 / 4.05 and 3.94 / 3.98).

    ⚠ ON A NON-UNIFORM GRID IT STAYS FIRST ORDER, at 2x the old accuracy:
    a variable-step multistep discrete adjoint draws a1, a2 from later steps
    and is the continuous adjoint's only to O(h) (four scalings measured, all
    halving).  Gear is refused by the modal spectra there anyway: its phase
    multiplier leaves the unit circle at O(h^2) (1 - 5e-05 at N = 400), as
    did trap's under its TR-BDF2 twin; radau's collocation solve keeps it at
    1 + 1e-11 and its modes reproduce the uniform grid to 1e-10, so
    modal_spectrum RUNS under radau on the 3:1 grid -- and under trap since
    its twin is radau by default (2026-09-24: 6e-12 off the circle).
    ⚠ Corrects b6e874a's 'the phase multiplier leaves 1 on a non-uniform
    grid': true of gear, not of radau.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)

    def fracs(n):
        f = 1.0 + 0.5 * np.sin(2 * np.pi * np.arange(n) / n)
        return f / f.sum()

    def solve(n, nonuniform, method):
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=4.0)
        cir['L'] = L('v', gnd, L=0.25)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
        pss = PSS(cir, method=method, reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / n, x0=np.array([2.0, 0.0]),
                      maxiterations=300, break_events=False,
                      grid=(fracs(n) if nonuniform else None))
        assert pss.converged
        return cir, pss

    def phase_off(pss):
        return min(abs(abs(m['lam']) - 1.0) for m in pss.floquet_modes(pss))

    def parts(cir, pss):
        pac = PAC(cir, toolkit=circuit.numeric)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            ms = pac.modal_spectrum(pss, np.array([1e-3, 3e-3]), 0, maxharmonics=8,
                                    maxsidebands=16)
        return np.concatenate([ms[k] for k in ('phase', 'orbital', 'correlation')])

    ## radau: on the unit circle on the 3:1 grid, runs, and agrees with the
    ## uniform grid to 2.6e-6 -- ORDER 5 AT 400 POINTS ON A GENUINE 3:1 GRID.
    ## ⚠ This used to be pinned at rtol 1e-8 (2026-09-20, earlier the same
    ## day) and read as "exact": the factored period then replayed on a
    ## UNIFORM grid whatever the solve used, so the two sides were the same
    ## replay.  Since the replay honours the solved grid the difference is
    ## radau's own fifth-order error on the coarse 3:1 steps.
    cu, pu = solve(400, False, 'radau')
    cn, pn = solve(400, True, 'radau')
    assert phase_off(pn) < 1e-9, phase_off(pn)
    np.testing.assert_allclose(parts(cn, pn), parts(cu, pu), rtol=1e-5, atol=0)
    ## gear and trap: O(h^2) off the circle on the same grid -- and they RUN
    ## (2026-09-20, Andreas: gear as a first-class choice on non-uniform
    ## grids): the phase mode is identified by its tangent alignment, its
    ## exponent forced to 0, and the departure is WARNED with its size.
    ## Gear's total converges to radau's at second order there: 7.6e-02 /
    ## 1.76e-02 / 4.2e-03 at N = 200 / 400 / 800 (ratios 4.35, 4.15).
    ## trap reads its twin, radau by default since 2026-09-24: on the circle
    ## as radau is (6e-12; under the TR-BDF2 twin it read 6.0e-6)
    _ct, pt = solve(400, True, 'trap')
    assert phase_off(pt) < 1e-9, phase_off(pt)
    for method in ('gear',):
        cir, pss = solve(400, True, method)
        off = phase_off(pss)
        ## ⚠ the lower bound was 1e-5 while the replay was UNIFORM whatever
        ## the solve used; gear's replay was always on the caller's grid
        assert 1e-6 < off < 1e-4, (method, off)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            PAC(cir, toolkit=circuit.numeric).modal_spectrum(
                pss, np.array([1e-3]), 0, maxharmonics=8)
        assert any('off the unit circle' in str(r.message) for r in rec), method
    g200 = parts(*solve(200, True, 'gear'))
    g400 = parts(*solve(400, True, 'gear'))
    ref = parts(cu, pu)                       # radau: uniform == 3:1 to 3e-6
    e200 = float(np.max(np.abs(g200 / ref - 1.0)))
    e400 = float(np.max(np.abs(g400 / ref - 1.0)))
    assert e400 < 4e-2, e400
    assert 2.5 < e200 / e400 < 5.5, (e200, e400)       # second order

    ## gear's invariant on a UNIFORM grid: second order now (it halved before)
    def drift(n):
        cir, pss = solve(n, False, 'gear')
        W = np.delete(np.asarray(pss.waveform[1], dtype=float), pss.irefnode, axis=0)
        out = []
        for md in pss.floquet_modes(pss):
            P, Q = np.asarray(md['p']), np.asarray(md['q'])
            K = min(P.shape[1], Q.shape[1], W.shape[1])   # the modes' own width
            a = np.abs([np.vdot(Q[:, j], np.asarray(pss._C_at(W[:, j]), dtype=float)
                                @ P[:, j]) for j in range(K)])
            out.append((a.max() - a.min()) / a.mean())
        return max(out)
    d200, d400 = drift(200), drift(400)
    assert 3e-4 < d200 < 3e-3, d200                   # was 9.9e-03 before the fix
    assert 3.3 < d200 / d400 < 4.8, (d200, d400)      # second order, not first


def test_gear_adjoint_modes_are_second_order_on_a_non_uniform_grid_and_orbital_correlation_refuses_an_off_circle_phase_mode():
    """The separately-discretised continuous adjoint (Andreas: "Do it"), gated
    on STRUCTURE: gear's two-step transpose on a grid whose step changes is a
    first-order scheme for the adjoint equation (it draws a1, a2 from later
    steps) and no rescaling lifts it -- four measured, all halving.  On such a
    grid `floquet_modes` now integrates `C^T dq/dt = G^T q` backwards with BDF2
    on the REVERSE grid's own step pair, matched to the forward multiplier.
    Invariant `q^T C p` spread at N = 200 / 400 / 800, 3:1 grid, measured:

        transpose (a0-scaled)   2.2e-02  1.1e-02  5.9e-03   halving
        continuous              1.8e-03  4.6e-04  1.1e-04   QUARTERING

    Uniform grids and one-step kinds never take this path (bit-identical).

    ⚠ And a latent defect the gate found: `orbital_correlation` swept an
    off-circle PHASE mode (gear's sits at 1 - 5e-05 on that grid) into its
    orbital sum and returned 146x .. 2449x radau's R, growing as N^2,
    silently.  It refuses now, naming radau, as `modal_spectrum` does.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)

    def fracs(n):
        f = 1.0 + 0.5 * np.sin(2 * np.pi * np.arange(n) / n)
        return f / f.sum()

    def solve(n, nonuniform, method):
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=4.0)
        cir['L'] = L('v', gnd, L=0.25)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
        pss = PSS(cir, method=method, reltol=1e-12)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / n, x0=np.array([2.0, 0.0]),
                      maxiterations=300, break_events=False,
                      grid=(fracs(n) if nonuniform else None))
        assert pss.converged
        return cir, pss

    def drift(pss):
        W = np.delete(np.asarray(pss.waveform[1], dtype=float), pss.irefnode, axis=0)
        out = []
        for md in pss.floquet_modes(pss):
            P, Q = np.asarray(md['p']), np.asarray(md['q'])
            K = min(P.shape[1], Q.shape[1], W.shape[1])
            a = np.abs([np.vdot(Q[:, j], np.asarray(pss._C_at(W[:, j]), dtype=float)
                                @ P[:, j]) for j in range(K)])
            out.append((a.max() - a.min()) / a.mean())
        return max(out)

    ## gear on the 3:1 grid: the continuous adjoint, second order
    d200 = drift(solve(200, True, 'gear')[1])
    d400 = drift(solve(400, True, 'gear')[1])
    assert 8e-4 < d200 < 4e-3, d200                    # 2.2e-02 on the transpose
    assert 3.3 < d200 / d400 < 4.8, (d200, d400)

    ## gear, trap and radau all RUN on the 3:1 grid now (the phase mode by its
    ## tangent alignment, never in the orbital sum): gear's R converges to
    ## radau's at ~x3.5 per doubling (2.3e-02 / 7.2e-03 / 2.0e-03 measured)
    Rr, _c = PAC(*[solve(400, True, 'radau')[0]], toolkit=circuit.numeric).orbital_correlation(
        solve(400, True, 'radau')[1], maxharmonics=8)
    for method in ('gear', 'trap'):
        cir, pss = solve(400, True, method)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            R, _c = PAC(cir, toolkit=circuit.numeric).orbital_correlation(pss, maxharmonics=8)
        assert np.all(np.isfinite(R))
        assert np.linalg.norm(R - Rr) / np.linalg.norm(Rr) < 3e-2, method


def _vdp_tank_noise(tau_over_T, npts=200):
    """van der Pol (Q = 8, gear -- a PAIR map) with one current source on
    the tank: white (`tau_over_T` None) or a Lorentzian `IS(noiseTau)` of
    the same level."""
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: mu * (u - u ** 3 / 3.0))
    kw = {'noiseTau': tau_over_T * T} if tau_over_T else {}
    c['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6, **kw)
    pss = PSS(c, method='gear', reltol=1e-12)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    with _w.catch_warnings():
        _w.simplefilter('ignore')
        pss.solve(period=T, timestep=T / npts, x0=x0, maxiterations=300)
    assert pss.converged
    return c, pss


def test_the_orbital_mode_weights_take_a_coloured_source_in_the_maps_own_space():
    """#17 B3 (2026-09-29): `orbital_mode_weights(colour_fmin=...)`.  The
    coloured part of the bounded covariance is built in the MAP's own space
    (the pair `(x_n, x_{n-1})` on gear) from the bordered solution at node 0,
    so the phase row of the coloured weights vanishes BY THE BORDER
    (`v_pair' wb = 0`): measured 1.1e-12 of their scale -- a stacking of
    node-projected responses leaks it at first order.  And a Lorentzian whose
    corner is far above the band (tau = 1e-4 T) is the white source of the
    same level: the orbital weight agrees with the white Lyapunov route's to
    1.75e-6 (measured), two routes sharing no integral."""
    cw_w, pw = _vdp_tank_noise(None)
    white, _iw = PAC(cw_w, toolkit=circuit.numeric).orbital_mode_weights(pw)
    cc, pc = _vdp_tank_noise(1e-4)
    f0 = 1.0 / float(pc.period)
    pac = PAC(cc, toolkit=circuit.numeric)
    with pytest.raises(TypeError, match='COLOURED'):
        pac.orbital_mode_weights(pc)
    cw, info = pac.orbital_mode_weights(pc, colour_fmin=1e-6 * f0)
    cwc = info['cw_coloured']
    assert np.array_equal(cw, cwc)           # every source coloured here
    kph, _orb = pac._phase_mode_split(pc, info['modes'], 'test')
    scale = np.max(np.abs(cwc))
    assert np.max(np.abs(cwc[kph, :])) < 1e-10 * scale
    assert np.max(np.abs(cwc[:, kph])) < 1e-10 * scale
    o = [k for k in range(cwc.shape[0]) if k != kph]
    rel = np.abs(cwc[np.ix_(o, o)] / white[np.ix_(o, o)] - 1.0)
    assert np.max(rel) < 1e-4, rel
    assert np.shape(info['K_coloured']) == np.shape(info['K']) ==         (2 * (cc.n - 1),) * 2            # the gear PAIR space


def test_the_phase_mode_is_found_with_a_dc_source_on_the_tank():
    """The review's D1 (2026-10-01): `_phase_mode_split` picks the phase mode
    by its alignment with the orbit tangent, ``C xdot = -(i + u)``, and left
    the SOURCE out.  A 1 A DC current on the tank node of a van der Pol
    (Q = 8) put `xdot(0)` 56 % off, the phase mode's alignment 0.83 under
    the 0.9 bar, and `orbital_correlation` REFUSED the oscillator ("no
    Floquet mode ... aligned with the orbit tangent").  With the source the
    alignment is 1.000 and the split goes through."""
    from pycircuit.circuit.tests._shooting_fixtures import _loss_osc
    _cir, pss, pac, _rs = _loss_osc('parallel', idc_node='v', idc=1.0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        modes = pss.floquet_modes()
        k, _rest = pac._phase_mode_split(pss, modes, 'test')
        R, _C = pac.orbital_correlation(pss)
    assert abs(abs(complex(modes[k]['lam'])) - 1.0) < 1e-9
    assert np.all(np.isfinite(R))
