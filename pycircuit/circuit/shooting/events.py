"""The event columns of a staged solve (`EventColumns`): the total map, the
Schur elimination and the costate injections, shared by every consumer.
"""
import numpy as np


class EventColumns(dict):
    """A staged solve's EVENT COLUMNS and the linear algebra every bordered
    consumer does with them.

    The dict keys are the ones the consumers always read -- `nodes` (the
    grid nodes the crossings land on), `W`, `c` (the event rows `W_k . x =
    c_k`), `P_end` (the columns at the period: `d x_N / d theta`, m x K or
    gear's 2m x K), `Pk_nodes` (the columns at every node), `P_nodes` (the
    homogeneous maps to every node), `G = W P_node` and `Gt = W Pk_node`
    (the event rows' derivatives) -- so every `ev['P_end']` read keeps
    working and a consumer's "unbordered" control stays
    `pss._event_columns = None`.  `dth` is `dtheta/dx_0 = -Gt^-1 G`.

    The derivations every consumer needs -- the TOTAL map `M + P_end dth`,
    the bordered adjoint's Schur elimination, the `-zeta_k W_k` injection
    -- are methods here; a consumer calls them rather than re-deriving
    them.

    History: `doc/shooting_history.md`, `EventColumns`.
    """

    @classmethod
    def build(cls, nodes, W, c, P_end, Pk_nodes, P_nodes, G=None,
              lazy_P_nodes=None):
        """The columns from their node arrays: `G`, `Gt` and `dth`.

        ⚠ A MATRIX-FREE STAGE HAS NO DENSE MAP TO EVERY NODE: it passes
        `P_nodes` None with `G` (the event rows' derivatives, one reverse
        replay each) and a builder the dense `P_nodes` comes from on first
        READ (`__missing__`) -- only the consumers that need a node's whole
        map (PAC's driven bordered solve, the covariance closure) pay for
        it, once."""
        K = len(nodes)
        Pk_nodes = np.asarray(Pk_nodes, dtype=float)
        Gt = np.zeros((K, K))
        if P_nodes is not None:
            P_nodes = np.asarray(P_nodes, dtype=float)
            G = np.zeros((K, P_nodes.shape[2]))
        for k, jn in enumerate(nodes):
            if P_nodes is not None:
                G[k] = W[k] @ P_nodes[jn]
            Gt[k] = W[k] @ Pk_nodes[jn]
        if P_nodes is not None:
            ev = cls(nodes=list(nodes), W=W, c=c, P_end=P_end,
                     Pk_nodes=Pk_nodes, P_nodes=P_nodes, G=G, Gt=Gt)
        else:
            G = np.asarray(G, dtype=float)
            ev = cls(nodes=list(nodes), W=W, c=c, P_end=P_end,
                     Pk_nodes=Pk_nodes, G=G, Gt=Gt)
        ev._lazy_P_nodes = lazy_P_nodes
        ev.dth = -np.linalg.solve(Gt, G)
        return ev

    def __missing__(self, key):
        ## the dense maps to every node, built on first read on a
        ## matrix-free stage (see `build`)
        lazy = getattr(self, '_lazy_P_nodes', None)
        if key == 'P_nodes' and lazy is not None:
            P = np.asarray(lazy(), dtype=float)
            self['P_nodes'] = P
            return P
        raise KeyError(key)

    @classmethod
    def from_capture(cls, captured, nsteps, m, P0, nodes, W, c, P_end):
        """The columns from a traversal's `_captured` nodes (`(x, P, [Pk])`
        per node 1..nsteps); `P0` is the map to node 0 -- the identity, or
        gear's `[I, 0]` on its `(x_0, x_-1)` pair."""
        K = np.asarray(P_end).shape[1]
        Pk_nodes = np.zeros((nsteps + 1, m, K))
        P_nodes = np.zeros((nsteps + 1, m, np.asarray(P0).shape[1]))
        P_nodes[0] = P0
        for j in range(1, nsteps + 1):
            _xj, Pj, Pkj = captured[j]
            P_nodes[j] = np.asarray(Pj, dtype=float)
            for k in range(K):
                Pk_nodes[j, :, k] = np.asarray(Pkj[k], dtype=float).ravel()
        return cls.build(nodes, W, c, P_end, Pk_nodes, P_nodes)

    @classmethod
    def from_capture_rows(cls, captured, nsteps, m, nodes, W, c, P_end, G,
                          lazy_P_nodes):
        """The columns of a MATRIX-FREE stage: the event columns at every
        node from `_captured` (their dense maps are None), the event rows'
        derivatives `G` given, `P_nodes` built on first read (see
        `build`)."""
        K = np.asarray(P_end).shape[1]
        Pk_nodes = np.zeros((nsteps + 1, m, K))
        for j in range(1, nsteps + 1):
            _xj, _Pj, Pkj = captured[j]
            for k in range(K):
                Pk_nodes[j, :, k] = np.asarray(Pkj[k], dtype=float).ravel()
        return cls.build(nodes, W, c, P_end, Pk_nodes, None, G=G,
                         lazy_P_nodes=lazy_P_nodes)

    @staticmethod
    def of(pss, n=None):
        """The solve's columns, or `None` when it is not staged -- or, with
        `n`, when they were built on a map of another width (a host that
        borrows a twin, trbdf2's covariance on gear)."""
        ev = getattr(pss, '_event_columns', None)
        if ev is None or (n is not None and not ev.fits(n)):
            return None
        return ev

    def fits(self, n):
        return np.asarray(self['P_end']).shape[0] == n

    def total_matrix(self, M):
        """The TOTAL monodromy `M + P_end dth`: the crossings move with the
        state."""
        return M + np.asarray(self['P_end'], dtype=float) @ np.asarray(self.dth, dtype=float)

    def total_matvec(self, mv, transposed=False):
        """`mv` (the fixed-grid map's product) made the total map's."""
        P = np.asarray(self['P_end'], dtype=float)
        D = np.asarray(self.dth, dtype=float)
        if transposed:
            return lambda z: np.asarray(mv(z)) + D.T @ (P.T @ np.asarray(z))
        return lambda z: np.asarray(mv(z)) + P @ (D @ np.asarray(z))

    def injection(self, zeta, nsteps, width, dtype=complex):
        """The event rows' costate term as a reverse pass's `inject` array:
        `-zeta_k W_k` at node `nd_k`."""
        inject = np.zeros((nsteps, width), dtype=dtype)
        for k, nd in enumerate(self['nodes']):
            if 0 <= nd < nsteps:
                wk = np.asarray(self['W'][k])
                inject[nd, :len(wk)] -= zeta[k] * wk
        return inject

    def injection_dict(self, zeta):
        """The same as `{node: vector}`, the stage passes' `extra` form."""
        return {nd: -zeta[k] * np.asarray(self['W'][k])
                for k, nd in enumerate(self['nodes'])}

    def costate_injection(self, v, nsteps, width):
        """The injection that samples a LEFT vector `v` of the total map
        along the orbit: `zeta = Gt^-T P_end^T v` -- the transpose of the
        saltation, carried by the reverse pass to every earlier node."""
        v = np.asarray(v)
        zeta = np.linalg.solve(np.asarray(self['Gt']).T,
                               np.asarray(self['P_end']).T @ v)
        return self.injection(zeta, nsteps, width, v.dtype)

    def forced_shift(self, f_nodes):
        """The crossings' motion driven by the SOURCE alone, `-Gt^-1 W
        f_node` (the forced response at each event node)."""
        W = np.asarray(self['W'])
        r = np.array([W[k] @ f_nodes[nd] for k, nd in enumerate(self['nodes'])])
        return -np.linalg.solve(np.asarray(self['Gt'], dtype=complex), r)

    @staticmethod
    def g_theta(cn, Pk_fixed, d, N):
        """The output's theta-sensitivity: `sum_n c_n d . Pk_fixed[n]`."""
        return np.array([np.sum(cn * (Pk_fixed[:N, :, kk] @ d))
                         for kk in range(Pk_fixed.shape[2])])

    def collapsed_zeta(self, g_theta, alpha, z):
        """`zeta = Gt^-T (g_theta + a P_end^T z)` -- the event rows' costate
        when `z` already solved the TOTAL operator (an oscillator's
        deflated adjoint, the sampled series)."""
        return np.linalg.solve(np.asarray(self['Gt']).T,
                               g_theta + alpha * (np.asarray(self['P_end'], dtype=complex).T
                                                  @ np.asarray(z)))

    def bordered_adjoint(self, solve, g, g_theta, alpha):
        """`B^T [z; zeta] = [g; g_theta]` with `B = [[I - aM, -a P_end],
        [G, Gt]]`, by block elimination: `z = z_g - Z_G zeta`, `Z_G = (I -
        aM)^-T G^T`, `zeta = (Gt^T + a P_end^T Z_G)^-1 (g_theta + a
        P_end^T z_g)`.  `solve(b, k)` solves `(I - aM)^T x = b` (`k` the
        event index for its failure message, `None` for `g`).  The
        transpose of `PAC.solve`'s bordered forward system, and dual-
        consistent with it (the suite's adjoint tests)."""
        G = np.asarray(self['G'])
        K = G.shape[0]
        z_g = solve(g, None)
        Z_G = np.column_stack([solve(np.asarray(G[k], dtype=complex), k)
                               for k in range(K)])
        Pth = np.asarray(self['P_end'], dtype=complex)
        Sb = np.asarray(self['Gt']).T + alpha * (Pth.T @ Z_G)
        zeta = np.linalg.solve(Sb, g_theta + alpha * (Pth.T @ z_g))
        return z_g - Z_G @ zeta, zeta
