
/* -- PCNR's refine (`_pcnr_c`): `pcnr.limit_block`'s rule for a device whose
   laws over-determine its unknowns (a redundant probe: SPICE's MOSFET set,
   `(b,d) = (b,s) - (d,s)`).  Follows limit.c (its `_sort_by`, `_find`) in
   every class's refine kernel. ----------------------------------------- */
#define _PT_MAX 32

/* the probes kept by Kruskal on decreasing |lim - raw| (ties by row pair,
   `_sort_by`'s stable order), each joining two components: `kept`, their
   count returned -- or -1 where a key is NaN (the Python answers) */
static int _pcnr_kept(int np, const int *A, const int *B, const double *raw,
                      const double *lim, int nrows, int *kept) {
  double k0[_PT_MAX];
  int order[_PT_MAX], parent[_PT_MAX], nk = 0;
  for (int q = 0; q < np; q++) {
    k0[q] = -fabs(lim[q] - raw[q]);
    if (isnan(k0[q])) return -1;
    order[q] = q;
  }
  _sort_by(order, np, k0, A, B);
  for (int r = 0; r < nrows; r++) parent[r] = r;
  for (int i = 0; i < np; i++) {
    int q = order[i];
    int ka = _find(parent, A[q]), kb = _find(parent, B[q]);
    if (ka == kb) continue;
    parent[ka] = kb;
    kept[nk++] = q;
  }
  return nk;
}

/* the unknowns from the kept tree's potentials, as the Python forms them:
   the nodes in the order a kept probe first names them (its dict), each
   one's neighbours in insertion order, a depth-first walk on a stack from
   each unvisited node at 0, `pot[w] = pot[u] + sgn` (sgn -lim from a, +lim
   from b), then `pot[a] - pot[b]` per tree probe */
static void _pcnr_potentials(int m, const int *A, const int *B, const double *lim,
                             const int *kept, int nk, int nrows, double *out) {
  int nodes[_PT_MAX], present[_PT_MAX], deg[_PT_MAX], has[_PT_MAX];
  int adj_w[_PT_MAX][2 * _PT_MAX], stack[_PT_MAX];
  double adj_s[_PT_MAX][2 * _PT_MAX], pot[_PT_MAX];
  int nn = 0;
  for (int r = 0; r < nrows; r++) { present[r] = 0; deg[r] = 0; has[r] = 0; pot[r] = 0.0; }
  for (int i = 0; i < nk; i++) {
    int q = kept[i], a = A[q], b = B[q];
    if (!present[a]) { present[a] = 1; nodes[nn++] = a; }
    adj_w[a][deg[a]] = b; adj_s[a][deg[a]++] = -lim[q];
    if (!present[b]) { present[b] = 1; nodes[nn++] = b; }
    adj_w[b][deg[b]] = a; adj_s[b][deg[b]++] = +lim[q];
  }
  for (int s = 0; s < nn; s++) {
    int start = nodes[s], sp = 0;
    if (has[start]) continue;
    pot[start] = 0.0; has[start] = 1;
    stack[sp++] = start;
    while (sp) {
      int u = stack[--sp];
      for (int e = 0; e < deg[u]; e++) {
        int w = adj_w[u][e];
        if (!has[w]) { pot[w] = pot[u] + adj_s[u][e]; has[w] = 1; stack[sp++] = w; }
      }
    }
  }
  for (int j = 0; j < m; j++) out[j] = pot[A[j]] - pot[B[j]];
}

/* the over-determined case of `limit_block`: the kept set, the tree itself
   (`sorted(kept) == range(m)`: the limited values as they are) or the
   potentials.  0, or 1 where a key is NaN (declined) */
static int _pcnr_tree(int m, int np, const int *A, const int *B, const double *raw,
                      const double *lim, int nrows, double *out) {
  int kept[_PT_MAX];
  int nk = _pcnr_kept(np, A, B, raw, lim, nrows, kept);
  if (nk < 0) return 1;
  int tree = nk == m;
  for (int i = 0; tree && i < nk; i++) if (kept[i] >= m) tree = 0;
  if (tree) {
    for (int j = 0; j < m; j++) out[j] = lim[j];
    return 0;
  }
  _pcnr_potentials(m, A, B, lim, kept, nk, nrows, out);
  return 0;
}
