
/* -- the hdl limiter: _limiting.py's laws in Python's own forms ----------- */
static double _pymax(double a, double b) { return (b > a) ? b : a; }
static double _pymin(double a, double b) { return (b < a) ? b : a; }

static double _lim_pnj(double vnew, double vold, double VT, double IS) {
  if (IS <= 0.0) return vnew;
  double vc = VT * log(VT / (IS * 1.414213562));
  if (fabs(vnew - vold) <= 2.0 * VT) return vnew;
  if (vnew > vc && vnew > 0.0) {
    if (vold > 0.0) {
      double arg = 1.0 + (vnew - vold) / VT;
      return (arg > 0.0) ? vold + VT * log(arg) : vc;
    }
    return VT * log(vnew / VT);
  }
  return vnew;
}

static double _lim_fet(double vnew, double vold, double vto) {
  double vtsthi = fabs(2.0 * (vold - vto)) + 2.0;
  double vtstlo = fabs(vold - vto) + 1.0;
  double vtox = vto + 3.5;
  double delv = vnew - vold;
  if (vold >= vto) {
    if (vold >= vtox) {
      if (delv <= 0.0) {
        if (vnew >= vtox) {
          if (-delv > vtstlo) vnew = vold - vtstlo;
        } else {
          vnew = _pymax(vnew, vto + 2.0);
        }
      } else {
        if (delv >= vtsthi) vnew = vold + vtsthi;
      }
    } else {
      if (delv <= 0.0) vnew = _pymax(vnew, vto - 0.5);
      else vnew = _pymin(vnew, vto + 4.0);
    }
  } else {
    if (delv <= 0.0) {
      if (-delv > vtsthi) vnew = vold - vtsthi;
    } else {
      double vtemp = vto + 0.5;
      if (vnew <= vtemp) {
        if (delv > vtstlo) vnew = vold + vtstlo;
      } else {
        vnew = vtemp;
      }
    }
  }
  return vnew;
}

static double _lim_vds(double vnew, double vold) {
  if (vold < 0.0) return -_lim_vds(-vnew, -vold);
  if (vold >= 3.5) {
    if (vnew > vold) vnew = _pymin(vnew, 3.0 * vold + 2.0);
    else if (vnew < 3.5) vnew = _pymax(vnew, 2.0);
  } else {
    if (vnew > vold) vnew = _pymin(vnew, 4.0);
    else vnew = _pymax(vnew, -0.5);
  }
  return vnew;
}

static double _lim_delta(double vnew, double vold, double vmax) {
  double d = vnew - vold;
  if (d > vmax) return vold + vmax;
  if (d < -vmax) return vold - vmax;
  return vnew;
}

/* `apply_limit`: pnj's pars are (IS, VT), the law takes (VT, IS) */
static double _lim(int kind, double vnew, double vold, double p0, double p1) {
  switch (kind) {
    case 0: return _lim_pnj(vnew, vold, p1, p0);
    case 1: return _lim_fet(vnew, vold, p0);
    case 2: return _lim_vds(vnew, vold);
    case 3: return _lim_delta(vnew, vold, p0);
    default: return vnew;
  }
}

/* Python's tuple `<` on (double, int, int) keys: the first element that
   is not `==` decides */
static int _tless(double a0, int a1, int a2, double b0, int b1, int b2) {
  if (!(a0 == b0)) return a0 < b0;
  if (a1 != b1) return a1 < b1;
  return a2 < b2;
}

/* `_tless` on the keys of order[i] and order[j] */
static int _klt(const int *o, int i, int j, const double *k0, const int *k1,
                const int *k2) {
  return _tless(k0[o[i]], k1[o[i]], k2[o[i]], k0[o[j]], k1[o[j]], k2[o[j]]);
}

static void _rev(int *o, int lo, int hi) {
  for (hi--; lo < hi; lo++, hi--) { int t = o[lo]; o[lo] = o[hi]; o[hi] = t; }
}

/* `sorted(order, key=(k0[j], k1[j], k2[j]))` as CPython 3.14's list.sort
   orders a list under 64 entries -- `count_run` (the ascending prefix, else
   a descending run with its equal stretches kept in order, reversed, then
   extended), then binary insertion of the rest -- on `_tless`'s compares
   (the keys are indexed by the entries of `order`).  A key holding a NaN,
   never `==` nor `<`, lands where Python's sort leaves it; keys without
   one, where any stable sort does.  (`_hdl_climit.NAN_KEYS`: off CPython
   3.14 a NaN key is declined before it gets here.) */
static void _sort_by(int *order, int n, const double *k0, const int *k1,
                     const int *k2) {
  int m = 1;
  if (n < 2) return;
  while (m < n && !_klt(order, m, m - 1, k0, k1, k2)) m++;
  if (m < n && !(m > 1 && _klt(order, 0, m - 1, k0, k1, k2))) {
    int neq = 0;
    if (m > 1) _rev(order, 0, m);
    for (m++; m < n; m++) {
      if (_klt(order, m, m - 1, k0, k1, k2)) {
        if (neq) { neq++; _rev(order, m - neq, m); neq = 0; }
      } else if (_klt(order, m - 1, m, k0, k1, k2)) {
        break;
      } else {
        neq++;
      }
    }
    if (neq) { neq++; _rev(order, m - neq, m); }
    _rev(order, 0, m);
    while (m < n && !_klt(order, m, m - 1, k0, k1, k2)) m++;
  }
  for (; m < n; m++) {
    int L = 0, R = m, pivot = order[m];
    do {
      int M = (L + R) >> 1;
      if (_tless(k0[pivot], k1[pivot], k2[pivot], k0[order[M]], k1[order[M]], k2[order[M]]))
        R = M;
      else
        L = M + 1;
    } while (L < R);
    for (int M = m; M > L; M--) order[M] = order[M - 1];
    order[L] = pivot;
  }
}

static int _find(int *parent, int n) {
  while (parent[n] != n) { parent[n] = parent[parent[n]]; n = parent[n]; }
  return n;
}

/* `device_writeback(out, targets, drift, pinned=moved)`: the rows it
   wrote are marked in `moved` AFTER the walk (the walk reads the set as it
   was).  Returns -1 where a Kruskal key is NaN and NaN keys are declined
   (PCX_NAN_DECLINE: off CPython 3.14, `_hdl_climit.NAN_KEYS`). */
#define _WB_MAX 32
static int _writeback(double *o, const int *tra, const int *trb,
                      const double *tvn, const double *tvl, int nt,
                      const double *drift, int *moved, int nx) {
  int any = 0;
  for (int i = 0; i < nt; i++) if (tvl[i] != tvn[i]) any = 1;
  if (!any) return 0;
  double k0[_WB_MAX];
  int order[_WB_MAX];
  for (int i = 0; i < nt; i++) {
    k0[i] = -fabs(tvl[i] - tvn[i]);
#ifdef PCX_NAN_DECLINE
    if (isnan(k0[i])) return -1;
#endif
    order[i] = i;
  }
  _sort_by(order, nt, k0, tra, trb);
  int parent[_WB_MAX], innode[_WB_MAX], adjn[_WB_MAX];
  int adjto[_WB_MAX][_WB_MAX];
  double adjs[_WB_MAX][_WB_MAX];
  for (int k = 0; k < nx; k++) { parent[k] = k; innode[k] = 0; adjn[k] = 0; }
  for (int m = 0; m < nt; m++) {
    int i = order[m];
    int ra = tra[i], rb = trb[i];
    double vl = tvl[i];
    int ka = _find(parent, ra), kb = _find(parent, rb);
    if (ka == kb) continue;
    parent[ka] = kb;
    adjto[ra][adjn[ra]] = rb; adjs[ra][adjn[ra]] = -vl; adjn[ra]++;
    adjto[rb][adjn[rb]] = ra; adjs[rb][adjn[rb]] = +vl; adjn[rb]++;
    innode[ra] = 1; innode[rb] = 1;
  }
  int written[_WB_MAX], seen[_WB_MAX], done[_WB_MAX], stack[_WB_MAX];
  for (int k = 0; k < nx; k++) { written[k] = 0; seen[k] = 0; done[k] = 0; }
  /* the components, their members in increasing row order (the set's) */
  for (int r = 0; r < nx; r++) {
    if (!innode[r]) continue;
    int root = _find(parent, r);
    if (done[root]) continue;
    done[root] = 1;
    /* the anchor: min by (0 if pinned else 1, drift, row) */
    int anchor = -1;
    int a0 = 0; double a1 = 0.0;
    for (int n = r; n < nx; n++) {
      if (!innode[n] || _find(parent, n) != root) continue;
      int b0 = moved[n] ? 0 : 1;
      double b1 = drift[n];
      int less;
      if (anchor < 0) less = 1;
      else if (b0 != a0) less = b0 < a0;
      else if (!(b1 == a1)) less = b1 < a1;
      else less = n < anchor;
      if (less) { anchor = n; a0 = b0; a1 = b1; }
    }
    int sp = 0;
    seen[anchor] = 1;
    stack[sp++] = anchor;
    while (sp > 0) {
      int u = stack[--sp];
      for (int e = 0; e < adjn[u]; e++) {
        int v = adjto[u][e];
        double s = adjs[u][e];
        if (seen[v]) continue;
        seen[v] = 1;
        if (moved[v]) continue;
        o[v] = o[u] + s;
        written[v] = 1;
        stack[sp++] = v;
      }
    }
  }
  for (int k = 0; k < nx; k++) if (written[k]) moved[k] = 1;
  return 0;
}
