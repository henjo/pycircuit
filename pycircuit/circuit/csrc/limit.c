
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

/* `sorted(order, key=(k0[j], k1[j], k2[j]))` for finite k0: a stable
   insertion sort (the keys are indexed by the entries of `order`) */
static void _sort_by(int *order, int n, const double *k0, const int *k1,
                     const int *k2) {
  for (int i = 1; i < n; i++) {
    int key = order[i];
    int j = i - 1;
    while (j >= 0 && _tless(k0[key], k1[key], k2[key],
                            k0[order[j]], k1[order[j]], k2[order[j]])) {
      order[j + 1] = order[j];
      j--;
    }
    order[j + 1] = key;
  }
}

static int _find(int *parent, int n) {
  while (parent[n] != n) { parent[n] = parent[parent[n]]; n = parent[n]; }
  return n;
}

/* `device_writeback(out, targets, drift, pinned=moved)`: the rows it
   wrote are marked in `moved` AFTER the walk (the walk reads the set as it
   was).  Returns -1 where a Kruskal key is NaN (declined). */
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
    if (isnan(k0[i])) return -1;
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
