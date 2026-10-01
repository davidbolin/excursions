#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

using namespace std;

/*
 Helpers for excursions.regions(): bivariate normal probabilities for the
 pairwise failure probabilities, connected components of a graph, and the
 growth order of a connected region.

 Graphs are symmetric sparse patterns in CSC format (0-based, both triangles,
 no diagonal).
*/

static double Phi(double x) {
  return Rf_pnorm5(x, 0.0, 1.0, 1, 0);
}

/*
 P(X > dh, Y > dk) for a standard bivariate normal vector with correlation r,
 by the method of Genz (2004), Numerical computation of rectangular bivariate
 and trivariate normal and t probabilities, Statistics and Computing 14,
 which is also used in the BVND routine of mvtnorm.
*/
static double bvnd(double dh, double dk, double r) {
  static const double w[3][10] = {
    {0.1713244923791705, 0.3607615730481384, 0.4679139345726904},
    {0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
     0.2031674267230659, 0.2334925365383547, 0.2491470458134029},
    {0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
     0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
     0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
     0.1527533871307259}};
  static const double x[3][10] = {
    {-0.9324695142031522, -0.6612093864662647, -0.2386191860831970},
    {-0.9815606342467191, -0.9041172563704750, -0.7699026741943050,
     -0.5873179542866171, -0.3678314989981802, -0.1252334085114692},
    {-0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
     -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
     -0.5108670019508271, -0.3737060887154196, -0.2277858511416451,
     -0.07652652113349733}};
  const double twopi = 2.0 * M_PI;

  int ng, lg;
  if (fabs(r) < 0.3) {
    ng = 0;
    lg = 3;
  } else if (fabs(r) < 0.75) {
    ng = 1;
    lg = 6;
  } else {
    ng = 2;
    lg = 10;
  }

  double h = dh, k = dk, hk = h * k, bvn = 0.0;
  if (fabs(r) < 0.925) {
    const double hs = (h * h + k * k) / 2.0, asr = asin(r);
    for (int i = 0; i < lg; i++) {
      double sn = sin(asr * (x[ng][i] + 1.0) / 2.0);
      bvn += w[ng][i] * exp((sn * hk - hs) / (1.0 - sn * sn));
      sn = sin(asr * (-x[ng][i] + 1.0) / 2.0);
      bvn += w[ng][i] * exp((sn * hk - hs) / (1.0 - sn * sn));
    }
    bvn = bvn * asr / (2.0 * twopi) + Phi(-h) * Phi(-k);
  } else {
    if (r < 0) {
      k = -k;
      hk = -hk;
    }
    if (fabs(r) < 1) {
      const double as = (1.0 - r) * (1.0 + r);
      double a = sqrt(as);
      const double bs = (h - k) * (h - k);
      const double c = (4.0 - hk) / 8.0, d = (12.0 - hk) / 16.0;
      bvn = a * exp(-(bs / as + hk) / 2.0) *
            (1.0 - c * (bs - as) * (1.0 - d * bs / 5.0) / 3.0 + c * d * as * as / 5.0);
      if (hk > -160) {
        const double b = sqrt(bs);
        bvn -= exp(-hk / 2.0) * sqrt(twopi) * Phi(-b / a) * b *
               (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
      }
      a = a / 2.0;
      for (int i = 0; i < lg; i++) {
        double xs = a * (x[ng][i] + 1.0);
        xs = xs * xs;
        double rs = sqrt(1.0 - xs);
        bvn += a * w[ng][i] *
               (exp(-bs / (2.0 * xs) - hk / (1.0 + rs)) / rs -
                exp(-(bs / xs + hk) / 2.0) * (1.0 + c * xs * (1.0 + d * xs)));
        xs = as * (-x[ng][i] + 1.0) * (-x[ng][i] + 1.0) / 4.0;
        rs = sqrt(1.0 - xs);
        bvn += a * w[ng][i] * exp(-(bs / xs + hk) / 2.0) *
               (exp(-hk * (1.0 - rs) / (2.0 * (1.0 + rs))) / rs -
                (1.0 + c * xs * (1.0 + d * xs)));
      }
      bvn = -bvn / twopi;
    }
    if (r > 0) {
      bvn += Phi(-max(h, k));
    } else {
      bvn = -bvn + max(0.0, Phi(-h) - Phi(-k));
    }
  }
  return min(max(bvn, 0.0), 1.0);
}

// P(X < h, Y < k) for standard bivariate normal vectors with correlations r.
extern "C" SEXP regions_bvn_lower(SEXP Rh, SEXP Rk, SEXP Rr) {
  const int m = Rf_length(Rh);
  if (Rf_length(Rk) != m || Rf_length(Rr) != m)
    Rf_error("regions_bvn_lower: arguments of different lengths.");
  SEXP out = PROTECT(Rf_allocVector(REALSXP, m));
  const double *h = REAL(Rh), *k = REAL(Rk), *r = REAL(Rr);
  double *p = REAL(out);
  for (int i = 0; i < m; i++) {
    p[i] = bvnd(-h[i], -k[i], max(-1.0, min(1.0, r[i])));
  }
  UNPROTECT(1);
  return out;
}

// Connected components of the subgraph induced by the nodes with mask != 0.
// Returns the component of each node, 1, 2, ..., and 0 outside the mask.
extern "C" SEXP regions_components(SEXP Rp, SEXP Ri, SEXP Rmask) {
  const int n = Rf_length(Rp) - 1;
  if (Rf_length(Rmask) != n)
    Rf_error("regions_components: mask has the wrong length.");
  const int *Gp = INTEGER(Rp), *Gi = INTEGER(Ri), *mask = INTEGER(Rmask);
  SEXP out = PROTECT(Rf_allocVector(INTSXP, n));
  int *lab = INTEGER(out);
  fill(lab, lab + n, 0);
  vector<int> stack;
  int nc = 0;
  for (int s = 0; s < n; s++) {
    if (!mask[s] || lab[s])
      continue;
    lab[s] = ++nc;
    stack.push_back(s);
    while (!stack.empty()) {
      const int v = stack.back();
      stack.pop_back();
      for (int q = Gp[v]; q < Gp[v + 1]; q++) {
        const int w = Gi[q];
        if (mask[w] && !lab[w]) {
          lab[w] = nc;
          stack.push_back(w);
        }
      }
    }
  }
  UNPROTECT(1);
  return out;
}

/*
 Order in which a connected region grows from the seed node within the nodes
 with mask != 0. Each step adds the neighbour of the region with the smallest
 score, and ties are broken by the largest z.

 With bound = TRUE, the score of node j is fail[j] - max q[i, j] over its
 neighbours i in the region, where fail[j] is the marginal failure
 probability and q[i, j] the pairwise failure probability of the edge. This
 is the decrease of the Hunter lower bound on the joint probability of the
 region, and the edges used form the maximum spanning tree of the bound
 (Prim's algorithm). With bound = FALSE, the score is -z[j], the marginal
 excursion probability on the standardised scale.

 Returns 1-based node indices.
*/
extern "C" SEXP regions_grow(SEXP Rp, SEXP Ri, SEXP Rq, SEXP Rfail, SEXP Rz,
                             SEXP Rmask, SEXP Rseed, SEXP Rbound) {
  const int n = Rf_length(Rp) - 1;
  if (Rf_length(Rmask) != n || Rf_length(Rfail) != n || Rf_length(Rz) != n)
    Rf_error("regions_grow: arguments have the wrong length.");
  const int *Gp = INTEGER(Rp), *Gi = INTEGER(Ri), *mask = INTEGER(Rmask);
  const double *q = REAL(Rq), *fail = REAL(Rfail), *z = REAL(Rz);
  const int seed = Rf_asInteger(Rseed) - 1;
  const bool bound = Rf_asLogical(Rbound) == TRUE;
  if (seed < 0 || seed >= n || !mask[seed])
    Rf_error("regions_grow: the seed is not in the mask.");

  // Entries (score, -z, node), smallest first. Scores only decrease, so an
  // entry is stale if its score differs from the current score of the node.
  typedef pair<pair<double, double>, int> entry;
  priority_queue<entry, vector<entry>, greater<entry> > heap;
  vector<double> best(n, 0.0), score(n, R_PosInf);
  vector<char> in(n, 0);
  vector<int> order;
  order.reserve(n);

  score[seed] = 0.0;
  heap.push(entry(make_pair(0.0, -z[seed]), seed));
  while (!heap.empty()) {
    const entry e = heap.top();
    heap.pop();
    const int v = e.second;
    if (in[v] || e.first.first != score[v])
      continue;
    in[v] = 1;
    order.push_back(v + 1);
    for (int t = Gp[v]; t < Gp[v + 1]; t++) {
      const int w = Gi[t];
      if (!mask[w] || in[w])
        continue;
      double s;
      if (bound) {
        best[w] = max(best[w], q[t]);
        s = fail[w] - best[w];
      } else {
        s = -z[w];
      }
      if (s < score[w]) {
        score[w] = s;
        heap.push(entry(make_pair(s, -z[w]), w));
      }
    }
  }

  SEXP out = PROTECT(Rf_allocVector(INTSXP, order.size()));
  copy(order.begin(), order.end(), INTEGER(out));
  UNPROTECT(1);
  return out;
}

static int uf_find(vector<int> &parent, int v) {
  while (parent[v] != v) {
    parent[v] = parent[parent[v]];
    v = parent[v];
  }
  return v;
}

/*
 Topographic prominence of the local maxima of z within a set of nodes. The
 nodes are given in decreasing order of z (ties by increasing index), and
 are added one at a time. A node without added neighbours starts a new
 peak. When a node joins the areas of several peaks, all but the highest
 peak end there, and their prominence is the height of the peak minus z at
 the node. The highest peak of each connected component has infinite
 prominence, and nodes that are not local maxima have prominence zero.
*/
extern "C" SEXP regions_prominence(SEXP Rp, SEXP Ri, SEXP Rz, SEXP Rorder) {
  const int n = Rf_length(Rp) - 1;
  const int m = Rf_length(Rorder);
  if (Rf_length(Rz) != n)
    Rf_error("regions_prominence: z has the wrong length.");
  const int *Gp = INTEGER(Rp), *Gi = INTEGER(Ri), *ord = INTEGER(Rorder);
  const double *z = REAL(Rz);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double *prom = REAL(out);
  fill(prom, prom + n, 0.0);

  vector<int> parent(n, -1), peak(n, -1);
  for (int t = 0; t < m; t++) {
    const int v = ord[t] - 1;
    parent[v] = v;
    peak[v] = v;
    for (int s = Gp[v]; s < Gp[v + 1]; s++) {
      const int w = Gi[s];
      if (parent[w] < 0)
        continue;
      int rv = uf_find(parent, v), rw = uf_find(parent, w);
      if (rv == rw)
        continue;
      // The root of v is v itself until v joins a peak, and v is then
      // merged into the area of that peak without ending it
      if (peak[rv] == v) {
        parent[rv] = rw;
        continue;
      }
      // Two peaks meet at v, and the lower one ends
      const int pv = peak[rv], pw = peak[rw];
      const bool v_lower = z[pv] < z[pw] || (z[pv] == z[pw] && pv > pw);
      const int lo = v_lower ? pv : pw;
      prom[lo] = z[lo] - z[v];
      if (v_lower) {
        parent[rv] = rw;
      } else {
        parent[rw] = rv;
      }
    }
    // A node that started a peak and joined nothing is a local maximum, and
    // keeps an infinite prominence unless its peak ends later
    if (peak[uf_find(parent, v)] == v) {
      prom[v] = R_PosInf;
    }
  }
  UNPROTECT(1);
  return out;
}
