// Independent brute-force enumeration of all two-sided terminal elements in
// finite E_n (n = 6, 7), paper labelling: chain 0..n-2, node n-1 attached to 2.
// Elements stored as integer matrices (columns = images of simple roots) in
// the simple-root basis.  BFS over right multiplication enumerates the group.
// Terminal test: for s in R(w) and t ~ s, w(alpha_s + alpha_t) > 0; left test
// uses w^{-1} = C^{-1} w^T C (computed with 2*C^{-1}, which is integral for E7;
// for E6 we use 3*C^{-1}).
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <array>
#include <string>
#include <algorithm>
using namespace std;

static int N;
static int C[8][8];
static vector<int> adj[8];
static long long DCinv[8][8]; // D * C^{-1}, integral
static long long D;

struct Mat { signed char a[8][8]; }; // a[j][i] = coordinate i of w(alpha_j)

struct MatHash { size_t operator()(const Mat& m) const {
    size_t h = 1469598103934665603ULL;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) { h ^= (unsigned char)m.a[j][i]; h *= 1099511628211ULL; }
    return h; } };
struct MatEq { bool operator()(const Mat& x, const Mat& y) const { return memcmp(x.a, y.a, sizeof(x.a)) == 0; } };

static Mat identity() { Mat m; memset(m.a, 0, sizeof(m.a)); for (int i = 0; i < N; i++) m.a[i][i] = 1; return m; }
static Mat rmul(const Mat& m, int s) {
    Mat r = m;
    for (int i = 0; i < N; i++) r.a[s][i] = -m.a[s][i];
    for (int t : adj[s]) for (int i = 0; i < N; i++) r.a[t][i] = m.a[t][i] + m.a[s][i];
    return r;
}
static bool isneg(const signed char* v) { bool nz = false; for (int i = 0; i < N; i++) { if (v[i] > 0) return false; if (v[i] < 0) nz = true; } return nz; }
static bool ispos(const long long* v) { bool nz = false; for (int i = 0; i < N; i++) { if (v[i] < 0) return false; if (v[i] > 0) nz = true; } return nz; }

static bool right_terminal(const Mat& m, int* desc, int& nd) {
    nd = 0;
    for (int s = 0; s < N; s++) if (isneg(m.a[s])) desc[nd++] = s;
    for (int k = 0; k < nd; k++) { int s = desc[k];
        for (int t : adj[s]) { long long v[8]; for (int i = 0; i < N; i++) v[i] = m.a[s][i] + m.a[t][i]; if (!ispos(v)) return false; } }
    return true;
}
// inverse: columns of w^{-1}: w^{-1}(alpha_j) = C^{-1} w^T C alpha_j.  Compute D*inverse then divide.
static Mat inverse(const Mat& m) {
    // W = matrix with W[i][j] = m.a[j][i].  inv = Cinv * W^T * C.
    long long WT_C[8][8]; // (W^T C)[i][j] = sum_k W[k][i] C[k][j] = sum_k m.a[i][k] C[k][j]
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) { long long s = 0; for (int k = 0; k < N; k++) s += (long long)m.a[i][k] * C[k][j]; WT_C[i][j] = s; }
    Mat r;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) { long long s = 0; for (int k = 0; k < N; k++) s += DCinv[i][k] * WT_C[k][j];
        if (s % D != 0) { fprintf(stderr, "non-integral inverse\n"); exit(1); } r.a[j][i] = (signed char)(s / D); }
    return r;
}

int main(int argc, char** argv) {
    N = atoi(argv[1]);
    FILE* dump = (argc > 2) ? fopen(argv[2], "w") : nullptr;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) C[i][j] = (i == j) ? 2 : 0;
    auto edge = [&](int a, int b) { C[a][b] = C[b][a] = -1; adj[a].push_back(b); adj[b].push_back(a); };
    for (int i = 0; i + 1 <= N - 2; i++) edge(i, i + 1);
    edge(2, N - 1);
    // compute D*C^{-1} by Gaussian elimination with fractions (small): use long double then round, verify.
    D = (N == 6) ? 3 : (N == 7 ? 2 : 1);
    { long double A[8][16];
      for (int i = 0; i < N; i++) { for (int j = 0; j < N; j++) A[i][j] = C[i][j]; for (int j = 0; j < N; j++) A[i][N + j] = (i == j); }
      for (int c = 0; c < N; c++) { int p = c; while (A[p][c] == 0) p++; for (int j = 0; j < 2 * N; j++) swap(A[c][j], A[p][j]);
        long double piv = A[c][c]; for (int j = 0; j < 2 * N; j++) A[c][j] /= piv;
        for (int r = 0; r < N; r++) if (r != c && A[r][c] != 0) { long double f = A[r][c]; for (int j = 0; j < 2 * N; j++) A[r][j] -= f * A[c][j]; } }
      for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) { long double v = A[i][N + j] * D; long long iv = llround(v); if (fabsl(v - iv) > 1e-9) { fprintf(stderr, "Cinv not integral with D\n"); return 1; } DCinv[i][j] = iv; }
      // verify DCinv * C = D I
      for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) { long long s = 0; for (int k = 0; k < N; k++) s += DCinv[i][k] * C[k][j]; if (s != (i == j ? D : 0)) { fprintf(stderr, "Cinv check failed\n"); return 1; } }
    }
    unordered_map<Mat, int, MatHash, MatEq> index; // element -> length
    vector<Mat> elems; vector<int> len; vector<int> parent; vector<signed char> lastgen;
    Mat e = identity(); index.emplace(e, 0); elems.push_back(e); len.push_back(0); parent.push_back(-1); lastgen.push_back(-1);
    for (size_t k = 0; k < elems.size(); k++) {
        Mat m = elems[k];
        for (int s = 0; s < N; s++) {
            if (isneg(m.a[s])) continue; // descent: shorter
            Mat r = rmul(m, s);
            if (index.find(r) == index.end()) { index.emplace(r, len[k] + 1); elems.push_back(r); len.push_back(len[k] + 1); parent.push_back((int)k); lastgen.push_back((signed char)s); }
        }
    }
    printf("E%d: group order %zu, max length %d\n", N, elems.size(), *max_element(len.begin(), len.end()));
    long long rt = 0, both = 0, noncomm = 0;
    vector<pair<int,string>> nc;
    for (size_t k = 0; k < elems.size(); k++) {
        int desc[8], nd;
        if (!right_terminal(elems[k], desc, nd)) continue;
        rt++;
        if (dump) { string w; int cur = (int)k; while (cur > 0) { w.push_back('0' + lastgen[cur]); cur = parent[cur]; } reverse(w.begin(), w.end()); fprintf(dump, "%s\n", w.c_str()); }
        Mat inv = inverse(elems[k]); int d2[8], nd2;
        if (!right_terminal(inv, d2, nd2)) continue;
        both++;
        // commuting product?  descents pairwise nonadjacent and length == nd
        bool indep = true; for (int a = 0; a < nd; a++) for (int b : adj[desc[a]]) for (int c = 0; c < nd; c++) if (desc[c] == b) indep = false;
        bool comm = indep && (len[k] == nd);
        if (comm) {
            // check product of the descents equals the element
            Mat p = identity(); for (int a = 0; a < nd; a++) p = rmul(p, desc[a]);
            if (memcmp(p.a, elems[k].a, sizeof(p.a)) != 0) comm = false;
        }
        if (!comm) { noncomm++;
            string w; int cur = (int)k; while (cur > 0) { w.push_back('0' + lastgen[cur]); cur = parent[cur]; } reverse(w.begin(), w.end());
            string L, R; for (int a = 0; a < nd; a++) R.push_back('0' + desc[a]); for (int a = 0; a < nd2; a++) L.push_back('0' + d2[a]);
            nc.push_back({len[k], w + " L=" + L + " R=" + R}); }
    }
    sort(nc.begin(), nc.end());
    printf("right-terminal: %lld, two-sided terminal: %lld, noncommuting terminal: %lld, commuting: %lld\n", rt, both, noncomm, both - noncomm);
    for (auto& p : nc) printf("  len %d  %s\n", p.first, p.second.c_str());
    return 0;
}
