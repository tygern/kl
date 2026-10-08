// Exact-integer Kazhdan-Lusztig computation on the Bruhat lower ideal of a
// given element w, in two independent group models of W(D_n):
//   model 1: signed permutations (Gern, thesis Example 1.1.6 / Prop 2.2.3)
//   model 2: geometric representation on the root lattice (Cartan matrix)
// The KL engine stores only extremal pairs (x,y) (R(y)<=R(x), L(y)<=L(x))
// and uses the standard recurrence with right descents (KL79 2.2.c, inverted).
// At the top element it recomputes the full vector P_{.,w} by the direct
// recurrence (no extremal shortcut) as an internal consistency check.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <functional>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
using namespace std;
typedef long long ll;
static const int MAXDEG = 20;
typedef array<ll, MAXDEG> Poly;   // coefficient of q^i at index i

struct Model {
    int n;                                // rank
    int N;                                // ideal size
    vector<array<int32_t, 8>> rmul, lmul; // -1 if outside ideal
    vector<uint8_t> len, rdes, ldes;      // descent bitmasks: bit g = generator g (0-based = Gern s_{g+1})
    vector<vector<int32_t>> coatoms;      // z < y with len(z)=len(y)-1
    int id_e = -1, id_top = -1, id_x = -1;
    vector<string> canon;                 // canonical (min-left-descent) reduced word
};

// ---------------- Model 1: signed permutations ----------------
namespace sp {
typedef array<int8_t, 8> Perm;
int n;
uint32_t key(const Perm& w) { uint32_t k = 0; for (int i = 0; i < n; i++) k |= (uint32_t)((abs(w[i]) - 1) | (w[i] < 0 ? 8 : 0)) << (4 * i); return k; }
int length(const Perm& w) { int l = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (w[i] > w[j]) l++; if (w[i] + w[j] < 0) l++; } return l; }
// right multiplication by Gern's s_{g+1}, g=0..n-1
Perm rmul(Perm w, int g) { if (g == 0) { int8_t a = w[0], b = w[1]; w[0] = -b; w[1] = -a; } else { swap(w[g - 1], w[g]); } return w; }
// left multiplication by Gern's s_{g+1}: acts on values
Perm lmul(Perm w, int g) {
    for (int i = 0; i < n; i++) {
        int v = w[i], a = abs(v), sg = v < 0 ? -1 : 1;
        if (g == 0) { if (a == 1) w[i] = -2 * sg; else if (a == 2) w[i] = -1 * sg; }
        else { if (a == g) w[i] = (g + 1) * sg; else if (a == g + 1) w[i] = g * sg; }
    }
    return w;
}
// right multiplication by reflection: swap positions i<j, optionally negate both
Perm refl(Perm w, int i, int j, bool neg) { swap(w[i], w[j]); if (neg) { w[i] = -w[i]; w[j] = -w[j]; } return w; }
Model build(int rank, const Perm& top, const vector<int>& xword) {
    n = rank; Model M; M.n = n;
    vector<Perm> elems; unordered_map<uint32_t, int32_t> idx;
    auto add = [&](const Perm& w) { auto k = key(w); auto it = idx.find(k); if (it != idx.end()) return it->second; int id = elems.size(); elems.push_back(w); idx[k] = id; return id; };
    add(top);
    for (size_t p = 0; p < elems.size(); p++) {
        Perm w = elems[p]; int lw = length(w);
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) for (int ng = 0; ng < 2; ng++) {
            Perm u = refl(w, i, j, ng); if (length(u) < lw) add(u);
        }
    }
    M.N = elems.size(); M.rmul.resize(M.N); M.lmul.resize(M.N); M.len.resize(M.N); M.rdes.assign(M.N, 0); M.ldes.assign(M.N, 0); M.coatoms.resize(M.N);
    for (int id = 0; id < M.N; id++) {
        const Perm& w = elems[id]; int lw = length(w); M.len[id] = lw;
        for (int g = 0; g < 8; g++) { M.rmul[id][g] = -1; M.lmul[id][g] = -1; }
        for (int g = 0; g < n; g++) {
            Perm u = rmul(w, g); auto it = idx.find(key(u)); if (it != idx.end()) M.rmul[id][g] = it->second;
            if (length(u) < lw) M.rdes[id] |= 1 << g;
            Perm v = lmul(w, g); it = idx.find(key(v)); if (it != idx.end()) M.lmul[id][g] = it->second;
            if (length(v) < lw) M.ldes[id] |= 1 << g;
        }
        // independent descent formulas (Gern Prop 2.2.4) as a check
        uint8_t rd = 0; for (int g = 0; g < n; g++) { if (g == 0 ? (w[0] + w[1] < 0) : (w[g - 1] > w[g])) rd |= 1 << g; }
        if (rd != M.rdes[id]) { fprintf(stderr, "descent formula mismatch\n"); exit(1); }
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) for (int ng = 0; ng < 2; ng++) {
            Perm u = refl(w, i, j, ng); if (length(u) == lw - 1) M.coatoms[id].push_back(idx[key(u)]);
        }
    }
    Perm e; for (int i = 0; i < n; i++) e[i] = i + 1;
    M.id_e = idx[key(e)]; M.id_top = idx[key(top)];
    Perm x = e; for (int g : xword) x = rmul(x, g - 1);
    M.id_x = idx.count(key(x)) ? idx[key(x)] : -1;
    return M;
}
} // namespace sp

// ---------------- Model 2: geometric representation ----------------
namespace geo {
int n;
int A[8][8];                 // Cartan matrix (symmetric, simply laced)
typedef array<int8_t, 64> Mat; // Mat[j*8+i] = coeff of alpha_i in w(alpha_j)
typedef array<int, 8> Vec;
vector<Vec> posroots;
struct H { size_t operator()(const Mat& m) const { return std::hash<string_view>()(string_view((const char*)m.data(), 64)); } };
Vec wapply(const Mat& w, const Vec& v) { Vec r{}; for (int j = 0; j < n; j++) if (v[j]) for (int i = 0; i < n; i++) r[i] += v[j] * w[j * 8 + i]; return r; }
int pairing(const Vec& v, const Vec& b) { int s = 0; for (int i = 0; i < n; i++) if (v[i]) for (int j = 0; j < n; j++) s += v[i] * A[i][j] * b[j]; return s; }
Vec reflectv(const Vec& v, const Vec& b) { int c = pairing(v, b); Vec r = v; for (int i = 0; i < n; i++) r[i] -= c * b[i]; return r; }
bool isneg(const Vec& v) { bool nz = false; for (int i = 0; i < n; i++) { if (v[i] > 0) return false; if (v[i] < 0) nz = true; } return nz; }
int length(const Mat& w) { int l = 0; for (auto& b : posroots) if (isneg(wapply(w, b))) l++; return l; }
Mat identity() { Mat m{}; for (int i = 0; i < n; i++) m[i * 8 + i] = 1; return m; }
// right multiplication by reflection s_b: (w s_b)(alpha_j) = w(alpha_j) - <alpha_j, b> w(b)
Mat rrefl(const Mat& w, const Vec& b) {
    Vec wb = wapply(w, b); Mat r = w;
    for (int j = 0; j < n; j++) { Vec ej{}; ej[j] = 1; int c = pairing(ej, b); if (c) for (int i = 0; i < n; i++) r[j * 8 + i] -= c * wb[i]; }
    return r;
}
Mat rmul(const Mat& w, int g) { Vec b{}; b[g] = 1; return rrefl(w, b); }
// left multiplication by s_g: apply s_g to each column
Mat lmul(const Mat& w, int g) { Vec b{}; b[g] = 1; Mat r = w; for (int j = 0; j < n; j++) { Vec col; for (int i = 0; i < n; i++) col[i] = w[j * 8 + i]; Vec c2 = reflectv(col, b); for (int i = 0; i < n; i++) r[j * 8 + i] = c2[i]; } return r; }
Model build(int rank, const vector<int>& topword, const vector<int>& xword) {
    n = rank;
    // Gern's D_n labeling: s1-s3, s2-s3, s3-s4-...-sn
    for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) A[i][j] = (i == j) ? 2 : 0;
    auto edge = [&](int a, int b) { A[a - 1][b - 1] = A[b - 1][a - 1] = -1; };
    edge(1, 3); edge(2, 3); for (int i = 3; i < n; i++) edge(i, i + 1);
    // positive roots by closure
    { vector<Vec> roots; unordered_map<string, int> seen;
      auto keyv = [&](const Vec& v) { return string((const char*)v.data(), sizeof(Vec)); };
      for (int i = 0; i < n; i++) { Vec v{}; v[i] = 1; roots.push_back(v); seen[keyv(v)] = 1; }
      for (size_t p = 0; p < roots.size(); p++) for (int g = 0; g < n; g++) { Vec b{}; b[g] = 1; Vec r = reflectv(roots[p], b); bool pos = true, nz = false; for (int i = 0; i < n; i++) { if (r[i] < 0) pos = false; if (r[i]) nz = true; } if (pos && nz && !seen.count(keyv(r))) { seen[keyv(r)] = 1; roots.push_back(r); } }
      posroots = roots; fprintf(stderr, "[geo] %zu positive roots\n", posroots.size()); }
    Mat top = identity(); for (int g : topword) top = rmul(top, g - 1);
    Model M; M.n = n;
    vector<Mat> elems; unordered_map<Mat, int32_t, H> idx;
    auto add = [&](const Mat& w) { auto it = idx.find(w); if (it != idx.end()) return it->second; int id = elems.size(); elems.push_back(w); idx[w] = id; return id; };
    add(top);
    for (size_t p = 0; p < elems.size(); p++) { Mat w = elems[p]; int lw = length(w); for (auto& b : posroots) { Mat u = rrefl(w, b); if (length(u) < lw) add(u); } }
    M.N = elems.size(); M.rmul.resize(M.N); M.lmul.resize(M.N); M.len.resize(M.N); M.rdes.assign(M.N, 0); M.ldes.assign(M.N, 0); M.coatoms.resize(M.N);
    for (int id = 0; id < M.N; id++) {
        const Mat& w = elems[id]; int lw = length(w); M.len[id] = lw;
        for (int g = 0; g < 8; g++) { M.rmul[id][g] = -1; M.lmul[id][g] = -1; }
        for (int g = 0; g < n; g++) {
            Mat u = rmul(w, g); auto it = idx.find(u); if (it != idx.end()) M.rmul[id][g] = it->second;
            if (length(u) < lw) M.rdes[id] |= 1 << g;
            Mat v = lmul(w, g); it = idx.find(v); if (it != idx.end()) M.lmul[id][g] = it->second;
            if (length(v) < lw) M.ldes[id] |= 1 << g;
        }
        // independent descent check: s_g right descent iff w(alpha_g) < 0
        uint8_t rd = 0; for (int g = 0; g < n; g++) { Vec col; for (int i = 0; i < n; i++) col[i] = w[g * 8 + i]; if (isneg(col)) rd |= 1 << g; }
        if (rd != M.rdes[id]) { fprintf(stderr, "geo descent mismatch\n"); exit(1); }
        for (auto& b : posroots) { Mat u = rrefl(w, b); if (length(u) == lw - 1) M.coatoms[id].push_back(idx[u]); }
    }
    M.id_e = idx[identity()]; M.id_top = idx[top];
    Mat x = identity(); for (int g : xword) x = rmul(x, g - 1);
    M.id_x = idx.count(x) ? idx[x] : -1;
    return M;
}
} // namespace geo

// canonical words
void canonical_words(Model& M) {
    vector<int> order(M.N); iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return M.len[a] < M.len[b]; });
    M.canon.assign(M.N, "");
    for (int id : order) {
        if (M.len[id] == 0) { M.canon[id] = ""; continue; }
        int g = __builtin_ctz(M.ldes[id]); int u = M.lmul[id][g];
        if (u < 0 || M.len[u] != M.len[id] - 1) { fprintf(stderr, "canon error\n"); exit(1); }
        M.canon[id] = string(1, char('1' + g)) + M.canon[u];
    }
}

// ---------------- KL engine ----------------
struct Entry { int32_t x; int32_t p; };
struct KL {
    const Model& M; int policy;
    vector<vector<Entry>> table;                // per y: sorted extremal entries
    vector<vector<pair<int32_t, ll>>> mulist;   // per y: (z, mu(z,y)) nonzero
    vector<Poly> polys; unordered_map<string, int32_t> polyidx; mutex pmx;
    vector<uint32_t> stampA; uint32_t curA = 0;
    KL(const Model& m, int pol) : M(m), policy(pol) { table.resize(M.N); mulist.resize(M.N); stampA.assign(M.N, 0); }
    int32_t intern(const Poly& p) {
        string k((const char*)p.data(), sizeof(Poly)); lock_guard<mutex> lk(pmx);
        auto it = polyidx.find(k); if (it != polyidx.end()) return it->second;
        int32_t id = polys.size(); polys.push_back(p); polyidx[k] = id; return id;
    }
    static int deg(const Poly& p) { for (int i = MAXDEG - 1; i >= 0; i--) if (p[i]) return i; return -1; }
    bool extremal(int x, int y) const { return (M.rdes[y] & ~M.rdes[x]) == 0 && (M.ldes[y] & ~M.ldes[x]) == 0; }
    // lift u by descents of z until extremal for z; returns -1 if leaves ideal
    int extremalize(int u, int z) const {
        while (true) {
            uint8_t r = M.rdes[z] & ~M.rdes[u];
            if (r) { u = M.rmul[u][__builtin_ctz(r)]; if (u < 0) return -1; continue; }
            uint8_t l = M.ldes[z] & ~M.ldes[u];
            if (l) { u = M.lmul[u][__builtin_ctz(l)]; if (u < 0) return -1; continue; }
            return u;
        }
    }
    const Poly* lookup_ext(int xe, int z) const { // xe already extremal for z
        const auto& t = table[z];
        auto it = lower_bound(t.begin(), t.end(), xe, [](const Entry& e, int v) { return e.x < v; });
        if (it == t.end() || it->x != xe) return nullptr; return &polys[it->p];
    }
    const Poly* lookup(int u, int z) const { if (u < 0) return nullptr; int xe = extremalize(u, z); if (xe < 0) return nullptr; return lookup_ext(xe, z); }
    // direct recurrence for P_{x,y} with descent s, v = ys (no shortcut on x)
    Poly compute(int x, int y, int s, int v) const {
        Poly P{}; int c = (M.rdes[x] >> s) & 1; int xs = M.rmul[x][s];
        if (const Poly* a = lookup(xs, v)) for (int i = 0; i + (1 - c) < MAXDEG; i++) P[i + 1 - c] += (*a)[i];
        if (const Poly* b = lookup(x, v)) for (int i = 0; i + c < MAXDEG; i++) P[i + c] += (*b)[i];
        for (auto& zm : mulist[v]) {
            int z = zm.first; if (!((M.rdes[z] >> s) & 1)) continue;
            const Poly* pz = lookup(x, z); if (!pz) continue;
            int sh = (M.len[y] - M.len[z]) / 2; ll mu = zm.second;
            for (int i = 0; i + sh < MAXDEG; i++) P[i + sh] -= mu * (*pz)[i];
        }
        return P;
    }
    int choose(int y) const { uint8_t r = M.rdes[y]; if (policy == 0) return __builtin_ctz(r); return 31 - __builtin_clz((unsigned)r); }
    // enumerate [e,v] from extremal entries by descending closure; returns list
    void lower_ideal(int v, vector<int32_t>& out) {
        curA++; out.clear();
        for (auto& e : table[v]) { if (stampA[e.x] != curA) { stampA[e.x] = curA; out.push_back(e.x); } }
        for (size_t p = 0; p < out.size(); p++) {
            int x = out[p];
            uint8_t r = M.rdes[v] & M.rdes[x]; while (r) { int g = __builtin_ctz(r); r &= r - 1; int u = M.rmul[x][g]; if (u >= 0 && stampA[u] != curA) { stampA[u] = curA; out.push_back(u); } }
            uint8_t l = M.ldes[v] & M.ldes[x]; while (l) { int g = __builtin_ctz(l); l &= l - 1; int u = M.lmul[x][g]; if (u >= 0 && stampA[u] != curA) { stampA[u] = curA; out.push_back(u); } }
        }
    }
    void process(int y, vector<int32_t>& buf, vector<pair<int32_t, Poly>>& res) {
        res.clear();
        if (M.len[y] == 0) { Poly one{}; one[0] = 1; res.push_back({y, one}); return; }
        int s = choose(y); int v = M.rmul[y][s];
        lower_ideal(v, buf);
        size_t nv = buf.size();
        // [e,y] = [e,v] u [e,v]s ; stampA currently marks [e,v]
        for (size_t p = 0; p < nv; p++) { int x = buf[p]; int xs = M.rmul[x][s]; if (xs >= 0 && stampA[xs] != curA) { stampA[xs] = curA; buf.push_back(xs); } }
        for (int x : buf) if (extremal(x, y)) res.push_back({x, compute(x, y, s, v)});
        sort(res.begin(), res.end(), [](auto& a, auto& b) { return a.first < b.first; });
    }
    void finish(int y, const vector<pair<int32_t, Poly>>& res, ll& npairs) {
        auto& t = table[y]; t.reserve(res.size());
        for (auto& r : res) { t.push_back({r.first, intern(r.second)}); }
        npairs += t.size();
        // mu list: coatoms plus extremal with odd codim and full degree
        vector<pair<int32_t, ll>> mu;
        for (int z : M.coatoms[y]) mu.push_back({z, 1});
        for (auto& e : t) { int cod = M.len[y] - M.len[e.x]; if (cod <= 1 || cod % 2 == 0) continue; const Poly& p = polys[e.p]; int d = (cod - 1) / 2; if (p[d]) mu.push_back({e.x, p[d]}); }
        sort(mu.begin(), mu.end()); mu.erase(unique(mu.begin(), mu.end()), mu.end());
        for (size_t i = 1; i < mu.size(); i++) if (mu[i].first == mu[i - 1].first) { fprintf(stderr, "mu dup conflict\n"); exit(1); }
        mulist[y] = mu;
    }
    void run(int nthreads) {
        vector<int> order(M.N); iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b) { return M.len[a] < M.len[b]; });
        int maxlen = M.len[order.back()];
        ll npairs = 0; auto t0 = chrono::steady_clock::now();
        size_t pos = 0;
        for (int L = 0; L <= maxlen; L++) {
            size_t start = pos; while (pos < order.size() && M.len[order[pos]] == L) pos++;
            vector<int> layer(order.begin() + start, order.begin() + pos);
            vector<vector<pair<int32_t, Poly>>> results(layer.size());
            atomic<size_t> next{0};
            auto worker = [&]() {
                // each thread needs its own stamp array -> use a thread-local copy of KL state for stamps
                vector<int32_t> buf; vector<uint32_t> stamps(M.N, 0); uint32_t cur = 0;
                while (true) {
                    size_t i = next.fetch_add(1); if (i >= layer.size()) break;
                    int y = layer[i]; auto& res = results[i]; res.clear();
                    if (M.len[y] == 0) { Poly one{}; one[0] = 1; res.push_back({y, one}); continue; }
                    int s = choose(y); int v = M.rmul[y][s];
                    cur++; buf.clear();
                    for (auto& e : table[v]) if (stamps[e.x] != cur) { stamps[e.x] = cur; buf.push_back(e.x); }
                    for (size_t p = 0; p < buf.size(); p++) {
                        int x = buf[p];
                        uint8_t r = M.rdes[v] & M.rdes[x]; while (r) { int g = __builtin_ctz(r); r &= r - 1; int u = M.rmul[x][g]; if (u >= 0 && stamps[u] != cur) { stamps[u] = cur; buf.push_back(u); } }
                        uint8_t l = M.ldes[v] & M.ldes[x]; while (l) { int g = __builtin_ctz(l); l &= l - 1; int u = M.lmul[x][g]; if (u >= 0 && stamps[u] != cur) { stamps[u] = cur; buf.push_back(u); } }
                    }
                    size_t nv = buf.size();
                    for (size_t p = 0; p < nv; p++) { int x = buf[p]; int xs = M.rmul[x][s]; if (xs >= 0 && stamps[xs] != cur) { stamps[xs] = cur; buf.push_back(xs); } }
                    for (int x : buf) if (extremal(x, y)) res.push_back({x, compute(x, y, s, v)});
                    sort(res.begin(), res.end(), [](auto& a, auto& b) { return a.first < b.first; });
                }
            };
            vector<thread> th; for (int t = 0; t < nthreads; t++) th.emplace_back(worker); for (auto& t : th) t.join();
            for (size_t i = 0; i < layer.size(); i++) finish(layer[i], results[i], npairs);
            double el = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
            fprintf(stderr, "  layer %2d: %7zu elements, cumulative extremal pairs %lld, distinct polys %zu, %.1fs\n", L, layer.size(), npairs, polys.size(), el);
        }
    }
};

static string polystr(const Poly& p) { string s; int d = KL::deg(p); if (d < 0) return "0"; for (int i = 0; i <= d; i++) { if (!p[i]) continue; if (!s.empty()) s += (p[i] < 0 ? " - " : " + "); ll a = llabs(p[i]); if (i == 0) s += to_string(a); else { if (a != 1) s += to_string(a); s += "q"; if (i > 1) s += "^" + to_string(i); } } return s; }

int main(int argc, char** argv) {
    // usage: kl <model:sp|geo> <n> <policy 0|1> <threads> <outfile>
    if (argc < 6) { fprintf(stderr, "usage\n"); return 1; }
    string model = argv[1]; int n = atoi(argv[2]); int policy = atoi(argv[3]); int nth = atoi(argv[4]); string outfile = argv[5];
    // Gern's w_n (n even) and x_n
    vector<int> topword, xword; // Gern labels 1..n
    auto br = [&](int j, int i) { vector<int> w; if (i == 0) for (int k = j; k >= 1; k--) w.push_back(k); else for (int k = j; k >= i; k--) w.push_back(k); return w; };
    { int k = n / 2 - 2; for (int j = 2; j <= n; j += 2) { auto b = br(j, 0); topword.insert(topword.end(), b.begin(), b.end()); }
      for (int i = 0; i <= k; i++) { auto b = br(n - k + i, n - 2 * k + 2 * i); topword.insert(topword.end(), b.begin(), b.end()); } }
    xword = {1, 2}; for (int j = 4; j <= n; j += 2) xword.push_back(j);
    bool override = false;
    if (argc > 6 && strlen(argv[6]) > 0) { override = true; topword.clear(); for (char* c = argv[6]; *c; c++) topword.push_back(*c - '0'); xword = {1}; }
    fprintf(stderr, "top word (Gern Lemma 2.3.4):"); for (int g : topword) fprintf(stderr, " %d", g); fprintf(stderr, "  (length %zu)\n", topword.size());
    Model M;
    if (model == "sp") {
        sp::Perm top{};
        if (override) { sp::n = n; for (int i = 0; i < n; i++) top[i] = i + 1; for (int g : topword) top = sp::rmul(top, g - 1); }
        else // Corollary 2.2.19, n even: w(1)=(-1)^{n/2}, odd i>1 -> i, even i -> -(n+2-i)
        for (int i = 1; i <= n; i++) top[i - 1] = (i == 1) ? ((n / 2) % 2 ? -1 : 1) : (i % 2 ? i : -(n + 2 - i));
        fprintf(stderr, "[sp] top signed permutation:"); for (int i = 0; i < n; i++) fprintf(stderr, " %d", top[i]); fprintf(stderr, "\n");
        M = sp::build(n, top, xword);
        // check Lemma 2.3.4 word gives the same element
        sp::Perm e; for (int i = 0; i < n; i++) e[i] = i + 1; sp::Perm w = e; for (int g : topword) w = sp::rmul(w, g - 1);
        fprintf(stderr, "[sp] Lemma 2.3.4 word evaluates to:"); for (int i = 0; i < n; i++) fprintf(stderr, " %d", w[i]); fprintf(stderr, " length %d; equals Corollary 2.2.19 form: %s\n", sp::length(w), w == top ? "yes" : "NO");
    } else {
        M = geo::build(n, topword, xword);
    }
    fprintf(stderr, "[%s] ideal size %d, len(top)=%d, len(x)=%d, id_x=%d\n", model.c_str(), M.N, M.len[M.id_top], M.id_x >= 0 ? M.len[M.id_x] : -1, M.id_x);
    canonical_words(M);
    fprintf(stderr, "[%s] canon(top)=%s canon(x)=%s\n", model.c_str(), M.canon[M.id_top].c_str(), M.id_x >= 0 ? M.canon[M.id_x].c_str() : "?");
    // rank generating function of ideal
    { vector<int> cnt(64, 0); for (int i = 0; i < M.N; i++) cnt[M.len[i]]++; fprintf(stderr, "[%s] ideal rank vector:", model.c_str()); for (int i = 0; i <= M.len[M.id_top]; i++) fprintf(stderr, " %d", cnt[i]); fprintf(stderr, "\n"); }
    KL kl(M, policy);
    kl.run(nth);
    int top = M.id_top;
    // final: full vector P_{x,top} for all x in ideal via direct recurrence for every right descent of top
    FILE* f = fopen(outfile.c_str(), "w");
    vector<int> order(M.N); iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return M.canon[a] < M.canon[b]; });
    uint8_t rd = M.rdes[top]; int ndesc = __builtin_popcount(rd);
    int mismatches = 0; vector<ll> intervalrank(64, 0);
    for (int x : order) {
        Poly ref{}; bool first = true; const Poly* ext = kl.lookup(x, top);
        for (int s = 0; s < n; s++) if ((rd >> s) & 1) {
            Poly P = kl.compute(x, top, s, M.rmul[top][s]);
            if (first) { ref = P; first = false; } else if (P != ref) mismatches++;
        }
        Poly extP{}; if (ext) extP = *ext;
        if (extP != ref) mismatches++;
        if (KL::deg(ref) >= 0) intervalrank[M.len[x]]++;
        fprintf(f, "%s\t%d\t%s\n", M.canon[x].c_str(), (int)M.len[x], polystr(ref).c_str());
    }
    fclose(f);
    fprintf(stderr, "[%s] top vector: %d descents checked per x, mismatches=%d\n", model.c_str(), ndesc, mismatches);
    if (argc > 7) { FILE* d = fopen(argv[7], "w"); for (int y = 0; y < M.N; y++) for (auto& e : kl.table[y]) fprintf(d, "%s\t%s\t%s\n", M.canon[e.x].c_str(), M.canon[y].c_str(), polystr(kl.polys[e.p]).c_str()); fclose(d); }
    // report
    const Poly* Px = kl.lookup(M.id_x, top);
    Poly Pe = kl.compute(M.id_e, top, __builtin_ctz(rd), M.rmul[top][__builtin_ctz(rd)]);
    printf("model=%s n=%d policy=%d\n", model.c_str(), n, policy);
    printf("ideal size = %d\n", M.N);
    printf("len(w) = %d, len(x) = %d\n", M.len[top], M.len[M.id_x]);
    printf("P(x,w) = %s\n", Px ? polystr(*Px).c_str() : "?");
    printf("P(e,w) [direct recurrence] = %s\n", polystr(Pe).c_str());
    int cod = M.len[top] - M.len[M.id_x];
    if (cod % 2 == 0) printf("mu(x,w) = 0 (codimension %d even)\n", cod);
    else { int md = (cod - 1) / 2; printf("mu(x,w) = coefficient of q^%d in P(x,w) = %lld\n", md, Px ? (*Px)[md] : -1); }
    { vector<ll> ir(64, 0); ll tot = 0; for (int z = 0; z < M.N; z++) if (kl.lookup(M.id_x, z)) { ir[M.len[z]]++; tot++; }
      printf("true interval [x,w] (z with P(x,z)!=0) rank vector:"); for (int i = M.len[M.id_x]; i <= M.len[top]; i++) printf(" %lld", ir[i]); printf("  (total %lld)\n", tot); }
    printf("deg P(e,w) = %d\n", KL::deg(Pe));
    printf("elements of [e,w] with length >= len(x), by length:"); { int lx = M.len[M.id_x]; ll tot = 0; for (int i = lx; i <= M.len[top]; i++) { printf(" %lld", intervalrank[i]); tot += intervalrank[i]; } printf("  (total %lld)\n", tot); }
    printf("distinct polynomials in table: %zu\n", kl.polys.size());
    // max coefficient / degree sanity
    ll maxc = 0; int maxd = 0; for (auto& p : kl.polys) { for (int i = 0; i < MAXDEG; i++) maxc = max(maxc, llabs(p[i])); maxd = max(maxd, KL::deg(p)); }
    printf("max |coefficient| = %lld, max degree = %d\n", maxc, maxd);
    // positivity and constant-term checks over all stored pairs
    ll bad = 0; for (auto& p : kl.polys) { if (p[0] != 1) bad++; for (int i = 0; i < MAXDEG; i++) if (p[i] < 0) bad++; }
    printf("polys with nonpositive coefficient or constant term != 1: %lld\n", bad);
    return 0;
}
