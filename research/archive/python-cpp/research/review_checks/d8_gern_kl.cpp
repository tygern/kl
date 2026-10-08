// Kazhdan-Lusztig polynomials on the Bruhat lower ideal of Gern's bad element
// w_n in W(D_n), n even, in two independent models of the group:
//   model sp:  signed permutations, Gern's conventions (thesis Example 1.1.6,
//              Proposition 2.2.4): s_1 acts on the right by (w_1,w_2) -> (-w_2,-w_1),
//              s_i (i >= 2) swaps positions i-1 and i;
//   model geo: the geometric representation on the root lattice (Cartan matrix
//              of D_n with Gern's labels: nodes 1 and 2 attached to 3, chain 3-4-...-n).
// The engine stores only extremal pairs (x,y) with R(y) subset R(x) and
// L(y) subset L(x) and uses the right-descent recurrence (Kazhdan-Lusztig
// 1979, (2.2.c) transported to the right).  At the top element it recomputes
// P_{x,w} for every x in the ideal by the direct recurrence for every right
// descent of w, and all values must agree with the stored ones.
//
// Certified here: D6 (P_{x_6,w_6} = 1+6q+11q^2+6q^3+q^4+q^5, used by the
// manuscript's main theorem; Gern obtained it with du Cloux's Coxeter) and
// D8 (P_{x_8,w_8} = P_{e,w_8} = 1+12q+59q^2+154q^3+233q^4+221q^5+147q^6+70q^7+20q^8+2q^9,
// mu(x_8,w_8) = 0, used in the archived note kl-transfer).  Gern's element
// w_n is taken from Corollary 2.2.19 (signed permutation) and checked against
// the bracket word of Lemma 2.3.4; x_n = s_1 s_2 s_4 s_6 ... s_n.
//
// Build: c++ -std=c++17 -O3 d8_gern_kl.cpp -o d8_gern_kl
// Run:   ./d8_gern_kl <output.json> [threads]
// The D8 ideal has 265,760 elements; each model takes a few seconds with
// several threads.  Exact 64-bit integers; coefficients stay far below 2^62.

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <memory>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
using namespace std;
typedef long long ll;
static const int MAXDEG = 14;  // degrees stay below (l(w)-1)/2 = 12 for D8
typedef array<ll, MAXDEG> Poly;

struct Model {
    int n = 0, N = 0;
    vector<array<int32_t, 8>> rmul, lmul;  // -1 if outside the ideal
    vector<uint8_t> len, rdes, ldes;       // descent masks: bit g = Gern's s_{g+1}
    vector<vector<int32_t>> coatoms;       // z < y with l(z) = l(y) - 1
    int id_e = -1, id_top = -1, id_x = -1;
    vector<string> canon;                  // canonical reduced word (smallest left descent first)
};

namespace sp {
typedef array<int8_t, 8> Perm;
int n;
uint32_t key(const Perm& w) { uint32_t k = 0; for (int i = 0; i < n; i++) k |= (uint32_t)((abs(w[i]) - 1) | (w[i] < 0 ? 8 : 0)) << (4 * i); return k; }
int length(const Perm& w) { int l = 0; for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) { if (w[i] > w[j]) l++; if (w[i] + w[j] < 0) l++; } return l; }
Perm rmul(Perm w, int g) { if (g == 0) { int8_t a = w[0], b = w[1]; w[0] = -b; w[1] = -a; } else swap(w[g - 1], w[g]); return w; }
Perm lmul(Perm w, int g) {
    for (int i = 0; i < n; i++) {
        int v = w[i], a = abs(v), sg = v < 0 ? -1 : 1;
        if (g == 0) { if (a == 1) w[i] = -2 * sg; else if (a == 2) w[i] = -1 * sg; }
        else { if (a == g) w[i] = (g + 1) * sg; else if (a == g + 1) w[i] = g * sg; }
    }
    return w;
}
Perm refl(Perm w, int i, int j, bool neg) { swap(w[i], w[j]); if (neg) { w[i] = -w[i]; w[j] = -w[j]; } return w; }
Model build(int rank, const Perm& top, const vector<int>& xword) {
    n = rank; Model M; M.n = n;
    vector<Perm> elems; unordered_map<uint32_t, int32_t> idx;
    auto add = [&](const Perm& w) { auto k = key(w); auto it = idx.find(k); if (it != idx.end()) return it->second; int id = elems.size(); elems.push_back(w); idx[k] = id; return id; };
    add(top);
    for (size_t p = 0; p < elems.size(); p++) {
        Perm w = elems[p]; int lw = length(w);
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) for (int ng = 0; ng < 2; ng++) { Perm u = refl(w, i, j, ng); if (length(u) < lw) add(u); }
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
        uint8_t rd = 0; for (int g = 0; g < n; g++) if (g == 0 ? (w[0] + w[1] < 0) : (w[g - 1] > w[g])) rd |= 1 << g;  // Gern Prop. 2.2.4
        if (rd != M.rdes[id]) { fprintf(stderr, "descent formula mismatch\n"); exit(1); }
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) for (int ng = 0; ng < 2; ng++) { Perm u = refl(w, i, j, ng); if (length(u) == lw - 1) M.coatoms[id].push_back(idx[key(u)]); }
    }
    Perm e; for (int i = 0; i < n; i++) e[i] = i + 1;
    M.id_e = idx[key(e)]; M.id_top = idx[key(top)];
    Perm x = e; for (int g : xword) x = rmul(x, g - 1);
    M.id_x = idx.count(key(x)) ? idx[key(x)] : -1;
    return M;
}
}  // namespace sp

namespace geo {
int n; int A[8][8];
typedef array<int8_t, 64> Mat;  // Mat[j*8+i] = coefficient of alpha_i in w(alpha_j)
typedef array<int, 8> Vec;
vector<Vec> posroots;
struct H { size_t operator()(const Mat& m) const { return std::hash<string_view>()(string_view((const char*)m.data(), 64)); } };
Vec wapply(const Mat& w, const Vec& v) { Vec r{}; for (int j = 0; j < n; j++) if (v[j]) for (int i = 0; i < n; i++) r[i] += v[j] * w[j * 8 + i]; return r; }
int pairing(const Vec& v, const Vec& b) { int s = 0; for (int i = 0; i < n; i++) if (v[i]) for (int j = 0; j < n; j++) s += v[i] * A[i][j] * b[j]; return s; }
Vec reflectv(const Vec& v, const Vec& b) { int c = pairing(v, b); Vec r = v; for (int i = 0; i < n; i++) r[i] -= c * b[i]; return r; }
bool isneg(const Vec& v) { bool nz = false; for (int i = 0; i < n; i++) { if (v[i] > 0) return false; if (v[i] < 0) nz = true; } return nz; }
int length(const Mat& w) { int l = 0; for (auto& b : posroots) if (isneg(wapply(w, b))) l++; return l; }
Mat identity() { Mat m{}; for (int i = 0; i < n; i++) m[i * 8 + i] = 1; return m; }
Mat rrefl(const Mat& w, const Vec& b) { Vec wb = wapply(w, b); Mat r = w; for (int j = 0; j < n; j++) { Vec ej{}; ej[j] = 1; int c = pairing(ej, b); if (c) for (int i = 0; i < n; i++) r[j * 8 + i] -= c * wb[i]; } return r; }
Mat rmul(const Mat& w, int g) { Vec b{}; b[g] = 1; return rrefl(w, b); }
Mat lmul(const Mat& w, int g) { Vec b{}; b[g] = 1; Mat r = w; for (int j = 0; j < n; j++) { Vec col; for (int i = 0; i < n; i++) col[i] = w[j * 8 + i]; Vec c2 = reflectv(col, b); for (int i = 0; i < n; i++) r[j * 8 + i] = c2[i]; } return r; }
Model build(int rank, const vector<int>& topword, const vector<int>& xword) {
    n = rank;
    for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) A[i][j] = (i == j) ? 2 : 0;
    auto edge = [&](int a, int b) { A[a - 1][b - 1] = A[b - 1][a - 1] = -1; };
    edge(1, 3); edge(2, 3); for (int i = 3; i < n; i++) edge(i, i + 1);
    { vector<Vec> roots; unordered_map<string, int> seen;
      auto keyv = [&](const Vec& v) { return string((const char*)v.data(), sizeof(Vec)); };
      for (int i = 0; i < n; i++) { Vec v{}; v[i] = 1; roots.push_back(v); seen[keyv(v)] = 1; }
      for (size_t p = 0; p < roots.size(); p++) for (int g = 0; g < n; g++) { Vec b{}; b[g] = 1; Vec r = reflectv(roots[p], b); bool pos = true, nz = false; for (int i = 0; i < n; i++) { if (r[i] < 0) pos = false; if (r[i]) nz = true; } if (pos && nz && !seen.count(keyv(r))) { seen[keyv(r)] = 1; roots.push_back(r); } }
      posroots = roots;
      if ((int)posroots.size() != n * (n - 1)) { fprintf(stderr, "wrong number of positive roots\n"); exit(1); } }
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
        uint8_t rd = 0; for (int g = 0; g < n; g++) { Vec col; for (int i = 0; i < n; i++) col[i] = w[g * 8 + i]; if (isneg(col)) rd |= 1 << g; }
        if (rd != M.rdes[id]) { fprintf(stderr, "geo descent mismatch\n"); exit(1); }
        for (auto& b : posroots) { Mat u = rrefl(w, b); if (length(u) == lw - 1) M.coatoms[id].push_back(idx[u]); }
    }
    M.id_e = idx[identity()]; M.id_top = idx[top];
    Mat x = identity(); for (int g : xword) x = rmul(x, g - 1);
    M.id_x = idx.count(x) ? idx[x] : -1;
    return M;
}
}  // namespace geo

void canonical_words(Model& M) {
    vector<int> order(M.N); iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return M.len[a] < M.len[b]; });
    M.canon.assign(M.N, "");
    for (int id : order) {
        if (M.len[id] == 0) continue;
        int g = __builtin_ctz(M.ldes[id]); int u = M.lmul[id][g];
        if (u < 0 || M.len[u] != M.len[id] - 1) { fprintf(stderr, "canonical word error\n"); exit(1); }
        M.canon[id] = string(1, char('1' + g)) + M.canon[u];
    }
}

struct Entry { int32_t x; int32_t p; };
// Append-only polynomial storage whose element addresses never move.  Worker
// threads look up polynomials of lower layers (interned before the preceding
// join) while new ones are appended under the mutex; a std::vector would
// reallocate and invalidate those addresses.  The chunk table has fixed size.
struct PolyStore {
    static const size_t CH = 1 << 12, MAXCH = 1 << 16;
    vector<unique_ptr<Poly[]>> chunks; size_t n = 0;
    PolyStore() : chunks(MAXCH) {}
    size_t size() const { return n; }
    const Poly& operator[](size_t i) const { return chunks[i / CH][i % CH]; }
    size_t push(const Poly& p) {
        if (n % CH == 0) { if (n / CH >= MAXCH) { fprintf(stderr, "too many polynomials\n"); exit(1); } chunks[n / CH].reset(new Poly[CH]); }
        chunks[n / CH][n % CH] = p; return n++;
    }
};
struct KL {
    const Model& M;
    vector<vector<Entry>> table;
    vector<vector<pair<int32_t, ll>>> mulist;
    PolyStore polys; unordered_map<string, int32_t> polyidx; mutex pmx;
    KL(const Model& m) : M(m) { table.resize(M.N); mulist.resize(M.N); }
    int32_t intern(const Poly& p) { string k((const char*)p.data(), sizeof(Poly)); lock_guard<mutex> lk(pmx); auto it = polyidx.find(k); if (it != polyidx.end()) return it->second; int32_t id = polys.push(p); polyidx[k] = id; return id; }
    static int deg(const Poly& p) { for (int i = MAXDEG - 1; i >= 0; i--) if (p[i]) return i; return -1; }
    bool extremal(int x, int y) const { return (M.rdes[y] & ~M.rdes[x]) == 0 && (M.ldes[y] & ~M.ldes[x]) == 0; }
    int extremalize(int u, int z) const {
        while (true) {
            uint8_t r = M.rdes[z] & ~M.rdes[u]; if (r) { u = M.rmul[u][__builtin_ctz(r)]; if (u < 0) return -1; continue; }
            uint8_t l = M.ldes[z] & ~M.ldes[u]; if (l) { u = M.lmul[u][__builtin_ctz(l)]; if (u < 0) return -1; continue; }
            return u;
        }
    }
    const Poly* lookup_ext(int xe, int z) const { const auto& t = table[z]; auto it = lower_bound(t.begin(), t.end(), xe, [](const Entry& e, int v) { return e.x < v; }); if (it == t.end() || it->x != xe) return nullptr; return &polys[it->p]; }
    const Poly* lookup(int u, int z) const { if (u < 0) return nullptr; int xe = extremalize(u, z); if (xe < 0) return nullptr; return lookup_ext(xe, z); }
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
    void finish(int y, vector<Entry>& res, ll& npairs) {
        auto& t = table[y]; t.swap(res);
        npairs += t.size();
        vector<pair<int32_t, ll>> mu;
        for (int z : M.coatoms[y]) mu.push_back({z, 1});
        for (auto& e : t) { int cod = M.len[y] - M.len[e.x]; if (cod <= 1 || cod % 2 == 0) continue; const Poly& p = polys[e.p]; int d = (cod - 1) / 2; if (p[d]) mu.push_back({e.x, p[d]}); }
        sort(mu.begin(), mu.end()); mu.erase(unique(mu.begin(), mu.end()), mu.end());
        for (size_t i = 1; i < mu.size(); i++) if (mu[i].first == mu[i - 1].first) { fprintf(stderr, "conflicting mu values\n"); exit(1); }
        mulist[y] = mu;
    }
    void run(int nthreads) {
        vector<int> order(M.N); iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b) { return M.len[a] < M.len[b]; });
        int maxlen = M.len[order.back()]; ll npairs = 0; size_t pos = 0;
        for (int L = 0; L <= maxlen; L++) {
            size_t start = pos; while (pos < order.size() && M.len[order[pos]] == L) pos++;
            vector<int> layer(order.begin() + start, order.begin() + pos);
            vector<vector<Entry>> results(layer.size());
            atomic<size_t> next{0};
            auto worker = [&]() {
                vector<int32_t> buf; vector<uint32_t> stamps(M.N, 0); uint32_t cur = 0;
                while (true) {
                    size_t i = next.fetch_add(1); if (i >= layer.size()) break;
                    int y = layer[i]; auto& res = results[i]; res.clear();
                    if (M.len[y] == 0) { Poly one{}; one[0] = 1; res.push_back({y, intern(one)}); continue; }
                    int s = __builtin_ctz(M.rdes[y]); int v = M.rmul[y][s];
                    cur++; buf.clear();
                    for (auto& e : table[v]) if (stamps[e.x] != cur) { stamps[e.x] = cur; buf.push_back(e.x); }
                    for (size_t p = 0; p < buf.size(); p++) {
                        int x = buf[p];
                        uint8_t r = M.rdes[v] & M.rdes[x]; while (r) { int g = __builtin_ctz(r); r &= r - 1; int u = M.rmul[x][g]; if (u >= 0 && stamps[u] != cur) { stamps[u] = cur; buf.push_back(u); } }
                        uint8_t l = M.ldes[v] & M.ldes[x]; while (l) { int g = __builtin_ctz(l); l &= l - 1; int u = M.lmul[x][g]; if (u >= 0 && stamps[u] != cur) { stamps[u] = cur; buf.push_back(u); } }
                    }
                    size_t nv = buf.size();
                    for (size_t p = 0; p < nv; p++) { int x = buf[p]; int xs = M.rmul[x][s]; if (xs >= 0 && stamps[xs] != cur) { stamps[xs] = cur; buf.push_back(xs); } }
                    for (int x : buf) if (extremal(x, y)) { Poly P = compute(x, y, s, v); int cod = M.len[y] - M.len[x]; if (deg(P) > (cod - 1) / 2 || deg(P) >= MAXDEG - 1) { fprintf(stderr, "degree bound violated\n"); exit(1); } res.push_back({x, intern(P)}); }
                    sort(res.begin(), res.end(), [](const Entry& a, const Entry& b) { return a.x < b.x; });
                }
            };
            vector<thread> th; for (int t = 0; t < nthreads; t++) th.emplace_back(worker); for (auto& t : th) t.join();
            for (size_t i = 0; i < layer.size(); i++) finish(layer[i], results[i], npairs);
        }
    }
};

static string polyjson(const Poly& p) { int d = KL::deg(p); string s = "["; for (int i = 0; i <= max(d, 0); i++) { if (i) s += ","; s += to_string(p[i]); } return s + "]"; }
static int failures = 0;
#define EXPECT(cond, what) do { if (!(cond)) { fprintf(stderr, "EXPECTATION FAILED: %s\n", what); failures++; } } while (0)

struct Result { int ideal, lx, lw, interval; vector<ll> rank_ideal, rank_interval; Poly Px, Pe; ll mu; int mismatches; size_t polys; ll maxc; string canon_top, canon_x; };

static Result run_model(const string& model, int n, int nth) {
    vector<int> topword, xword;  // Gern labels 1..n
    auto br = [&](int j, int i) { vector<int> w; if (i == 0) for (int k = j; k >= 1; k--) w.push_back(k); else for (int k = j; k >= i; k--) w.push_back(k); return w; };
    { int k = n / 2 - 2; for (int j = 2; j <= n; j += 2) { auto b = br(j, 0); topword.insert(topword.end(), b.begin(), b.end()); }
      for (int i = 0; i <= k; i++) { auto b = br(n - k + i, n - 2 * k + 2 * i); topword.insert(topword.end(), b.begin(), b.end()); } }
    xword = {1, 2}; for (int j = 4; j <= n; j += 2) xword.push_back(j);
    Model M;
    if (model == "sp") {
        sp::Perm top{};
        for (int i = 1; i <= n; i++) top[i - 1] = (i == 1) ? ((n / 2) % 2 ? -1 : 1) : (i % 2 ? i : -(n + 2 - i));  // Corollary 2.2.19
        M = sp::build(n, top, xword);
        sp::Perm e; for (int i = 0; i < n; i++) e[i] = i + 1; sp::Perm w = e; for (int g : topword) w = sp::rmul(w, g - 1);
        if (w != top) { fprintf(stderr, "Lemma 2.3.4 word and Corollary 2.2.19 disagree\n"); exit(1); }
        if (sp::length(top) != (3 * n * n + 2 * n) / 8) { fprintf(stderr, "length of w_n is not 3n^2/8 + n/4\n"); exit(1); }
    } else M = geo::build(n, topword, xword);
    if ((int)topword.size() != M.len[M.id_top]) { fprintf(stderr, "bracket word not reduced\n"); exit(1); }
    canonical_words(M);
    KL kl(M); kl.run(nth);
    int top = M.id_top; uint8_t rd = M.rdes[top];
    Result R; R.ideal = M.N; R.lx = M.len[M.id_x]; R.lw = M.len[top]; R.canon_top = M.canon[top]; R.canon_x = M.canon[M.id_x];
    R.rank_ideal.assign(R.lw + 1, 0); R.rank_interval.assign(R.lw - R.lx + 1, 0); R.mismatches = 0;
    for (int x = 0; x < M.N; x++) {
        Poly ref{}; bool first = true;
        for (int s = 0; s < n; s++) if ((rd >> s) & 1) { Poly P = kl.compute(x, top, s, M.rmul[top][s]); if (first) { ref = P; first = false; } else if (P != ref) R.mismatches++; }
        Poly extP{}; if (const Poly* e = kl.lookup(x, top)) extP = *e;
        if (extP != ref) R.mismatches++;
        R.rank_ideal[M.len[x]]++;
    }
    R.interval = 0;
    for (int z = 0; z < M.N; z++) if (kl.lookup(M.id_x, z)) { R.rank_interval[M.len[z] - R.lx]++; R.interval++; }
    R.Px = *kl.lookup(M.id_x, top);
    R.Pe = kl.compute(M.id_e, top, __builtin_ctz(rd), M.rmul[top][__builtin_ctz(rd)]);
    int cod = R.lw - R.lx; R.mu = (cod % 2 == 0) ? 0 : R.Px[(cod - 1) / 2];
    R.polys = kl.polys.size(); R.maxc = 0;
    for (size_t j = 0; j < kl.polys.size(); j++) { const Poly& p = kl.polys[j]; if (p[0] != 1) { fprintf(stderr, "constant term != 1\n"); exit(1); } for (int i = 0; i < MAXDEG; i++) { if (p[i] < 0) { fprintf(stderr, "negative coefficient\n"); exit(1); } R.maxc = max(R.maxc, p[i]); } }
    fprintf(stderr, "D%d %s: ideal %d, l(w)=%d, l(x)=%d, interval %d, P(x,w)=%s, mismatches %d\n", n, model.c_str(), R.ideal, R.lw, R.lx, R.interval, polyjson(R.Px).c_str(), R.mismatches);
    return R;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: d8_gern_kl <output.json> [threads] [sp|geo|both]\n"); return 2; }
    int nth = argc > 2 && atoi(argv[2]) > 0 ? atoi(argv[2]) : (int)min<unsigned>(8, max<unsigned>(1, thread::hardware_concurrency()));
    string which = argc > 3 ? argv[3] : "both";
    FILE* f = fopen(argv[1], "w"); if (!f) { fprintf(stderr, "cannot open output\n"); return 2; }
    fprintf(f, "{\n  \"group_conventions\": \"W(D_n), Gern's labels 1..n: s_1 acts on the right by (w_1,w_2)->(-w_2,-w_1), s_i (i>=2) swaps positions i-1,i; nodes 1 and 2 attached to 3, chain 3-...-n\",\n");
    fprintf(f, "  \"models\": [\"signed permutations\", \"geometric representation on the root lattice\"],\n  \"ranks\": {\n");
    vector<ll> expectedD6 = {1, 6, 11, 6, 1, 1};
    vector<ll> expectedD8 = {1, 12, 59, 154, 233, 221, 147, 70, 20, 2};
    for (int n : {6, 8}) {
        Result a = run_model(which == "geo" ? "geo" : "sp", n, nth);
        Result b = which == "both" ? run_model("geo", n, nth) : a;
        EXPECT(a.Px == b.Px && a.Pe == b.Pe && a.ideal == b.ideal && a.interval == b.interval && a.rank_ideal == b.rank_ideal && a.rank_interval == b.rank_interval, "the two models agree");
        EXPECT(a.mismatches == 0 && b.mismatches == 0, "direct recurrence for every right descent agrees with the stored table");
        vector<ll> px(a.Px.begin(), a.Px.begin() + KL::deg(a.Px) + 1), pe(a.Pe.begin(), a.Pe.begin() + KL::deg(a.Pe) + 1);
        if (n == 6) {
            EXPECT(px == expectedD6, "P(x_6,w_6) = 1+6q+11q^2+6q^3+q^4+q^5");
            EXPECT(a.lw == 15 && a.lx == 4 && a.ideal == 3184 && a.interval == 1676 && a.mu == 1, "D6: l(w)=15, l(x)=4, ideal 3184, interval 1676, mu = 1");
            EXPECT(a.rank_interval == vector<ll>({1, 12, 55, 140, 248, 339, 360, 287, 162, 59, 12, 1}), "D6 interval rank vector");
        } else {
            EXPECT(px == expectedD8 && pe == expectedD8, "P(x_8,w_8) = P(e,w_8) = 1+12q+59q^2+154q^3+233q^4+221q^5+147q^6+70q^7+20q^8+2q^9");
            EXPECT(a.lw == 26 && a.lx == 5 && a.ideal == 265760 && a.interval == 163724 && a.mu == 0, "D8: l(w)=26, l(x)=5, ideal 265,760, interval 163,724, mu = 0");
        }
        fprintf(f, "    \"D%d\": {\n", n);
        fprintf(f, "      \"w_canonical_word\": \"%s\", \"x_canonical_word\": \"%s\",\n", a.canon_top.c_str(), a.canon_x.c_str());
        fprintf(f, "      \"w_length\": %d, \"x_length\": %d, \"lower_ideal_size\": %d, \"interval_size\": %d,\n", a.lw, a.lx, a.ideal, a.interval);
        fprintf(f, "      \"lower_ideal_rank_vector\": ["); for (size_t i = 0; i < a.rank_ideal.size(); i++) fprintf(f, "%s%lld", i ? "," : "", a.rank_ideal[i]); fprintf(f, "],\n");
        fprintf(f, "      \"interval_rank_vector\": ["); for (size_t i = 0; i < a.rank_interval.size(); i++) fprintf(f, "%s%lld", i ? "," : "", a.rank_interval[i]); fprintf(f, "],\n");
        fprintf(f, "      \"P_x_w_ascending\": %s,\n      \"P_e_w_ascending\": %s,\n      \"degree\": %d,\n      \"mu_x_w\": %lld,\n", polyjson(a.Px).c_str(), polyjson(a.Pe).c_str(), KL::deg(a.Px), a.mu);
        fprintf(f, "      \"distinct_polynomials\": %zu, \"max_coefficient\": %lld,\n      \"models_agree\": true, \"direct_recurrence_mismatches\": %d\n    }%s\n", a.polys, a.maxc, a.mismatches + b.mismatches, n == 6 ? "," : "");
    }
    fprintf(f, "  },\n  \"x_n\": \"s_1 s_2 s_4 s_6 ... s_n\",\n  \"w_n\": \"Gern Corollary 2.2.19, equal to the bracket word of Lemma 2.3.4\",\n  \"status\": \"%s\"\n}\n", failures ? "FAILED" : "passed");
    fclose(f);
    printf("D6 and D8 Gern-element Kazhdan-Lusztig certificate: %s\n", failures ? "FAILED" : "passed");
    return failures ? 1 : 0;
}
