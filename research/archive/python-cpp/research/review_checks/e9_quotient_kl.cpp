// Kazhdan-Lusztig polynomials P_{x,w} for a fixed w in a simply laced Coxeter
// group, computed in the parabolic quotient W/W_I for a commuting set
// I subset R(w), in exact integer arithmetic with overflow checks.  Written
// for the affine E8 (= E9) certificate of the manuscript: the length-33
// reflection r_beta, beta = (1,2,3,3,2,2,1,1,2), and its two length-34 right
// extensions w s_0 and w s_1, whose eligible fully commutative lower endpoints
// have odd gaps, so parity gives no information and the polynomials must be
// computed.
//
// Method.  Elements are integer matrices on the span of the simple roots
// (row j is w(alpha_j) in the simple-root basis); right multiplication acts on
// rows, left multiplication on the coordinates.  For I subset R(w) write
// w = w' w_I with w' in W^I (minimal left coset representatives of W/W_I; w_I
// is the longest element of W_I, here a product of commuting generators).  The
// lower interval [e,w] is a union of cosets x' W_I, and Pq(x,y) := P_{x w_I, y w_I}
// for x,y in W^I satisfies the recursion obtained by restricting the ordinary
// left recursion (Kazhdan-Lusztig 1979, (2.2.c)) to the left ideal
// H C'_{w_I}: for y = s u with s in L(y), u in W^I,
//   Pq(x,y) = Pq(sx,u) + q Pq(x,u)                     if sx < x, sx in W^I,
//   Pq(x,y) = (1+q) Pq(x,u)                             if sx notin W^I (sx = xt, t in I),
//   minus  sum_{z in W^I, sz<z or sz notin W^I, z<u} mu(z w_I, u w_I) q^{(l(y)-l(z))/2} Pq(x,z).
// The quotient is built from the suffixes of a reduced word of w' by left
// multiplication, descents and lengths are read off the matrices, and every
// polynomial is checked for constant term 1, nonnegative coefficients and the
// degree bound.  Before the E9 computation the same engine is validated
// against a naive ordinary Kazhdan-Lusztig recursion on random elements of
// A4, D4 and E6 with random commuting I subset R(w).
//
// Build: c++ -std=c++17 -O3 e9_quotient_kl.cpp -o e9_quotient_kl
// Run:   ./e9_quotient_kl certify <output.json>
// Memory for the length-33 element is about 200 MB; total time a few seconds.
// E_n labelling: chain 0-1-...-(n-2), node n-1 attached to node 2.  A string
// of labels denotes the product of the simple reflections in the order
// written; W acts on the left and the matrix of w has columns w(alpha_j).

#include <algorithm>
#include <array>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
using namespace std;
typedef long long ll;
typedef array<int8_t, 81> Mat;  // m[j*9+i] = coefficient of alpha_i in w(alpha_j)
static int N;
static vector<vector<int>> ADJ;

struct MatHash { size_t operator()(const Mat& m) const { uint64_t h = 1469598103934665603ULL; for (int i = 0; i < N * 9; i++) { h ^= (uint8_t)m[i]; h *= 1099511628211ULL; } return h; } };
static Mat identity() { Mat m; m.fill(0); for (int j = 0; j < N; j++) m[j * 9 + j] = 1; return m; }
static inline int8_t chk(int v) { if (v > 120 || v < -120) { fprintf(stderr, "int8 overflow\n"); exit(3); } return (int8_t)v; }
static Mat rmul(const Mat& w, int s) { Mat r = w; for (int i = 0; i < N; i++) r[s * 9 + i] = chk(-w[s * 9 + i]); for (int t : ADJ[s]) for (int i = 0; i < N; i++) r[t * 9 + i] = chk(w[t * 9 + i] + w[s * 9 + i]); return r; }
static Mat lmul(const Mat& w, int s) { Mat r = w; for (int j = 0; j < N; j++) { int p = 2 * w[j * 9 + s]; for (int t : ADJ[s]) p -= w[j * 9 + t]; r[j * 9 + s] = chk(w[j * 9 + s] - p); } return r; }
static bool negcol(const Mat& w, int j) { bool nz = false; for (int i = 0; i < N; i++) { if (w[j * 9 + i] > 0) return false; if (w[j * 9 + i] < 0) nz = true; } if (!nz) { fprintf(stderr, "zero column\n"); exit(3); } return true; }
static vector<int> rword(Mat w) { vector<int> out; Mat e = identity(); while (w != e) { int s = -1; for (int j = 0; j < N; j++) if (negcol(w, j)) { s = j; break; } if (s < 0) { fprintf(stderr, "no descent\n"); exit(3); } w = rmul(w, s); out.push_back(s); } reverse(out.begin(), out.end()); return out; }
static Mat fromword(const vector<int>& wd) { Mat w = identity(); for (int s : wd) w = rmul(w, s); return w; }
static Mat inverse(const Mat& w) { vector<int> wd = rword(w); reverse(wd.begin(), wd.end()); return fromword(wd); }
static int Rmask(const Mat& w) { int m = 0; for (int j = 0; j < N; j++) if (negcol(w, j)) m |= 1 << j; return m; }
static int Lmask(const Mat& w) { return Rmask(inverse(w)); }
static int pairing(const vector<int>& a, const vector<int>& b) { int p = 0; for (int i = 0; i < N; i++) { p += 2 * a[i] * b[i]; for (int t : ADJ[i]) p -= a[i] * b[t]; } return p; }
static Mat reflection(const vector<int>& beta) { Mat m = identity(); for (int j = 0; j < N; j++) { vector<int> aj(N, 0); aj[j] = 1; int p = pairing(aj, beta); for (int i = 0; i < N; i++) m[j * 9 + i] = chk((i == j) - p * beta[i]); } return m; }
// Full commutativity: w is FC iff every element of its right weak lower
// interval has pairwise commuting right descents (Stembridge 1996, Prop. 2.1:
// a reduced word with a factor sts, m(s,t)=3, exists iff some such element has
// two adjacent right descents).
static unordered_map<Mat, bool, MatHash> FCMEMO;
static bool isFC(const Mat& w) { auto it = FCMEMO.find(w); if (it != FCMEMO.end()) return it->second; int R = Rmask(w); bool ok = true; for (int s = 0; s < N && ok; s++) if (R >> s & 1) { for (int t : ADJ[s]) if (R >> t & 1) ok = false; if (ok && !isFC(rmul(w, s))) ok = false; } FCMEMO[w] = ok; return ok; }

typedef vector<ll> Poly;
static void trim(Poly& p) { while (!p.empty() && p.back() == 0) p.pop_back(); }
static inline ll addchk(ll a, ll b) { ll c; if (__builtin_add_overflow(a, b, &c)) { fprintf(stderr, "overflow\n"); exit(3); } return c; }
static void addto(Poly& a, const Poly& b, int shift, ll scale) { if (a.size() < b.size() + shift) a.resize(b.size() + shift, 0); for (size_t i = 0; i < b.size(); i++) { ll t; if (__builtin_mul_overflow(b[i], scale, &t)) { fprintf(stderr, "overflow\n"); exit(3); } a[i + shift] = addchk(a[i + shift], t); } }
static string pjson(const Poly& p) { string s = "["; for (size_t i = 0; i < p.size(); i++) { if (i) s += ","; s += to_string(p[i]); } return s + "]"; }
static string wjson(const vector<int>& w) { string s = "["; for (size_t i = 0; i < w.size(); i++) { if (i) s += ","; s += to_string(w[i]); } return s + "]"; }
static string wstr(const vector<int>& w) { string s; for (int x : w) s += to_string(x); return s; }
static vector<int> masklist(int m) { vector<int> v; for (int i = 0; i < N; i++) if (m >> i & 1) v.push_back(i); return v; }

struct Quot {
    int I; Mat wI; int lwI;
    vector<Mat> el; vector<int> len; vector<int> Lq;
    vector<array<int, 9>> left;  // index of s x in the quotient, -1 if sx notin W^I (sx = x t), -2 if not in the ideal
    unordered_map<Mat, int, MatHash> idx;
    vector<vector<int>> bucket;  // elements by left descent mask, sorted by length
    bool inWI(const Mat& x) { for (int t = 0; t < N; t++) if (I >> t & 1) if (negcol(x, t)) return false; return true; }
    void build(const Mat& w, int Imask) {
        I = Imask;
        wI = identity();
        while (true) { bool moved = false; for (int t = 0; t < N; t++) if (I >> t & 1) if (!negcol(wI, t)) { wI = rmul(wI, t); moved = true; } if (!moved) break; }
        lwI = rword(wI).size();
        if (lwI != __builtin_popcount(I)) { fprintf(stderr, "I is not a commuting set\n"); exit(3); }
        Mat wp = w; for (int t : rword(wI)) wp = rmul(wp, t);
        if (!inWI(wp)) { fprintf(stderr, "I is not contained in R(w)\n"); exit(3); }
        vector<int> wd = rword(wp);
        unordered_map<Mat, int, MatHash> seen; vector<Mat> cur; Mat e = identity(); seen[e] = 0; cur.push_back(e);
        for (int k = (int)wd.size() - 1; k >= 0; k--) { int s = wd[k]; size_t sz = cur.size(); for (size_t i = 0; i < sz; i++) { Mat y = lmul(cur[i], s); if (!inWI(y) || seen.count(y)) continue; seen[y] = cur.size(); cur.push_back(y); } }
        vector<pair<int, int>> order; order.reserve(cur.size());
        for (size_t i = 0; i < cur.size(); i++) order.push_back({(int)rword(cur[i]).size(), (int)i});
        sort(order.begin(), order.end());
        el.resize(cur.size()); len.resize(cur.size());
        for (size_t i = 0; i < order.size(); i++) { el[i] = cur[order[i].second]; len[i] = order[i].first; idx[el[i]] = i; }
        left.resize(el.size()); Lq.assign(el.size(), 0);
        for (size_t i = 0; i < el.size(); i++) for (int s = 0; s < N; s++) {
            Mat y = lmul(el[i], s);
            if (!inWI(y)) { left[i][s] = -1; Lq[i] |= 1 << s; continue; }
            auto it = idx.find(y);
            if (it == idx.end()) { left[i][s] = -2; continue; }
            left[i][s] = it->second;
            if (len[it->second] < len[i]) Lq[i] |= 1 << s;
            else if (len[it->second] != len[i] + 1) { fprintf(stderr, "length mismatch\n"); exit(3); }
        }
        bucket.assign(1 << N, {});
        for (size_t i = 0; i < el.size(); i++) bucket[Lq[i]].push_back(i);
    }
};

struct KLQ {
    Quot& Q;
    vector<Poly> polys; unordered_map<string, int> pid;
    vector<vector<pair<int, int>>> lst; vector<char> done;
    size_t pairs = 0, computed = 0;
    KLQ(Quot& q) : Q(q) { lst.resize(Q.el.size()); done.assign(Q.el.size(), 0); }
    int id(const Poly& p) { string key((const char*)p.data(), p.size() * sizeof(ll)); auto it = pid.find(key); if (it != pid.end()) return it->second; int i = polys.size(); polys.push_back(p); pid[key] = i; return i; }
    const Poly* get(int x, int z) {
        static Poly zero;
        if (Q.len[x] > Q.len[z]) return &zero;
        int Lz = Q.Lq[z];
        while (true) {
            int d = Lz & ~Q.Lq[x]; if (!d) break;
            int t = __builtin_ctz(d); int nx = Q.left[x][t];
            if (nx < 0) { if (nx == -2) return &zero; fprintf(stderr, "impossible reduction\n"); exit(3); }
            x = nx; if (Q.len[x] > Q.len[z]) return &zero;
        }
        auto& v = lst[z]; auto it = lower_bound(v.begin(), v.end(), make_pair(x, INT_MIN));
        if (it == v.end() || it->first != x) return &zero;
        return &polys[it->second];
    }
    void compute(int y) {
        if (done[y]) return;
        if (Q.len[y] == 0) { lst[y] = {{y, id(Poly{1})}}; done[y] = 1; return; }
        int s = -1;
        for (int c = 0; c < N; c++) if ((Q.Lq[y] >> c & 1) && Q.left[y][c] >= 0) { if (s < 0) s = c; if (done[Q.left[y][c]]) { s = c; break; } }
        if (s < 0) { fprintf(stderr, "no quotient descent\n"); exit(3); }
        int u = Q.left[y][s]; compute(u);
        vector<tuple<int, ll, int>> corr;
        for (auto& pr : lst[u]) {
            int z = pr.first; if (z == u || !(Q.Lq[z] >> s & 1)) continue;
            int gap = Q.len[u] - Q.len[z]; if (gap % 2 == 0) continue;
            const Poly& p = polys[pr.second]; size_t top = (gap - 1) / 2;
            if (p.size() > top + 1) { fprintf(stderr, "degree violation\n"); exit(3); }
            if (p.size() == top + 1 && p[top] != 0) corr.emplace_back(z, p[top], (Q.len[y] - Q.len[z]) / 2);
        }
        for (int t = 0; t < N; t++) if ((Q.Lq[u] >> t & 1) && Q.left[u][t] >= 0) { int z = Q.left[u][t]; if (Q.Lq[z] >> s & 1) corr.emplace_back(z, 1, 1); }
        for (auto& c : corr) compute(std::get<0>(c));
        vector<pair<int, int>> out;
        int Ly = Q.Lq[y]; int full = (1 << N) - 1; int freebits = full & ~Ly;
        for (int sub = freebits;; sub = (sub - 1) & freebits) {
            int D = Ly | sub; auto& b = Q.bucket[D];
            for (int x : b) {
                if (Q.len[x] >= Q.len[y]) break;
                Poly p; int sx = Q.left[x][s];
                if (sx == -1) { const Poly* a = get(x, u); addto(p, *a, 0, 1); addto(p, *a, 1, 1); }
                else { if (sx < 0 || Q.len[sx] >= Q.len[x]) { fprintf(stderr, "candidate without s-descent\n"); exit(3); } addto(p, *get(sx, u), 0, 1); addto(p, *get(x, u), 1, 1); }
                for (auto& c : corr) { int z = std::get<0>(c); if (Q.len[x] > Q.len[z]) continue; const Poly* pz = get(x, z); if (!pz->empty()) addto(p, *pz, std::get<2>(c), -std::get<1>(c)); }
                trim(p); if (p.empty()) continue;
                int gap = Q.len[y] - Q.len[x];
                if ((int)p.size() - 1 > (gap - 1) / 2) { fprintf(stderr, "degree bound violated\n"); exit(3); }
                for (ll c : p) if (c < 0) { fprintf(stderr, "negative coefficient\n"); exit(3); }
                if (p[0] != 1) { fprintf(stderr, "constant term != 1\n"); exit(3); }
                out.push_back({x, id(p)});
            }
            if (sub == 0) break;
        }
        out.push_back({y, id(Poly{1})}); sort(out.begin(), out.end()); lst[y] = out; done[y] = 1; pairs += out.size(); computed++;
    }
};

// Naive ordinary Kazhdan-Lusztig recursion on the lower interval of w, for validation.
struct Ord {
    vector<Mat> el; vector<int> len; unordered_map<Mat, int, MatHash> idx; vector<array<int, 9>> left; vector<int> L;
    map<pair<int, int>, Poly> memo;
    void build(const Mat& w) {
        vector<int> wd = rword(w); unordered_map<Mat, int, MatHash> seen; vector<Mat> cur; Mat e = identity(); seen[e] = 0; cur.push_back(e);
        for (int s : wd) { size_t sz = cur.size(); for (size_t i = 0; i < sz; i++) { Mat y = rmul(cur[i], s); if (seen.count(y)) continue; seen[y] = cur.size(); cur.push_back(y); } }
        vector<pair<int, int>> order; for (size_t i = 0; i < cur.size(); i++) order.push_back({(int)rword(cur[i]).size(), (int)i}); sort(order.begin(), order.end());
        el.resize(cur.size()); len.resize(cur.size()); for (size_t i = 0; i < order.size(); i++) { el[i] = cur[order[i].second]; len[i] = order[i].first; idx[el[i]] = i; }
        left.resize(el.size()); L.assign(el.size(), 0);
        for (size_t i = 0; i < el.size(); i++) for (int s = 0; s < N; s++) { Mat y = lmul(el[i], s); auto it = idx.find(y); if (it == idx.end()) { left[i][s] = -2; continue; } left[i][s] = it->second; if (len[it->second] < len[i]) L[i] |= 1 << s; }
    }
    bool leq(int x, int y) { while (x != y) { if (len[x] >= len[y]) return false; int s = __builtin_ctz(L[y]); if (L[x] >> s & 1) x = left[x][s]; y = left[y][s]; } return true; }
    Poly P(int x, int y) {
        if (x == y) return Poly{1};
        if (!leq(x, y)) return Poly{};
        auto key = make_pair(x, y); auto it = memo.find(key); if (it != memo.end()) return it->second;
        int s = __builtin_ctz(L[y]); int u = left[y][s]; Poly p; int c = (L[x] >> s & 1); int sx = left[x][s];
        if (sx == -2) { fprintf(stderr, "ord: sx missing\n"); exit(3); }
        addto(p, P(sx, u), 1 - c, 1); addto(p, P(x, u), c, 1);
        for (int z = 0; z < (int)el.size(); z++) { if (len[z] >= len[u]) break; if (!(L[z] >> s & 1)) continue; int gap = len[u] - len[z]; if (gap % 2 == 0 || !leq(x, z)) continue; Poly pz = P(z, u); size_t top = (gap - 1) / 2; if (pz.size() == top + 1 && pz[top]) addto(p, P(x, z), (len[y] - len[z]) / 2, -pz[top]); }
        trim(p); memo[key] = p; return p;
    }
};

static void setup(const string& type) {
    ADJ.clear();
    auto add = [&](int a, int b) { ADJ[a].push_back(b); ADJ[b].push_back(a); };
    if (type[0] == 'E') { N = stoi(type.substr(1)); ADJ.assign(N, {}); for (int i = 0; i + 1 < N - 1; i++) add(i, i + 1); add(2, N - 1); }
    else if (type[0] == 'D') { N = stoi(type.substr(1)); ADJ.assign(N, {}); for (int i = 0; i + 1 < N - 1; i++) add(i, i + 1); add(N - 3, N - 1); }
    else if (type[0] == 'A') { N = stoi(type.substr(1)); ADJ.assign(N, {}); for (int i = 0; i + 1 < N; i++) add(i, i + 1); }
    else { fprintf(stderr, "unknown type\n"); exit(2); }
    if (N > 9) { fprintf(stderr, "rank at most 9\n"); exit(2); }
}

// Validation: random elements w, random commuting I subset R(w); every quotient
// pair must agree with the naive recursion.  Returns the number of pairs compared.
static long validate(const string& type, unsigned seed, int count, int maxlen) {
    setup(type); mt19937 rng(seed); long pairs = 0; int tested = 0;
    for (int it = 0; it < count; it++) {
        int L = 1 + rng() % maxlen; Mat w = identity();
        for (int i = 0; i < L; i++) { int s = rng() % N; if (!negcol(w, s)) w = rmul(w, s); }
        int R = Rmask(w); if (!R) continue;
        int I = 0;
        for (int s = 0; s < N; s++) if ((R >> s & 1) && (rng() % 2)) { bool adj = false; for (int t : ADJ[s]) if (I >> t & 1) adj = true; if (!adj) I |= 1 << s; }
        if (!I) I = 1 << __builtin_ctz(R);
        Ord o; o.build(w); if (o.el.size() > 6000) continue;
        Quot q; q.build(w, I); KLQ k(q);
        Mat wp = w; for (int t : rword(q.wI)) wp = rmul(wp, t);
        int yq = q.idx[wp]; k.compute(yq);
        int ow = o.idx[w];
        for (size_t xi = 0; xi < q.el.size(); xi++) {
            Mat xm = q.el[xi]; for (int t : rword(q.wI)) xm = rmul(xm, t);
            auto f = o.idx.find(xm); Poly ref = (f == o.idx.end()) ? Poly{} : o.P(f->second, ow);
            if (ref != *k.get(xi, yq)) { fprintf(stderr, "validation mismatch in %s\n", type.c_str()); exit(1); }
            pairs++;
        }
        tested++;
    }
    fprintf(stderr, "validated %s: %d elements, %ld quotient pairs agree with the naive recursion\n", type.c_str(), tested, pairs);
    return pairs;
}

struct Target { string name; vector<int> ext; };
static int failures = 0;
#define EXPECT(cond, what) do { if (!(cond)) { fprintf(stderr, "EXPECTATION FAILED: %s\n", what); failures++; } } while (0)

static string certify_one(const Target& tg, const vector<int>& beta, const vector<int>& baseword) {
    setup("E9");
    Mat w = reflection(beta);
    EXPECT(w == fromword(baseword), "r_beta equals the manuscript's reduced word");
    for (int s : tg.ext) { if (negcol(w, s)) { fprintf(stderr, "extension letter is a descent\n"); exit(3); } w = rmul(w, s); }
    vector<int> wd = rword(w); int lw = wd.size();
    int R = Rmask(w), L = Lmask(w); int I = R;
    Mat winv = inverse(w); bool invol = (winv == w);
    bool term = true;
    for (auto v : {w, winv}) for (int s = 0; s < N; s++) if (negcol(v, s)) { Mat vs = rmul(v, s); for (int t : ADJ[s]) if (negcol(vs, t)) term = false; }
    Quot q; q.build(w, I);
    size_t qsize = q.el.size(); long long ideal = (long long)qsize * (1LL << __builtin_popcount(I));
    // Eligible fully commutative lower endpoints: x <= w, L(w) subset L(x), R(w) subset R(x), x FC.
    vector<int> elig;
    for (size_t xi = 0; xi < qsize; xi++) { Mat xm = q.el[xi]; for (int t : rword(q.wI)) xm = rmul(xm, t); if ((Lmask(xm) & L) != L || (Rmask(xm) & R) != R) continue; if (!isFC(xm)) continue; elig.push_back(xi); }
    Mat wp = w; for (int t : rword(q.wI)) wp = rmul(wp, t); int yq = q.idx[wp];
    KLQ k(q); k.compute(yq);
    string js = "    {\"name\": \"" + tg.name + "\", \"word\": " + wjson(wd) + ", \"word_string\": \"" + wstr(wd) + "\", \"length\": " + to_string(lw);
    js += ", \"L\": " + wjson(masklist(L)) + ", \"R\": " + wjson(masklist(R)) + ", \"involution\": " + (invol ? "true" : "false") + ", \"terminal\": " + (term ? "true" : "false");
    js += ", \"I\": " + wjson(masklist(I)) + ", \"W_I_order\": " + to_string(1 << __builtin_popcount(I)) + ", \"quotient_size\": " + to_string(qsize) + ", \"lower_ideal_size\": " + to_string(ideal);
    js += ", \"quotient_elements_computed\": " + to_string(k.computed) + ", \"stored_pairs\": " + to_string(k.pairs) + ", \"distinct_polynomials\": " + to_string(k.polys.size());
    js += ",\n     \"eligible_fully_commutative_bottoms\": [";
    vector<tuple<string, int, int, ll, Poly>> rows;
    for (int xi : elig) { Mat xm = q.el[xi]; for (int t : rword(q.wI)) xm = rmul(xm, t); auto xw = rword(xm); const Poly* p = k.get(xi, yq); int gap = lw - xw.size(); ll mu = 0; if (gap % 2 == 1) { size_t top = (gap - 1) / 2; if (p->size() == top + 1) mu = (*p)[top]; } rows.emplace_back(wstr(xw), (int)xw.size(), gap, mu, *p); }
    sort(rows.begin(), rows.end());
    for (size_t i = 0; i < rows.size(); i++) { auto& [xs, xl, gap, mu, p] = rows[i]; js += string(i ? ",\n       " : "\n       ") + "{\"x\": \"" + xs + "\", \"length\": " + to_string(xl) + ", \"gap\": " + to_string(gap) + ", \"P_ascending\": " + pjson(p) + ", \"degree\": " + to_string((int)p.size() - 1) + ", \"mu\": " + to_string(mu) + "}"; }
    js += "],\n";
    // Histogram of nonzero mu over all extremal lower endpoints in the quotient.
    ll maxmu = 0; map<ll, int> hist;
    for (auto& pr : k.lst[yq]) { int x = pr.first; if (x == yq) continue; int gap = q.len[yq] - q.len[x]; if (gap % 2 == 0) continue; const Poly& p = k.polys[pr.second]; size_t top = (gap - 1) / 2; if (p.size() == top + 1 && p[top]) { hist[p[top]]++; maxmu = max(maxmu, p[top]); } }
    js += "     \"nonzero_mu_histogram_over_quotient_representatives\": {"; { bool first = true; for (auto& h : hist) { js += (first ? "" : ", ") + string("\"") + to_string(h.first) + "\": " + to_string(h.second); first = false; } } js += "}";
    js += ", \"max_mu_over_quotient_representatives\": " + to_string(maxmu);
    size_t checked = 0, bad = 0, missing = 0;
    if (invol) { for (auto& pr : k.lst[yq]) { Mat X = q.el[pr.first]; for (int t : rword(q.wI)) X = rmul(X, t); Mat Xi = inverse(X); for (int t : rword(q.wI)) Xi = rmul(Xi, t); auto it = q.idx.find(Xi); if (it == q.idx.end()) { missing++; continue; } if (*k.get(it->second, yq) != k.polys[pr.second]) bad++; checked++; } }
    js += ", \"inversion_symmetry\": {\"pairs_checked\": " + to_string(checked) + ", \"mismatches\": " + to_string(bad) + ", \"inverse_outside_quotient\": " + to_string(missing) + "}";
    ll maxc = 0; size_t maxdeg = 0; for (auto& p : k.polys) { for (ll c : p) maxc = max(maxc, c); maxdeg = max(maxdeg, p.size()); }
    js += ", \"max_coefficient\": " + to_string(maxc) + ", \"max_degree\": " + to_string(maxdeg ? maxdeg - 1 : 0) + "}";
    EXPECT(bad == 0, "P(x,w) = P(x^-1,w) for the involution w");
    // Expected values (reviewer computation of 7 October 2026).
    auto expect_row = [&](const string& xs, int gap, const Poly& P, ll mu) {
        bool found = false;
        for (auto& r : rows) if (get<0>(r) == xs) { found = true; EXPECT(get<2>(r) == gap, ("gap of " + xs).c_str()); EXPECT(get<4>(r) == P, ("polynomial of " + xs).c_str()); EXPECT(get<3>(r) == mu, ("mu of " + xs).c_str()); }
        EXPECT(found, ("eligible bottom " + xs + " present").c_str());
    };
    if (tg.name == "w33") {
        EXPECT(lw == 33 && invol && term, "length 33 terminal involution");
        EXPECT(masklist(R) == vector<int>({3, 5, 7, 8}) && masklist(L) == vector<int>({3, 5, 7, 8}), "L = R = {3,5,7,8}");
        EXPECT(qsize == 364156 && ideal == 5826496, "quotient 364,156 cosets, ideal 5,826,496 elements");
        EXPECT(rows.size() == 5, "exactly five eligible FC bottoms");
        expect_row("8753", 29, {1, 23, 236, 1391, 5298, 13861, 25666, 33996, 31954, 20820, 8985, 2380, 346, 23}, 0);
        expect_row("87531", 28, {1, 23, 236, 1389, 5262, 13609, 24701, 31740, 28612, 17622, 7052, 1654, 194, 9}, 0);
        expect_row("87530", 28, {1, 23, 235, 1370, 5128, 13088, 23433, 29667, 26255, 15772, 6082, 1350, 144, 5}, 0);
        expect_row("875301", 27, {1, 23, 234, 1348, 4931, 12136, 20641, 24398, 19796, 10646, 3551, 621, 40, 1}, 1);
        expect_row("875310", 27, {1, 23, 234, 1348, 4931, 12136, 20641, 24398, 19796, 10646, 3551, 621, 40, 1}, 1);
        EXPECT(set<int>(wd.begin(), wd.end()).size() == 9, "full support");
    } else if (tg.name == "w34_s0") {
        EXPECT(lw == 34 && invol && term, "length 34 terminal involution");
        EXPECT(masklist(R) == vector<int>({0, 3, 5, 7, 8}) && masklist(L) == vector<int>({0, 3, 5, 7, 8}), "L = R = {0,3,5,7,8}");
        EXPECT(qsize == 216990 && ideal == 6943680, "quotient 216,990 cosets, ideal 6,943,680 elements");
        EXPECT(rows.size() == 1, "unique eligible FC bottom");
        expect_row("87530", 29, {1, 18, 146, 712, 2342, 5456, 9200, 11264, 9947, 6220, 2649, 721, 116, 10, 1}, 1);
    } else if (tg.name == "w34_s1") {
        EXPECT(lw == 34 && invol && term, "length 34 terminal involution");
        EXPECT(masklist(R) == vector<int>({1, 3, 5, 7, 8}) && masklist(L) == vector<int>({1, 3, 5, 7, 8}), "L = R = {1,3,5,7,8}");
        EXPECT(qsize == 207866 && ideal == 6651712, "quotient 207,866 cosets, ideal 6,651,712 elements");
        EXPECT(rows.size() == 1, "unique eligible FC bottom");
        expect_row("87531", 29, {1, 18, 140, 639, 1957, 4280, 6858, 8114, 7084, 4497, 1996, 579, 100, 10, 1}, 1);
    }
    fprintf(stderr, "%s: length %d, quotient %zu, ideal %lld, %zu eligible FC bottoms\n", tg.name.c_str(), lw, qsize, ideal, rows.size());
    FCMEMO.clear();
    return js;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: e9_quotient_kl certify <output.json> | e9_quotient_kl validate <type> <seed> <count> <maxlen>\n"); return 2; }
    string mode = argv[1];
    if (mode == "validate") { if (argc < 6) return 2; long p = validate(argv[2], atoi(argv[3]), atoi(argv[4]), atoi(argv[5])); printf("validated: %ld pairs\n", p); return 0; }
    if (mode != "certify" || argc < 3) return 2;
    long va = validate("A4", 1, 200, 10), vd = validate("D4", 2, 200, 12), ve = validate("E6", 3, 120, 14);
    vector<int> beta = {1, 2, 3, 3, 2, 2, 1, 1, 2};
    vector<int> baseword = {8, 7, 5, 6, 3, 4, 5, 2, 3, 4, 1, 2, 3, 0, 1, 2, 8, 2, 3, 4, 1, 2, 3, 0, 1, 2, 8, 5, 6, 7, 4, 5, 3};
    vector<Target> targets = {{"w33", {}}, {"w34_s0", {0}}, {"w34_s1", {1}}};
    vector<string> parts;
    for (auto& t : targets) parts.push_back(certify_one(t, beta, baseword));
    FILE* f = fopen(argv[2], "w"); if (!f) { fprintf(stderr, "cannot open output\n"); return 2; }
    fprintf(f, "{\n  \"group\": \"affine E8 = E9, chain 0-1-...-7, node 8 attached to node 2\",\n");
    fprintf(f, "  \"method\": \"Kazhdan-Lusztig recursion in the parabolic quotient W/W_I, I = R(w); exact integers with overflow checks\",\n");
    fprintf(f, "  \"beta\": [1,2,3,3,2,2,1,1,2],\n  \"base_element\": \"r_beta, the manuscript's length-33 reflection\",\n");
    fprintf(f, "  \"engine_validation\": {\"A4_pairs\": %ld, \"D4_pairs\": %ld, \"E6_pairs\": %ld, \"all_equal_to_naive_recursion\": true},\n", va, vd, ve);
    fprintf(f, "  \"elements\": [\n");
    for (size_t i = 0; i < parts.size(); i++) fprintf(f, "%s%s", parts[i].c_str(), i + 1 < parts.size() ? ",\n" : "\n");
    fprintf(f, "  ],\n  \"all_eligible_mu_values_in\": [0, 1],\n  \"status\": \"%s\"\n}\n", failures ? "FAILED" : "passed");
    fclose(f);
    printf("E9 odd-gap certificate: %s\n", failures ? "FAILED" : "passed");
    return failures ? 1 : 0;
}
