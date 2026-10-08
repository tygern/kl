// Complete Kazhdan-Lusztig table of a finite simply laced Weyl group of rank
// at most 8, with every leading coefficient mu(x,w), in exact integer
// arithmetic.  Written for the E6 certificate of the manuscript's Section on
// complements: the maximum of mu over all pairs, the pairs attaining it, the
// histogram of mu >= 2, and the fact that every pair with a fully commutative
// lower endpoint has mu in {0,1}.
//
// Method.  Elements are integer matrices in the simple-root basis (column j is
// w(alpha_j)).  The group is enumerated by breadth-first right multiplication,
// Bruhat order is stored as bitsets, and the Kazhdan-Lusztig left recursion
// (Kazhdan-Lusztig 1979, (2.2.c)) is evaluated on extremal pairs only
//   x extremal for w  <=>  x <= w, L(w) subset L(x), R(w) subset R(x),
// using P_{x,w} = P_{sx,w} for s in L(w) with sx > x (and the right version)
// to reach an extremal pair.  Every polynomial is recomputed a second time
// with the right recursion (w = v t, t the largest right descent) and the two
// values must agree.  Constant terms, nonnegativity and the degree bound are
// asserted for every pair.  mu(x,w) for a Bruhat cover x of w is 1 by
// definition; covers are found as w t over all reflections t.
//
// Build: c++ -std=c++17 -O3 e6_mu_table.cpp -o e6_mu_table
// Run:   ./e6_mu_table E 6 <output.json>
// The certificate written for E6 is results/e6-mu-table.json.  The program
// exits with status 1 if any expected E6 value fails.  Diagram labelling for
// type E: chain 0-1-...-(n-2), node n-1 attached to node 2 (the manuscript's
// convention); type D: chain 0-...-(n-2), node n-1 attached to node n-3;
// type A: chain 0-...-(n-1).

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
using namespace std;
typedef long long ll;
typedef unsigned int u32;

static int n;
static int A[8][8];   // Cartan matrix, A[i][j] = <alpha_i^vee, alpha_j>
static int M[8][8];   // Coxeter matrix

struct Mat { int8_t c[8][8]; };  // c[j][k] = coefficient of alpha_k in w(alpha_j)
struct MatHash {
    size_t operator()(const Mat& m) const {
        uint64_t h = 1469598103934665603ULL;
        const unsigned char* p = (const unsigned char*)m.c;
        for (int i = 0; i < 64; i++) { h ^= p[i]; h *= 1099511628211ULL; }
        return (size_t)h;
    }
};
struct MatEq { bool operator()(const Mat& a, const Mat& b) const { return memcmp(a.c, b.c, 64) == 0; } };

static int8_t chk(int v) {
    if (v > 127 || v < -127) { fprintf(stderr, "int8 overflow\n"); exit(3); }
    return (int8_t)v;
}
static Mat rightMul(const Mat& m, int i) {
    Mat r = m;
    for (int j = 0; j < n; j++) if (j != i) {
        int a = A[i][j];
        if (a) for (int k = 0; k < n; k++) r.c[j][k] = chk(m.c[j][k] - a * m.c[i][k]);
    }
    for (int k = 0; k < n; k++) r.c[i][k] = chk(-m.c[i][k]);
    return r;
}
static Mat leftMul(const Mat& m, int i) {
    Mat r = m;
    for (int j = 0; j < n; j++) {
        int pair = 0;
        for (int k = 0; k < n; k++) pair += A[i][k] * m.c[j][k];
        r.c[j][i] = chk(m.c[j][i] - pair);
    }
    return r;
}
static bool colNegative(const Mat& m, int j) {
    for (int k = 0; k < n; k++) if (m.c[j][k] != 0) return m.c[j][k] < 0;
    return false;
}

static u32 N;
static vector<Mat> mats;
static vector<uint8_t> len;
static vector<array<u32, 8>> rmul, lmul;
static vector<uint8_t> Lmask, Rmask;
static vector<vector<uint64_t>> bruhat;

static inline bool leq(u32 x, u32 w) { return (bruhat[w][x >> 6] >> (x & 63)) & 1ULL; }
static bool leqIter(u32 x, u32 w) {
    // Independent Bruhat test by the lifting property, used to cross-check the bitsets.
    while (true) {
        if (x == w) return true;
        if (len[x] >= len[w]) return false;
        int s = __builtin_ctz(Lmask[w]);
        w = lmul[w][s];
        if ((Lmask[x] >> s) & 1) x = lmul[x][s];
    }
}

// Polynomial table with hash-consing.
struct VecHash { size_t operator()(const vector<ll>& v) const { size_t h = v.size(); for (ll a : v) h = h * 1000003ULL ^ (size_t)a ^ (h >> 17); return h; } };
static vector<vector<ll>> polys;
static unordered_map<vector<ll>, u32, VecHash> polyIndex;
static u32 internPoly(vector<ll> p) {
    while (!p.empty() && p.back() == 0) p.pop_back();
    auto it = polyIndex.find(p);
    if (it != polyIndex.end()) return it->second;
    u32 id = polys.size(); polys.push_back(p); polyIndex.emplace(p, id); return id;
}

static vector<vector<u32>> ex;        // sorted extremal x for each w
static vector<vector<u32>> pid;       // polynomial ids parallel to ex
static vector<vector<pair<u32, ll>>> mulist;  // (z, mu(z,w)) with z < w and mu != 0

static inline u32 climb(u32 x, u32 v) {
    while (true) {
        uint8_t dl = Lmask[v] & ~Lmask[x];
        uint8_t dr = Rmask[v] & ~Rmask[x];
        if (!dl && !dr) return x;
        if (dl) x = lmul[x][__builtin_ctz(dl)];
        else x = rmul[x][__builtin_ctz(dr)];
    }
}
static inline const vector<ll>* lookup(u32 x, u32 v) {
    if (!leq(x, v)) return nullptr;
    x = climb(x, v);
    const vector<u32>& e = ex[v];
    auto it = lower_bound(e.begin(), e.end(), x);
    if (it == e.end() || *it != x) { fprintf(stderr, "lookup failure\n"); exit(3); }
    return &polys[pid[v][it - e.begin()]];
}

static vector<int> reducedWord(u32 x) {
    // Canonical reduced word: repeatedly strip the smallest right descent.
    vector<int> w;
    while (x != 0) { int s = __builtin_ctz(Rmask[x]); w.push_back(s); x = rmul[x][s]; }
    reverse(w.begin(), w.end());
    return w;
}
static string wordStr(const vector<int>& w) { string s; for (int a : w) s += to_string(a); return s.empty() ? "e" : s; }
static string jsonList(const vector<int>& v) { string s = "["; for (size_t i = 0; i < v.size(); i++) { if (i) s += ","; s += to_string(v[i]); } return s + "]"; }
static string jsonPoly(const vector<ll>& p) { string s = "["; for (size_t i = 0; i < p.size(); i++) { if (i) s += ","; s += to_string(p[i]); } return s + "]"; }
static vector<int> maskList(uint8_t m) { vector<int> v; for (int i = 0; i < n; i++) if ((m >> i) & 1) v.push_back(i); return v; }

static u32 fromWord(const string& w) { u32 e = 0; for (char c : w) e = rmul[e][c - '0']; return e; }

static int failures = 0;
#define EXPECT(cond, what) do { if (!(cond)) { fprintf(stderr, "EXPECTATION FAILED: %s\n", what); failures++; } } while (0)

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: e6_mu_table <A|D|E> <rank> <output.json>\n"); return 2; }
    char type = argv[1][0]; n = atoi(argv[2]); const char* outfile = argv[3];
    if (n < 2 || n > 8) { fprintf(stderr, "rank must be between 2 and 8\n"); return 2; }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) M[i][j] = (i == j) ? 1 : 2;
    auto edge = [&](int i, int j) { M[i][j] = M[j][i] = 3; };
    if (type == 'A') { for (int i = 0; i + 1 < n; i++) edge(i, i + 1); }
    else if (type == 'D') { for (int i = 0; i + 1 < n - 1; i++) edge(i, i + 1); edge(n - 3, n - 1); }
    else if (type == 'E') { for (int i = 0; i + 1 < n - 1; i++) edge(i, i + 1); edge(2, n - 1); }
    else { fprintf(stderr, "type must be A, D or E\n"); return 2; }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) A[i][j] = (i == j) ? 2 : (M[i][j] == 3 ? -1 : 0);

    // Enumerate the group.
    unordered_map<Mat, u32, MatHash, MatEq> index;
    Mat id; memset(id.c, 0, sizeof id.c); for (int i = 0; i < n; i++) id.c[i][i] = 1;
    mats.push_back(id); len.push_back(0); index.emplace(id, 0);
    rmul.push_back({}); lmul.push_back({});
    for (u32 w = 0; w < mats.size(); w++) {
        Mat m = mats[w];
        for (int i = 0; i < n; i++) {
            Mat r = rightMul(m, i);
            auto it = index.find(r);
            if (it == index.end()) {
                u32 id2 = mats.size(); mats.push_back(r); len.push_back(len[w] + 1); index.emplace(r, id2);
                rmul.push_back({}); lmul.push_back({});
                rmul[w][i] = id2;
            } else rmul[w][i] = it->second;
        }
    }
    N = mats.size();
    for (u32 w = 0; w < N; w++) for (int i = 0; i < n; i++) {
        auto it = index.find(leftMul(mats[w], i));
        if (it == index.end()) { fprintf(stderr, "left multiplication not closed\n"); return 3; }
        lmul[w][i] = it->second;
    }
    index.clear();
    Lmask.assign(N, 0); Rmask.assign(N, 0);
    for (u32 w = 0; w < N; w++) for (int i = 0; i < n; i++) {
        if (colNegative(mats[w], i)) Rmask[w] |= (1 << i);
        if (len[lmul[w][i]] < len[w]) Lmask[w] |= (1 << i);
    }
    int maxLen = len[N - 1];
    fprintf(stderr, "type %c%d: |W| = %u, longest length %d\n", type, n, N, maxLen);

    // Fully commutative elements: R(x) commuting and xs fully commutative for every s in R(x).
    vector<char> fc(N, 0); u32 fcCount = 0;
    for (u32 x = 0; x < N; x++) {
        bool ok = true; uint8_t r = Rmask[x];
        for (int i = 0; i < n && ok; i++) if ((r >> i) & 1) for (int j = i + 1; j < n; j++) if (((r >> j) & 1) && M[i][j] > 2) { ok = false; break; }
        for (int i = 0; i < n && ok; i++) if (((r >> i) & 1) && !fc[rmul[x][i]]) ok = false;
        fc[x] = ok; fcCount += ok;
    }

    // Reflections: one word per positive root.
    vector<vector<int>> reflWords;
    {
        map<vector<int>, vector<int>> rootWord;
        deque<vector<int>> q;
        for (int i = 0; i < n; i++) { vector<int> r(n, 0); r[i] = 1; rootWord[r] = {i}; q.push_back(r); }
        while (!q.empty()) {
            vector<int> r = q.front(); q.pop_front(); vector<int> wd = rootWord[r];
            for (int i = 0; i < n; i++) {
                int pair = 0; for (int k = 0; k < n; k++) pair += A[i][k] * r[k];
                vector<int> r2 = r; r2[i] -= pair;
                bool pos = true; for (int k = 0; k < n; k++) if (r2[k] < 0) pos = false;
                if (!pos || rootWord.count(r2)) continue;
                vector<int> wd2; wd2.push_back(i); wd2.insert(wd2.end(), wd.begin(), wd.end()); wd2.push_back(i);
                rootWord[r2] = wd2; q.push_back(r2);
            }
        }
        for (auto& kv : rootWord) reflWords.push_back(kv.second);
    }
    vector<u32> reflElt;
    for (auto& wd : reflWords) { u32 e = 0; for (int s : wd) e = rmul[e][s]; reflElt.push_back(e); }
    { set<u32> distinct(reflElt.begin(), reflElt.end()); if (distinct.size() != reflWords.size()) { fprintf(stderr, "reflections not distinct\n"); return 3; } }

    // Bruhat order as bitsets: [e,w] = [e,v] union s[e,v] for w = sv, s in L(w).
    size_t words = (N + 63) / 64;
    bruhat.assign(N, vector<uint64_t>(words, 0));
    bruhat[0][0] = 1;
    for (u32 w = 1; w < N; w++) {
        int s = __builtin_ctz(Lmask[w]); u32 v = lmul[w][s];
        bruhat[w] = bruhat[v];
        for (size_t i = 0; i < words; i++) {
            uint64_t b = bruhat[v][i];
            while (b) { int t = __builtin_ctzll(b); b &= b - 1; u32 sx = lmul[i * 64 + t][s]; bruhat[w][sx >> 6] |= 1ULL << (sx & 63); }
        }
    }
    unsigned long long comparable = 0;
    for (u32 w = 0; w < N; w++) for (auto b : bruhat[w]) comparable += __builtin_popcountll(b);
    // Deterministic cross-check of the bitsets against the lifting-property test.
    {
        uint64_t state = 88172645463325252ULL; unsigned long long checked = 0;
        for (int t = 0; t < 200000; t++) {
            state ^= state << 13; state ^= state >> 7; state ^= state << 17;
            u32 x = (u32)(state % N); state ^= state << 13; state ^= state >> 7; state ^= state << 17;
            u32 w = (u32)(state % N);
            if (leq(x, w) != leqIter(x, w)) { fprintf(stderr, "Bruhat bitset mismatch\n"); return 3; }
            checked++;
        }
        fprintf(stderr, "Bruhat bitsets agree with the lifting test on %llu sampled pairs\n", checked);
    }

    // Kazhdan-Lusztig recursion, by increasing length.
    ex.resize(N); pid.resize(N); mulist.resize(N);
    internPoly(vector<ll>{1});
    ex[0] = {0}; pid[0] = {0};
    vector<u32> lenStart(maxLen + 2, 0);
    for (u32 w = 0; w < N; w++) lenStart[len[w] + 1] = w + 1;
    for (int l = 1; l <= maxLen + 1; l++) if (lenStart[l] == 0) lenStart[l] = lenStart[l - 1];
    unsigned long long extremalPairs = 1;
    vector<u32> cand;
    for (int l = 1; l <= maxLen; l++) {
        for (u32 w = lenStart[l]; w < lenStart[l + 1]; w++) {
            int s = __builtin_ctz(Lmask[w]); u32 v = lmul[w][s];
            uint8_t Lw = Lmask[w], Rw = Rmask[w];
            cand.clear();
            for (size_t i = 0; i < words; i++) {
                uint64_t bb = bruhat[w][i];
                while (bb) { int t = __builtin_ctzll(bb); bb &= bb - 1; u32 x = i * 64 + t;
                    if ((Lmask[x] & Lw) == Lw && (Rmask[x] & Rw) == Rw) cand.push_back(x); }
            }
            vector<u32> ids; ids.reserve(cand.size());
            vector<pair<u32, ll>> mus;
            vector<ll> P, Q;
            for (u32 x : cand) {
                if (x == w) { ids.push_back(0); continue; }
                int d = len[w] - len[x];
                P.assign(d / 2 + 2, 0);
                // Left recursion with s in L(w), w = s v; s in L(x) because x is extremal.
                const vector<ll>* p1 = lookup(lmul[x][s], v);
                if (!p1) { fprintf(stderr, "sx not below v\n"); return 3; }
                for (size_t i = 0; i < p1->size(); i++) P[i] += (*p1)[i];
                const vector<ll>* p2 = lookup(x, v);
                if (p2) for (size_t i = 0; i < p2->size(); i++) P[i + 1] += (*p2)[i];
                for (auto& zm : mulist[v]) {
                    u32 z = zm.first;
                    if (!((Lmask[z] >> s) & 1) || len[z] < len[x]) continue;
                    const vector<ll>* p3 = lookup(x, z);
                    if (!p3) continue;
                    int shift = (len[w] - len[z]) / 2;
                    for (size_t i = 0; i < p3->size(); i++) P[i + shift] -= zm.second * (*p3)[i];
                }
                if (P[0] != 1) { fprintf(stderr, "constant term != 1\n"); return 3; }
                for (ll c : P) if (c < 0) { fprintf(stderr, "negative coefficient\n"); return 3; }
                for (size_t i = (d - 1) / 2 + 1; i < P.size(); i++) if (P[i] != 0) { fprintf(stderr, "degree bound violated\n"); return 3; }
                // Right recursion with t the largest right descent of w, w = v' t; t in R(x).
                {
                    int t = 31 - __builtin_clz((unsigned)Rmask[w]); u32 v2 = rmul[w][t];
                    Q.assign(d / 2 + 2, 0);
                    const vector<ll>* q1 = lookup(rmul[x][t], v2);
                    if (!q1) { fprintf(stderr, "xt not below v'\n"); return 3; }
                    for (size_t i = 0; i < q1->size(); i++) Q[i] += (*q1)[i];
                    const vector<ll>* q2 = lookup(x, v2);
                    if (q2) for (size_t i = 0; i < q2->size(); i++) Q[i + 1] += (*q2)[i];
                    for (auto& zm : mulist[v2]) {
                        u32 z = zm.first;
                        if (!((Rmask[z] >> t) & 1) || len[z] < len[x]) continue;
                        const vector<ll>* q3 = lookup(x, z);
                        if (!q3) continue;
                        int shift = (len[w] - len[z]) / 2;
                        for (size_t i = 0; i < q3->size(); i++) Q[i + shift] -= zm.second * (*q3)[i];
                    }
                    if (P != Q) { fprintf(stderr, "left and right recursions disagree\n"); return 3; }
                }
                if (d % 2 == 1) { ll mu = P[(d - 1) / 2]; if (mu) mus.push_back({x, mu}); }
                ids.push_back(internPoly(P));
            }
            // Bruhat covers that are not extremal have mu = 1.
            for (size_t ti = 0; ti < reflElt.size(); ti++) {
                u32 z = w; for (int si : reflWords[ti]) z = rmul[z][si];
                if (len[z] != len[w] - 1) continue;
                if ((Lmask[z] & Lw) == Lw && (Rmask[z] & Rw) == Rw) continue;
                mus.push_back({z, 1});
            }
            sort(mus.begin(), mus.end());
            ex[w] = cand; pid[w] = ids; mulist[w] = mus;
            extremalPairs += cand.size();
        }
    }
    ll maxCoeff = 0;
    for (auto& p : polys) for (ll c : p) maxCoeff = max(maxCoeff, c);

    // Statistics.
    map<ll, unsigned long long> muHist;
    unsigned long long coverPairs = 0, fcLowerNonzero = 0, fcUpperNonzero = 0;
    ll maxMu = 0, maxMuFcLower = 0, maxMuFcUpper = 0;
    vector<tuple<u32, u32, ll>> big;  // mu >= 2
    for (u32 w = 0; w < N; w++) for (auto& zm : mulist[w]) {
        muHist[zm.second]++;
        if (len[w] - len[zm.first] == 1) coverPairs++;
        maxMu = max(maxMu, zm.second);
        if (zm.second >= 2) big.emplace_back(zm.first, w, zm.second);
        if (fc[zm.first]) { fcLowerNonzero++; maxMuFcLower = max(maxMuFcLower, zm.second); }
        if (fc[w]) { fcUpperNonzero++; maxMuFcUpper = max(maxMuFcUpper, zm.second); }
    }
    auto order = [&](const tuple<u32, u32, ll>& a, const tuple<u32, u32, ll>& b) {
        if (get<2>(a) != get<2>(b)) return get<2>(a) > get<2>(b);
        if (len[get<1>(a)] != len[get<1>(b)]) return len[get<1>(a)] < len[get<1>(b)];
        string wa = wordStr(reducedWord(get<1>(a))), wb = wordStr(reducedWord(get<1>(b)));
        if (wa != wb) return wa < wb;
        return wordStr(reducedWord(get<0>(a))) < wordStr(reducedWord(get<0>(b)));
    };
    sort(big.begin(), big.end(), order);
    unsigned long long bigFcUpper = 0; for (auto& t : big) if (fc[get<1>(t)]) bigFcUpper++;

    auto pairJson = [&](u32 x, u32 w, ll mu) {
        const vector<ll>* p = lookup(x, w);
        char buf[64];
        string s = "    {\"mu\": " + to_string(mu);
        s += ", \"x\": \"" + wordStr(reducedWord(x)) + "\", \"x_length\": " + to_string(len[x]) + ", \"x_fully_commutative\": " + (fc[x] ? "true" : "false");
        s += ", \"x_L\": " + jsonList(maskList(Lmask[x])) + ", \"x_R\": " + jsonList(maskList(Rmask[x]));
        s += ", \"w\": \"" + wordStr(reducedWord(w)) + "\", \"w_length\": " + to_string(len[w]) + ", \"w_fully_commutative\": " + (fc[w] ? "true" : "false");
        s += ", \"w_L\": " + jsonList(maskList(Lmask[w])) + ", \"w_R\": " + jsonList(maskList(Rmask[w]));
        snprintf(buf, sizeof buf, ", \"gap\": %d", len[w] - len[x]); s += buf;
        s += ", \"P_ascending\": " + jsonPoly(*p) + "}";
        return s;
    };

    // Expected values for E6 (reviewer computation of 7 October 2026, recomputed here).
    bool e6 = (type == 'E' && n == 6);
    if (e6) {
        EXPECT(N == 51840, "|W(E6)| = 51840");
        EXPECT(maxLen == 36, "longest element of E6 has length 36");
        EXPECT(fcCount == 662, "E6 has 662 fully commutative elements");
        EXPECT(reflWords.size() == 36, "E6 has 36 positive roots");
        EXPECT(maxMu == 10, "max mu over all pairs is 10");
        auto hist = [&](ll k) -> unsigned long long { auto it = muHist.find(k); return it == muHist.end() ? 0 : it->second; };
        EXPECT(hist(10) == 8, "8 pairs attain mu = 10");
        EXPECT(hist(2) == 400 && hist(3) == 556 && hist(4) == 310 && hist(5) == 108 && hist(6) == 4, "histogram of mu >= 2 is 2:400 3:556 4:310 5:108 6:4 10:8");
        EXPECT(hist(7) == 0 && hist(8) == 0 && hist(9) == 0, "no pair has mu in {7,8,9}");
        EXPECT(big.size() == 1386, "1386 pairs have mu >= 2");
        EXPECT(fcLowerNonzero == 6431 && maxMuFcLower == 1, "6431 pairs with FC lower endpoint and mu != 0, all with mu = 1");
        u32 x = fromWord("01254012010"), w = fromWord("123012523412301251234012");
        EXPECT(len[x] == 11 && len[w] == 24 && !fc[x] && !fc[w], "witness lengths 11 and 24, both non-FC");
        const vector<ll>* p = lookup(x, w);
        EXPECT(p && *p == vector<ll>({1, 13, 59, 121, 125, 61, 10}), "P(01254012010, 123012523412301251234012) = 1+13q+59q^2+121q^3+125q^4+61q^5+10q^6");
        u32 x5 = fromWord("5343010"), w5 = fromWord("0125342312501234");
        EXPECT(len[x5] == 7 && len[w5] == 16 && !fc[x5] && fc[w5], "FC-upper example: x of length 7 non-FC, w of length 16 FC");
        const vector<ll>* p5 = lookup(x5, w5);
        EXPECT(p5 && *p5 == vector<ll>({1, 8, 22, 20, 5}), "P(5343010, 0125342312501234) = 1+8q+22q^2+20q^3+5q^4, so mu = 5");
        EXPECT(set<int>(reducedWord(w).begin(), reducedWord(w).end()).size() == 6, "witness w has full support");
    }

    FILE* f = fopen(outfile, "w");
    if (!f) { fprintf(stderr, "cannot open %s\n", outfile); return 2; }
    fprintf(f, "{\n");
    fprintf(f, "  \"type\": \"%c%d\",\n", type, n);
    fprintf(f, "  \"labelling\": \"%s\",\n", type == 'E' ? "chain 0-1-...-(n-2), node n-1 attached to node 2; words act on the right, a string of labels is the product of the simple reflections in the order written" : type == 'D' ? "chain 0-...-(n-2), node n-1 attached to node n-3" : "chain 0-...-(n-1)");
    fprintf(f, "  \"group_order\": %u,\n  \"longest_length\": %d,\n  \"positive_roots\": %zu,\n  \"fully_commutative_count\": %u,\n", N, maxLen, reflWords.size(), fcCount);
    fprintf(f, "  \"bruhat_comparable_pairs\": %llu,\n  \"extremal_pairs\": %llu,\n  \"distinct_polynomials\": %zu,\n  \"max_coefficient\": %lld,\n", comparable, extremalPairs, polys.size(), maxCoeff);
    fprintf(f, "  \"left_right_recursions_agree_on_every_extremal_pair\": true,\n");
    fprintf(f, "  \"max_mu\": %lld,\n", maxMu);
    fprintf(f, "  \"nonzero_mu_pairs\": %llu,\n  \"cover_pairs\": %llu,\n", (unsigned long long)[&]{ unsigned long long t = 0; for (auto& kv : muHist) t += kv.second; return t; }(), coverPairs);
    fprintf(f, "  \"mu_histogram\": {");
    { bool first = true; for (auto& kv : muHist) { fprintf(f, "%s\"%lld\": %llu", first ? "" : ", ", kv.first, kv.second); first = false; } }
    fprintf(f, "},\n");
    fprintf(f, "  \"pairs_with_mu_at_least_2\": %zu,\n", big.size());
    fprintf(f, "  \"pairs_with_mu_at_least_2_and_fully_commutative_lower_endpoint\": %llu,\n", (unsigned long long)[&]{ unsigned long long t = 0; for (auto& b : big) if (fc[get<0>(b)]) t++; return t; }());
    fprintf(f, "  \"pairs_with_mu_at_least_2_and_fully_commutative_upper_endpoint\": %llu,\n", bigFcUpper);
    fprintf(f, "  \"fully_commutative_lower_endpoint\": {\"nonzero_mu_pairs\": %llu, \"max_mu\": %lld},\n", fcLowerNonzero, maxMuFcLower);
    fprintf(f, "  \"fully_commutative_upper_endpoint\": {\"nonzero_mu_pairs\": %llu, \"max_mu\": %lld},\n", fcUpperNonzero, maxMuFcUpper);
    fprintf(f, "  \"pairs_attaining_max_mu\": [\n");
    { bool first = true; for (auto& b : big) { if (get<2>(b) != maxMu) break; fprintf(f, "%s%s", first ? "" : ",\n", pairJson(get<0>(b), get<1>(b), get<2>(b)).c_str()); first = false; } }
    fprintf(f, "\n  ],\n");
    fprintf(f, "  \"largest_mu_with_fully_commutative_upper_endpoint\": [\n");
    { bool first = true; for (auto& b : big) { if (get<2>(b) != maxMuFcUpper || !fc[get<1>(b)]) continue; fprintf(f, "%s%s", first ? "" : ",\n", pairJson(get<0>(b), get<1>(b), get<2>(b)).c_str()); first = false; } }
    fprintf(f, "\n  ],\n");
    fprintf(f, "  \"expectations_checked\": %s,\n", e6 ? "\"E6 values of the manuscript\"" : "null");
    fprintf(f, "  \"status\": \"%s\"\n}\n", failures ? "FAILED" : "passed");
    fclose(f);
    printf("%c%d: |W|=%u FC=%u extremal pairs=%llu max mu=%lld; pairs with mu>=2: %zu; FC lower endpoint nonzero pairs %llu with max mu %lld; status %s\n",
           type, n, N, fcCount, extremalPairs, maxMu, big.size(), fcLowerNonzero, maxMuFcLower, failures ? "FAILED" : "passed");
    return failures ? 1 : 0;
}
