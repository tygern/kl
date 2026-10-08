// klmu.cpp -- independent computation of all Kazhdan-Lusztig polynomials and
// all mu-values of a finite Coxeter group, using exact integer arithmetic.
//
// Method: faithful integral geometric representation (matrices in the simple-root
// basis), BFS enumeration of the group, Bruhat order either as bitsets (small
// groups) or via the iterative lifting-property test (large groups), and the
// Kazhdan-Lusztig left recursion restricted to extremal pairs
//   (x extremal for w  <=>  x <= w, L(w) subset L(x), R(w) subset R(x)),
// using P_{x,w} = P_{sx,w} (s in L(w), sx > x) and P_{x,w} = P_{xs,w} to climb.
//
// Recursion (Kazhdan-Lusztig 1979, (2.2.c)), w = s v, l(w) = l(v)+1, s in L(x):
//   P_{x,w} = P_{sx,v} + q P_{x,v} - sum_{z < v, s z < z} mu(z,v) q^{(l(w)-l(z))/2} P_{x,z}
//
// Usage: klmu <type A|B|D|E|F|H> <rank> [bitset|desc] [dump <file>] [threads T]
// E_n numbering follows the manuscript: chain 0-1-...-(n-2), node n-1 attached to node 2.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <array>
#include <string>
#include <map>
#include <set>
#include <deque>
#include <unordered_map>
#include <algorithm>
#include <thread>
#include <mutex>
#include <atomic>
#include <chrono>
#include <random>
#include <tuple>
using namespace std;
typedef long long ll;
typedef unsigned int u32;

static int n;                       // rank
static int A[8][8];                 // Cartan matrix A[i][j] = <alpha_i^vee, alpha_j>
static int M[8][8];                 // Coxeter matrix

struct Mat { int8_t c[8][8]; };     // c[j][k] = coefficient of alpha_k in w(alpha_j)
struct MatHash {
    size_t operator()(const Mat& m) const {
        uint64_t h = 1469598103934665603ULL;
        const unsigned char* p = (const unsigned char*)m.c;
        for (int i = 0; i < 64; i++) { h ^= p[i]; h *= 1099511628211ULL; }
        return (size_t)h;
    }
};
struct MatEq { bool operator()(const Mat& a, const Mat& b) const { return memcmp(a.c, b.c, 64) == 0; } };

static Mat rightMul(const Mat& m, int i) {
    Mat r = m;
    for (int j = 0; j < n; j++) if (j != i) {
        int a = A[i][j];
        if (a) for (int k = 0; k < n; k++) r.c[j][k] = (int8_t)(m.c[j][k] - a * m.c[i][k]);
    }
    for (int k = 0; k < n; k++) r.c[i][k] = (int8_t)(-m.c[i][k]);
    return r;
}
static Mat leftMul(const Mat& m, int i) {
    Mat r = m;
    for (int j = 0; j < n; j++) {
        int pair = 0;
        for (int k = 0; k < n; k++) pair += A[i][k] * m.c[j][k];
        r.c[j][i] = (int8_t)(m.c[j][i] - pair);
    }
    return r;
}
static bool colNegative(const Mat& m, int j) {
    for (int k = 0; k < n; k++) if (m.c[j][k] != 0) return m.c[j][k] < 0;
    return false;
}

// ---------- group data ----------
static u32 N;
static vector<Mat> mats;
static vector<uint8_t> len;
static vector<array<u32,8>> rmul, lmul;
static vector<uint8_t> Lmask, Rmask;

static vector<vector<uint64_t>> bruhat;   // bitset mode
static bool useBitset = true;

static inline bool leqBit(u32 x, u32 w) { return (bruhat[w][x >> 6] >> (x & 63)) & 1ULL; }
static bool leqIter(u32 x, u32 w) {
    while (true) {
        if (x == w) return true;
        if (len[x] >= len[w]) return false;
        int s = __builtin_ctz(Lmask[w]);
        w = lmul[w][s];
        if ((Lmask[x] >> s) & 1) x = lmul[x][s];
    }
}
static inline bool leq(u32 x, u32 w) { return useBitset ? leqBit(x, w) : leqIter(x, w); }

// ---------- polynomial table (chunked, references stay valid on append) ----------
struct VecHash { size_t operator()(const vector<ll>& v) const { size_t h = v.size(); for (ll a : v) h = h * 1000003ULL ^ (size_t)a ^ (h >> 17); return h; } };
static const u32 CHUNK = 1 << 16;
static vector<ll>* polyChunks[1 << 14];
static u32 polyCount = 0;
static unordered_map<vector<ll>, u32, VecHash> polyIndex;
static mutex polyMutex;
static inline const vector<ll>& getPoly(u32 id) { return polyChunks[id / CHUNK][id % CHUNK]; }
static u32 internPoly(vector<ll>& p) {
    while (!p.empty() && p.back() == 0) p.pop_back();
    lock_guard<mutex> g(polyMutex);
    auto it = polyIndex.find(p);
    if (it != polyIndex.end()) return it->second;
    u32 id = polyCount++;
    if (id % CHUNK == 0) polyChunks[id / CHUNK] = new vector<ll>[CHUNK];
    polyChunks[id / CHUNK][id % CHUNK] = p; polyIndex.emplace(p, id); return id;
}

// per-element KL data
static vector<vector<u32>> ex;       // sorted extremal x for w
static vector<vector<u32>> pid;      // polynomial ids, parallel to ex
static vector<vector<pair<u32,ll>>> mulist; // (z, mu(z,w)) with z<w, mu != 0

static inline u32 climb(u32 x, u32 v) {
    // x <= v assumed; return extremal x' with P_{x,v} = P_{x',v}
    while (true) {
        uint8_t dl = Lmask[v] & ~Lmask[x];
        uint8_t dr = Rmask[v] & ~Rmask[x];
        if (!dl && !dr) return x;
        if (dl) x = lmul[x][__builtin_ctz(dl)];
        else x = rmul[x][__builtin_ctz(dr)];
    }
}
static inline const vector<ll>* lookup(u32 x, u32 v) {
    // returns pointer to P_{x,v} or nullptr if x not <= v
    if (!leq(x, v)) return nullptr;
    x = climb(x, v);
    const vector<u32>& e = ex[v];
    auto it = lower_bound(e.begin(), e.end(), x);
    if (it == e.end() || *it != x) { fprintf(stderr, "lookup failure x=%u v=%u\n", x, v); abort(); }
    return &getPoly(pid[v][it - e.begin()]);
}

// reflections (words) for cover enumeration
static vector<vector<int>> reflWords;

static vector<int> reducedWord(u32 x) {
    vector<int> w;
    while (x != 0) { int s = __builtin_ctz(Rmask[x]); w.push_back(s); x = rmul[x][s]; }
    reverse(w.begin(), w.end());
    return w;
}
static string wordStr(const vector<int>& w) { string s; for (int a : w) s += to_string(a); return s.empty() ? "e" : s; }
static string maskStr(uint8_t m) { string s = "{"; for (int i = 0; i < n; i++) if ((m >> i) & 1) { if (s.size() > 1) s += ","; s += to_string(i); } return s + "}"; }

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: klmu <type> <rank> [bitset|desc] [dump file] [threads T]\n"); return 1; }
    char type = argv[1][0]; n = atoi(argv[2]);
    string dumpFile; int threads = 1; bool verify = false;
    for (int i = 3; i < argc; i++) {
        string a = argv[i];
        if (a == "bitset") useBitset = true; else if (a == "desc") useBitset = false;
        else if (a == "dump") dumpFile = argv[++i];
        else if (a == "threads") threads = atoi(argv[++i]);
        else if (a == "verify") verify = true;
    }
    // Coxeter matrix
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) M[i][j] = (i == j) ? 1 : 2;
    auto edge = [&](int i, int j, int m) { M[i][j] = M[j][i] = m; };
    if (type == 'A') { for (int i = 0; i + 1 < n; i++) edge(i, i + 1, 3); }
    else if (type == 'B') { for (int i = 0; i + 1 < n; i++) edge(i, i + 1, 3); edge(n - 2, n - 1, 4); }
    else if (type == 'D') { for (int i = 0; i + 1 < n - 1; i++) edge(i, i + 1, 3); edge(n - 3, n - 1, 3); }
    else if (type == 'E') { for (int i = 0; i + 1 < n - 1; i++) edge(i, i + 1, 3); edge(2, n - 1, 3); }
    else if (type == 'F') { edge(0, 1, 3); edge(1, 2, 4); edge(2, 3, 3); }
    else if (type == 'H') { edge(0, 1, 5); for (int i = 1; i + 1 < n; i++) edge(i, i + 1, 3); }
    else { fprintf(stderr, "unknown type\n"); return 1; }
    // Cartan matrix (crystallographic for A,B,D,E,F; H handled via non-integral -> not supported here)
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) {
        if (i == j) A[i][j] = 2;
        else if (M[i][j] == 2) A[i][j] = 0;
        else if (M[i][j] == 3) A[i][j] = -1;
        else if (M[i][j] == 4) A[i][j] = 0; // set below
        else { fprintf(stderr, "non-crystallographic not supported\n"); return 1; }
    }
    if (type == 'B') { A[n - 2][n - 1] = -2; A[n - 1][n - 2] = -1; }   // alpha_{n-1} short
    if (type == 'F') { A[1][2] = -2; A[2][1] = -1; }

    // ---- enumerate group by BFS on right multiplication ----
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
        Mat l = leftMul(mats[w], i);
        auto it = index.find(l); if (it == index.end()) { fprintf(stderr, "left mult not closed\n"); return 1; }
        lmul[w][i] = it->second;
    }
    index.clear(); index.rehash(0);
    Lmask.assign(N, 0); Rmask.assign(N, 0);
    for (u32 w = 0; w < N; w++) for (int i = 0; i < n; i++) {
        if (colNegative(mats[w], i)) Rmask[w] |= (1 << i);
        if (len[lmul[w][i]] < len[w]) Lmask[w] |= (1 << i);
    }
    int maxLen = len[N - 1];
    fprintf(stderr, "type %c%d: |W| = %u, longest length %d\n", type, n, N, maxLen);

    // FC elements: fc[x] = R(x) independent and fc[xs] for all s in R(x)
    vector<char> fc(N, 0); u32 fcCount = 0;
    for (u32 x = 0; x < N; x++) {
        bool ok = true; uint8_t r = Rmask[x];
        for (int i = 0; i < n && ok; i++) if ((r >> i) & 1) for (int j = i + 1; j < n; j++) if (((r >> j) & 1) && M[i][j] > 2) { ok = false; break; }
        for (int i = 0; i < n && ok; i++) if (((r >> i) & 1) && !fc[rmul[x][i]]) ok = false;
        fc[x] = ok; fcCount += ok;
    }
    fprintf(stderr, "fully commutative elements: %u\n", fcCount);

    // reflections: positive roots via orbit of simple roots; word for t_beta
    {
        map<vector<int>, vector<int>> rootWord; // root -> word of its reflection
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
        fprintf(stderr, "positive roots: %zu\n", reflWords.size());
    }
    // verify reflections have odd length consistent: apply word to identity
    vector<u32> reflElt;
    for (auto& wd : reflWords) { u32 e = 0; for (int s : wd) e = rmul[e][s]; reflElt.push_back(e); }
    { set<u32> distinct(reflElt.begin(), reflElt.end()); if (distinct.size() != reflWords.size()) { fprintf(stderr, "reflection words not distinct\n"); return 1; } }

    // ---- Bruhat order ----
    if (useBitset) {
        size_t words = (N + 63) / 64;
        if ((double)N * words * 8 > 12e9) { fprintf(stderr, "bitset too large, use desc mode\n"); return 1; }
        bruhat.assign(N, vector<uint64_t>(words, 0));
        bruhat[0][0] = 1;
        for (u32 w = 1; w < N; w++) {
            int s = __builtin_ctz(Lmask[w]); u32 v = lmul[w][s];
            bruhat[w] = bruhat[v];
            for (size_t i = 0; i < words; i++) {
                uint64_t b = bruhat[v][i];
                while (b) { int t = __builtin_ctzll(b); b &= b - 1; u32 x = i * 64 + t; u32 sx = lmul[x][s]; bruhat[w][sx >> 6] |= 1ULL << (sx & 63); }
            }
        }
        // sanity: count pairs
        unsigned long long pairs = 0; for (u32 w = 0; w < N; w++) for (auto b : bruhat[w]) pairs += __builtin_popcountll(b);
        fprintf(stderr, "Bruhat comparable pairs (x<=w): %llu\n", pairs);
        // cross-check iterative test on a sample
        mt19937 rng(12345);
        for (int t = 0; t < 200000; t++) { u32 x = rng() % N, w = rng() % N; if (leqBit(x, w) != leqIter(x, w)) { fprintf(stderr, "Bruhat test mismatch\n"); return 1; } }
    }
    // descent buckets for desc mode
    vector<vector<u32>> bucket(1 << (2 * n));
    if (!useBitset) for (u32 x = 0; x < N; x++) bucket[(Lmask[x] << n) | Rmask[x]].push_back(x);

    // ---- KL computation ----
    ex.resize(N); pid.resize(N); mulist.resize(N);
    { vector<ll> one{1}; internPoly(one); }
    ex[0] = {0}; pid[0] = {0};
    // group elements by length (BFS order is nondecreasing length)
    vector<u32> lenStart(maxLen + 2, 0);
    for (u32 w = 0; w < N; w++) lenStart[len[w] + 1] = w + 1;
    for (int l = 1; l <= maxLen + 1; l++) if (lenStart[l] == 0) lenStart[l] = lenStart[l - 1];

    unsigned long long extremalPairs = 1; ll maxCoeff = 1; ll maxMu = 0;
    map<ll, unsigned long long> muHist;
    vector<tuple<u32,u32,ll>> bigMu; // (x, w, mu) with mu >= 2
    mutex statMutex;
    auto t0 = chrono::steady_clock::now();

    for (int l = 1; l <= maxLen; l++) {
        u32 a = lenStart[l], b = lenStart[l + 1];
        atomic<u32> next(a);
        auto worker = [&]() {
            vector<u32> cand;
            while (true) {
                u32 w = next.fetch_add(1); if (w >= b) break;
                int s = __builtin_ctz(Lmask[w]); u32 v = lmul[w][s];
                uint8_t Lw = Lmask[w], Rw = Rmask[w];
                cand.clear();
                if (useBitset) {
                    for (size_t i = 0; i < bruhat[w].size(); i++) {
                        uint64_t bb = bruhat[w][i];
                        while (bb) { int t = __builtin_ctzll(bb); bb &= bb - 1; u32 x = i * 64 + t;
                            if ((Lmask[x] & Lw) == Lw && (Rmask[x] & Rw) == Rw) cand.push_back(x); }
                    }
                } else {
                    // iterate over supersets of Lw and Rw
                    int full = (1 << n) - 1;
                    int freeL = full & ~Lw, freeR = full & ~Rw;
                    for (int subL = freeL;; subL = (subL - 1) & freeL) {
                        int LM = Lw | subL;
                        for (int subR = freeR;; subR = (subR - 1) & freeR) {
                            int RM = Rw | subR;
                            for (u32 x : bucket[(LM << n) | RM]) if (len[x] <= len[w] && leqIter(x, w)) cand.push_back(x);
                            if (subR == 0) break;
                        }
                        if (subL == 0) break;
                    }
                    sort(cand.begin(), cand.end());
                }
                vector<u32> ids; ids.reserve(cand.size());
                vector<pair<u32,ll>> mus;
                vector<ll> P;
                for (u32 x : cand) {
                    if (x == w) { ids.push_back(0); continue; }
                    int d = len[w] - len[x];
                    P.assign(d / 2 + 2, 0);
                    // P_{sx,v}
                    u32 sx = lmul[x][s];
                    const vector<ll>* p1 = lookup(sx, v);
                    if (!p1) { fprintf(stderr, "sx not <= v\n"); abort(); }
                    for (size_t i = 0; i < p1->size(); i++) P[i] += (*p1)[i];
                    // q P_{x,v}
                    const vector<ll>* p2 = lookup(x, v);
                    if (p2) for (size_t i = 0; i < p2->size(); i++) P[i + 1] += (*p2)[i];
                    // subtract mu terms
                    for (auto& zm : mulist[v]) {
                        u32 z = zm.first;
                        if (!((Lmask[z] >> s) & 1)) continue;
                        if (len[z] < len[x]) continue;
                        const vector<ll>* p3 = lookup(x, z);
                        if (!p3) continue;
                        int shift = (len[w] - len[z]) / 2;
                        for (size_t i = 0; i < p3->size(); i++) P[i + shift] -= zm.second * (*p3)[i];
                    }
                    // sanity checks: nonnegative, constant term 1, degree bound
                    if (P[0] != 1) { fprintf(stderr, "constant term != 1 at x=%u w=%u\n", x, w); abort(); }
                    for (size_t i = 0; i < P.size(); i++) if (P[i] < 0) { fprintf(stderr, "negative coefficient at x=%u w=%u\n", x, w); abort(); }
                    for (size_t i = (d - 1) / 2 + 1; i < P.size(); i++) if (P[i] != 0) { fprintf(stderr, "degree bound violated at x=%u w=%u\n", x, w); abort(); }
                    if (d % 2 == 1) { ll mu = P[(d - 1) / 2]; if (mu) mus.push_back({x, mu}); }
                    if (verify) {
                        // recompute with RIGHT recursion w = v' t, t = last right descent of w (t in R(x))
                        int t = 31 - __builtin_clz((unsigned)Rmask[w]); u32 v2 = rmul[w][t];
                        vector<ll> Q(d / 2 + 2, 0);
                        u32 xt = rmul[x][t];
                        const vector<ll>* q1 = lookup(xt, v2); if (!q1) { fprintf(stderr, "xt not <= v2\n"); abort(); }
                        for (size_t i = 0; i < q1->size(); i++) Q[i] += (*q1)[i];
                        const vector<ll>* q2 = lookup(x, v2);
                        if (q2) for (size_t i = 0; i < q2->size(); i++) Q[i + 1] += (*q2)[i];
                        for (auto& zm : mulist[v2]) {
                            u32 z = zm.first;
                            if (!((Rmask[z] >> t) & 1)) continue;
                            if (len[z] < len[x]) continue;
                            const vector<ll>* q3 = lookup(x, z); if (!q3) continue;
                            int shift = (len[w] - len[z]) / 2;
                            for (size_t i = 0; i < q3->size(); i++) Q[i + shift] -= zm.second * (*q3)[i];
                        }
                        for (size_t i = 0; i < P.size(); i++) if (P[i] != Q[i]) { fprintf(stderr, "LEFT/RIGHT recursion mismatch at x=%u w=%u\n", x, w); abort(); }
                    }
                    ids.push_back(internPoly(P));
                }
                // covers that are not extremal: mu = 1
                for (size_t ti = 0; ti < reflElt.size(); ti++) {
                    // z = w t ; compute via word
                    u32 z = w; for (int si : reflWords[ti]) z = rmul[z][si];
                    if (len[z] != len[w] - 1) continue;
                    if ((Lmask[z] & Lw) == Lw && (Rmask[z] & Rw) == Rw) continue; // extremal, already included
                    mus.push_back({z, 1});
                }
                sort(mus.begin(), mus.end());
                ex[w] = cand; pid[w] = ids; mulist[w] = mus;
                {
                    lock_guard<mutex> g(statMutex);
                    extremalPairs += cand.size();
                    for (auto& zm : mus) { muHist[zm.second]++; if (zm.second > maxMu) maxMu = zm.second; if (zm.second >= 2) bigMu.emplace_back(zm.first, w, zm.second); }
                }
            }
        };
        vector<thread> pool; int T = max(1, threads);
        for (int t = 0; t < T; t++) pool.emplace_back(worker);
        for (auto& th : pool) th.join();
        if (l % 4 == 0 || l == maxLen) {
            double secs = chrono::duration<double>(chrono::steady_clock::now() - t0).count();
            fprintf(stderr, "  length %d done (%u elements), extremal pairs so far %llu, distinct polys %u, %.1fs\n", l, b - a, extremalPairs, polyCount, secs);
        }
    }
    for (u32 i = 0; i < polyCount; i++) for (ll c : getPoly(i)) maxCoeff = max(maxCoeff, c);
    // covers have mu=1 and were counted; add count of cover pairs separately
    unsigned long long coverPairs = 0, nonCoverMu1 = 0;
    for (u32 w = 0; w < N; w++) for (auto& zm : mulist[w]) { if (len[w] - len[zm.first] == 1) coverPairs++; else if (zm.second == 1) nonCoverMu1++; }

    printf("RESULT type=%c%d |W|=%u longest=%d FC=%u\n", type, n, N, maxLen, fcCount);
    printf("extremal_pairs=%llu distinct_KL_polynomials=%u max_coefficient=%lld\n", extremalPairs, polyCount, maxCoeff);
    printf("max_mu=%lld\n", maxMu);
    printf("nonzero mu pairs: covers=%llu, noncover mu=1: %llu", coverPairs, nonCoverMu1);
    for (auto& kv : muHist) if (kv.first >= 2) printf(", mu=%lld: %llu", kv.first, kv.second);
    printf("\n");
    // check: P_{e,w0} and longest element
    sort(bigMu.begin(), bigMu.end(), [&](auto& p, auto& q) { if (get<2>(p) != get<2>(q)) return get<2>(p) > get<2>(q); return len[get<1>(p)] < len[get<1>(q)]; });
    size_t show = min<size_t>(bigMu.size(), 40);
    for (size_t i = 0; i < show; i++) {
        u32 x = get<0>(bigMu[i]), w = get<1>(bigMu[i]); ll mu = get<2>(bigMu[i]);
        const vector<ll>* p = lookup(x, w);
        string ps; for (size_t k = 0; k < p->size(); k++) { if ((*p)[k] == 0) continue; if (!ps.empty()) ps += " + "; ps += to_string((*p)[k]); if (k) ps += "q^" + to_string(k); }
        printf("mu=%lld  x=%s (len %d, FC=%d, L=%s R=%s)  w=%s (len %d, FC=%d, L=%s R=%s)  P_{x,w}=%s\n",
               mu, wordStr(reducedWord(x)).c_str(), len[x], (int)fc[x], maskStr(Lmask[x]).c_str(), maskStr(Rmask[x]).c_str(),
               wordStr(reducedWord(w)).c_str(), len[w], (int)fc[w], maskStr(Lmask[w]).c_str(), maskStr(Rmask[w]).c_str(), ps.c_str());
    }
    if (bigMu.size() > show) printf("... (%zu pairs with mu>=2 in total)\n", bigMu.size());
    // FC-bottom statistics: among pairs with mu != 0 and x FC, the mu values
    {
        ll maxMuFC = 0; unsigned long long cntFC = 0;
        for (u32 w = 0; w < N; w++) for (auto& zm : mulist[w]) if (fc[zm.first]) { cntFC++; maxMuFC = max(maxMuFC, zm.second); }
        printf("pairs with FC lower endpoint and mu!=0: %llu, max mu among them: %lld\n", cntFC, maxMuFC);
    }
    if (!dumpFile.empty()) {
        FILE* f = fopen(dumpFile.c_str(), "w");
        for (u32 w = 0; w < N; w++) {
            string ww = wordStr(reducedWord(w));
            for (u32 x = 0; x < N; x++) {
                if (!leq(x, w)) continue;
                const vector<ll>* p = lookup(x, w);
                fprintf(f, "%s %s", wordStr(reducedWord(x)).c_str(), ww.c_str());
                for (ll c : *p) fprintf(f, " %lld", c);
                fprintf(f, "\n");
            }
        }
        fclose(f);
    }
    return 0;
}
