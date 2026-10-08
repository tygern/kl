// Complete enumeration of the fully commutative elements of the generalized
// Coxeter group E_n (chain 0-1-...-(n-2), node n-1 attached to node 2), with
// the count and the maximum length.  Written for the manuscript's remark that
// a length-based exclusion of fully commutative Bruhat covers would not work
// for the rank-uniform family: the maximum FC lengths of E_10, ..., E_13 are
// 55, 66, 78, 92, while l(b_{3,0}) = 75 in E_13.
//
// Representation.  An element w is stored by the integer vector
// h(w)_i = height of w(alpha_i) (sum of the simple-root coordinates of the
// column w(alpha_i)).  Since h(w)_i = rho(w alpha_i) for the linear functional
// rho with rho(alpha_k) = 1 for all k, and rho lies in the interior of the
// fundamental chamber of the contragredient action, on which W acts freely,
// h(w) determines w.  Right multiplication by s_j acts on h as on the columns:
// h_i += h_j for i adjacent to j, then h_j = -h_j; s_j is a right descent iff
// h_j < 0.  Full commutativity is decided by the recurrence (Stembridge 1996,
// Proposition 2.1 with the heap description)
//   w is FC  <=>  R(w) is pairwise commuting and w t is FC for every t in R(w),
// applied level by level: the FC elements of length l+1 are the elements w s
// with w FC of length l, s not a right descent of w, such that R(ws) is
// commuting and every (ws)t with t in R(ws) is among the FC elements of length
// l.  Only two consecutive levels are kept in memory.  Independent of the
// repository's FC catalogue program, which uses full integer matrices.
//
// Build: c++ -std=c++17 -O3 fc_maxima.cpp -o fc_maxima
// Run:   ./fc_maxima <output.json> <rank> [rank ...]   (ranks between 6 and 15)
// Known values asserted: E6 662/16, E7 2670/27, E8 10846/29 (Stembridge 1998,
// Table 1, counts; maximum lengths from the repository catalogues), E9 44199/44
// (Biagioli-Jouhet-Nadeau, Section 5.2), E10 180438/55, E11 737762/66,
// E12 3021000/78, E13 12387990/92 (review computation of 7 October 2026).
// With several threads E13 takes a few seconds; E14 and E15 are much larger.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>
using namespace std;

static int n;
static int T = 8;
static vector<vector<int>> adj;
static bool A[20][20];

struct Key {
    int16_t h[16];
    bool operator==(const Key& o) const { return memcmp(h, o.h, sizeof(h)) == 0; }
};
struct KH {
    size_t operator()(const Key& k) const {
        uint64_t x = 1469598103934665603ULL;
        for (int i = 0; i < 16; i++) { x ^= (uint16_t)k.h[i]; x *= 1099511628211ULL; x ^= x >> 29; }
        return x;
    }
};
static inline uint64_t shardOf(const Key& k) { return KH()(k) * 0x9E3779B97F4A7C15ULL >> 40; }
static inline void rmul(Key& w, int j) {
    for (int i : adj[j]) { int v = w.h[i] + w.h[j]; if (v > 32000 || v < -32000) { fprintf(stderr, "height overflow\n"); exit(3); } w.h[i] = (int16_t)v; }
    w.h[j] = (int16_t)(-w.h[j]);
}

struct Outcome { long long count; int maxlen; vector<long long> levels; vector<string> longest; };

static Outcome enumerate(int rank) {
    n = rank;
    adj.assign(n, {});
    for (int i = 0; i < n - 2; i++) { adj[i].push_back(i + 1); adj[i + 1].push_back(i); }
    adj[n - 1].push_back(2); adj[2].push_back(n - 1);
    memset(A, 0, sizeof(A));
    for (int i = 0; i < n; i++) for (int j : adj[i]) A[i][j] = true;
    Key e; memset(e.h, 0, sizeof(e.h)); for (int i = 0; i < n; i++) e.h[i] = 1;
    vector<unordered_set<Key, KH>> prev(T), nxt(T);
    prev[shardOf(e) % T].insert(e);
    Outcome out; out.count = 1; out.maxlen = 0; out.levels = {1};
    while (true) {
        vector<vector<vector<Key>>> buckets(T, vector<vector<Key>>(T));
        vector<thread> th;
        for (int t = 0; t < T; t++) th.emplace_back([&, t]() {
            for (const Key& u : prev[t]) {
                for (int s = 0; s < n; s++) {
                    if (u.h[s] < 0) continue;  // s in R(u)
                    Key v = u; rmul(v, s);
                    int R[20], nr = 0;
                    for (int j = 0; j < n; j++) if (v.h[j] < 0) R[nr++] = j;
                    bool ok = true;
                    for (int a = 0; a < nr && ok; a++) for (int b = a + 1; b < nr; b++) if (A[R[a]][R[b]]) { ok = false; break; }
                    if (!ok) continue;
                    for (int a = 0; a < nr; a++) {
                        if (R[a] == s) continue;  // v s = u is FC by construction
                        Key x = v; rmul(x, R[a]);
                        if (!prev[shardOf(x) % T].count(x)) { ok = false; break; }
                    }
                    if (!ok) continue;
                    buckets[t][shardOf(v) % T].push_back(v);
                }
            }
        });
        for (auto& x : th) x.join();
        th.clear();
        for (int t = 0; t < T; t++) th.emplace_back([&, t]() {
            nxt[t].clear();
            size_t tot = 0; for (int s = 0; s < T; s++) tot += buckets[s][t].size();
            nxt[t].reserve(tot);
            for (int s = 0; s < T; s++) { for (const Key& k : buckets[s][t]) nxt[t].insert(k); vector<Key>().swap(buckets[s][t]); }
        });
        for (auto& x : th) x.join();
        size_t cnt = 0; for (int t = 0; t < T; t++) cnt += nxt[t].size();
        if (cnt == 0) break;
        out.maxlen++; out.count += cnt; out.levels.push_back((long long)cnt);
        prev.swap(nxt);
    }
    // Reduced words of the longest elements (strip the smallest right descent repeatedly).
    for (int t = 0; t < T; t++) for (const Key& k : prev[t]) {
        Key w = k; vector<int> word;
        while (true) { int d = -1; for (int j = 0; j < n; j++) if (w.h[j] < 0) { d = j; break; } if (d < 0) break; rmul(w, d); word.push_back(d); }
        reverse(word.begin(), word.end());
        string s; for (size_t i = 0; i < word.size(); i++) { if (i) s += " "; s += to_string(word[i]); }
        out.longest.push_back(s);
    }
    sort(out.longest.begin(), out.longest.end());
    return out;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: fc_maxima <output.json> <rank> [rank ...]\n"); return 2; }
    T = (int)min<unsigned>(8, max<unsigned>(1, thread::hardware_concurrency()));
    map<int, pair<long long, int>> known = {{6, {662, 16}}, {7, {2670, 27}}, {8, {10846, 29}}, {9, {44199, 44}},
                                            {10, {180438, 55}}, {11, {737762, 66}}, {12, {3021000, 78}}, {13, {12387990, 92}}};
    int failures = 0;
    FILE* f = fopen(argv[1], "w"); if (!f) { fprintf(stderr, "cannot open output\n"); return 2; }
    fprintf(f, "{\n  \"labelling\": \"E_n: chain 0-1-...-(n-2), node n-1 attached to node 2\",\n");
    fprintf(f, "  \"method\": \"level-wise closure under right ascents with the recurrence: w FC iff R(w) commuting and wt FC for every t in R(w); elements stored as height vectors h(w)_i = ht(w(alpha_i))\",\n");
    fprintf(f, "  \"ranks\": {\n");
    for (int a = 2; a < argc; a++) {
        int rank = atoi(argv[a]);
        if (rank < 6 || rank > 15) { fprintf(stderr, "rank out of range\n"); return 2; }
        Outcome o = enumerate(rank);
        bool ok = true;
        if (known.count(rank)) { ok = (o.count == known[rank].first && o.maxlen == known[rank].second); if (!ok) { fprintf(stderr, "EXPECTATION FAILED for E%d: %lld / %d\n", rank, o.count, o.maxlen); failures++; } }
        fprintf(f, "    \"E%d\": {\"fully_commutative_count\": %lld, \"maximum_length\": %d, \"longest_elements\": %zu,\n", rank, o.count, o.maxlen, o.longest.size());
        fprintf(f, "      \"length_distribution\": ["); for (size_t i = 0; i < o.levels.size(); i++) fprintf(f, "%s%lld", i ? "," : "", o.levels[i]); fprintf(f, "],\n");
        fprintf(f, "      \"longest_reduced_words\": ["); for (size_t i = 0; i < o.longest.size(); i++) fprintf(f, "%s\"%s\"", i ? ", " : "", o.longest[i].c_str()); fprintf(f, "],\n");
        fprintf(f, "      \"matches_expected\": %s}%s\n", known.count(rank) ? (ok ? "true" : "false") : "null", a + 1 < argc ? "," : "");
        printf("E%d: %lld fully commutative elements, maximum length %d%s\n", rank, o.count, o.maxlen, known.count(rank) ? (ok ? " (as expected)" : " (UNEXPECTED)") : "");
        fflush(stdout);
    }
    fprintf(f, "  },\n  \"status\": \"%s\"\n}\n", failures ? "FAILED" : "passed");
    fclose(f);
    return failures ? 1 : 0;
}
