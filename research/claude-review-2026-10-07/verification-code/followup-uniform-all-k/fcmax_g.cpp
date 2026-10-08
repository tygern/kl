// Independent FC catalogue of generalized E_n (chain 0-1-...-(n-2), node n-1 attached to 2).
// Elements are represented by the integer vector h(w)_i = height of w(alpha_i) (sum of the
// simple-root coordinates of column i).  Since h(w)_i = rho(w alpha_i) for the functional rho
// with rho(alpha_k)=1 for all k, and the contragredient W-action is free on the interior of the
// fundamental chamber, h(w) determines w.  Right multiplication by s_j acts on h exactly as on
// the columns: h_i += h_j for i ~ j, then h_j = -h_j.  s_j in R(w) iff h_j < 0.
// Level-wise closure under right ascents with the recurrence
//   v in FC  <=>  R(v) pairwise commuting and v t in FC for every t in R(v),
// keeping two levels, sharded hash sets, multithreaded.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <thread>
#include <unordered_set>
#include <algorithm>
#include <string>
using namespace std;

static int n;
static int T = 8;                 // threads / shards
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
    for (int i : adj[j]) w.h[i] += w.h[j];
    w.h[j] = -w.h[j];
}

int main(int argc, char** argv) {
    n = atoi(argv[1]);
    if (argc > 2) T = atoi(argv[2]);
    adj.assign(n, {});
    for (int i = 0; i < n - 2; i++) { adj[i].push_back(i + 1); adj[i + 1].push_back(i); }
    int br = (argc > 3 && argv[3][0]==(char)68) ? n - 3 : 2; adj[n - 1].push_back(br); adj[br].push_back(n - 1);
    memset(A, 0, sizeof(A));
    for (int i = 0; i < n; i++) for (int j : adj[i]) A[i][j] = true;

    Key e; memset(e.h, 0, sizeof(e.h)); for (int i = 0; i < n; i++) e.h[i] = 1;
    vector<unordered_set<Key, KH>> prev(T), nxt(T);
    prev[shardOf(e) % T].insert(e);
    long long total = 1; int len = 0;
    vector<long long> levelCounts{1};
    while (true) {
        // phase 1: each thread processes its shard of prev, producing candidate buckets per target shard
        vector<vector<vector<Key>>> buckets(T, vector<vector<Key>>(T));
        vector<thread> th;
        for (int t = 0; t < T; t++) th.emplace_back([&, t]() {
            for (const Key& u : prev[t]) {
                for (int s = 0; s < n; s++) {
                    if (u.h[s] < 0) continue;      // s in R(u): ws shorter
                    Key v = u; rmul(v, s);
                    int R[20], nr = 0;
                    for (int j = 0; j < n; j++) if (v.h[j] < 0) R[nr++] = j;
                    bool ok = true;
                    for (int a = 0; a < nr && ok; a++) for (int b = a + 1; b < nr; b++) if (A[R[a]][R[b]]) { ok = false; break; }
                    if (!ok) continue;
                    for (int a = 0; a < nr; a++) {
                        if (R[a] == s) continue;   // v s = u is in prev by construction
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
        // phase 2: thread t merges all candidate buckets of shard t
        for (int t = 0; t < T; t++) th.emplace_back([&, t]() {
            nxt[t].clear();
            size_t tot = 0; for (int s = 0; s < T; s++) tot += buckets[s][t].size();
            nxt[t].reserve(tot);
            for (int s = 0; s < T; s++) { for (const Key& k : buckets[s][t]) nxt[t].insert(k); vector<Key>().swap(buckets[s][t]); }
        });
        for (auto& x : th) x.join();
        size_t cnt = 0; for (int t = 0; t < T; t++) cnt += nxt[t].size();
        if (cnt == 0) break;
        len++; total += cnt; levelCounts.push_back((long long)cnt);
        fprintf(stderr, "E%d level %d: %zu\n", n, len, cnt);
        prev.swap(nxt);
    }
    printf("E%d FC count %lld max length %d\n", n, total, len);
    // print the top-level elements as reduced words (descent stripping) with multiplicity profiles
    int shown = 0;
    for (int t = 0; t < T; t++) for (const Key& k : prev[t]) {
        Key w = k; vector<int> word;
        while (true) { int d = -1; for (int j = 0; j < n; j++) if (w.h[j] < 0) { d = j; break; } if (d < 0) break; rmul(w, d); word.push_back(d); }
        reverse(word.begin(), word.end());
        vector<int> mult(n, 0); for (int x : word) mult[x]++;
        printf("  maxlen element %d: word=", shown);
        for (int x : word) printf("%d%s", x, x >= 10 ? "." : "");
        printf("  mult=");
        for (int i = 0; i < n; i++) printf("%d%s", mult[i], i + 1 < n ? "," : "");
        printf("\n");
        if (++shown >= 40) break;
    }
    printf("  levels:");
    for (size_t i = 0; i < levelCounts.size(); i++) printf(" %lld", levelCounts[i]);
    printf("\n");
    return 0;
}
