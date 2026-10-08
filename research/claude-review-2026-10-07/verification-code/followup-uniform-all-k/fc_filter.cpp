// Same level-wise FC closure as fcmax.cpp (height-vector representation), restricted to
// reporting: FC elements x of E_n with |R(x) ∩ I| >= |I|-1 and |L(x) ∩ I| >= |I|-1 and
// (R(x) ⊇ I or L(x) ⊇ I), for a given independent set I.  Prints all of them with
// L, R, length and a reduced word.  Usage: fc_filter n T i1,i2,...
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <thread>
#include <mutex>
#include <unordered_set>
#include <algorithm>
using namespace std;
static int n; static int T = 8;
static vector<vector<int>> adj; static bool A[20][20];
struct Key { int16_t h[16]; bool operator==(const Key& o) const { return memcmp(h, o.h, sizeof(h)) == 0; } };
struct KH { size_t operator()(const Key& k) const { uint64_t x = 1469598103934665603ULL; for (int i = 0; i < 16; i++) { x ^= (uint16_t)k.h[i]; x *= 1099511628211ULL; x ^= x >> 29; } return x; } };
static inline uint64_t shardOf(const Key& k) { return KH()(k) * 0x9E3779B97F4A7C15ULL >> 40; }
static inline void rmul(Key& w, int j) { for (int i : adj[j]) w.h[i] += w.h[j]; w.h[j] = -w.h[j]; }
static vector<int> wordOf(Key w) { vector<int> word; while (true) { int d = -1; for (int j = 0; j < n; j++) if (w.h[j] < 0) { d = j; break; } if (d < 0) break; rmul(w, d); word.push_back(d); } reverse(word.begin(), word.end()); return word; }
static mutex outm;
int main(int argc, char** argv) {
    n = atoi(argv[1]); T = atoi(argv[2]);
    vector<int> I; { char* p = strtok(argv[3], ","); while (p) { I.push_back(atoi(p)); p = strtok(NULL, ","); } }
    bool inI[20] = {false}; for (int i : I) inI[i] = true;
    adj.assign(n, {});
    for (int i = 0; i < n - 2; i++) { adj[i].push_back(i + 1); adj[i + 1].push_back(i); }
    adj[n - 1].push_back(2); adj[2].push_back(n - 1);
    memset(A, 0, sizeof(A)); for (int i = 0; i < n; i++) for (int j : adj[i]) A[i][j] = true;
    Key e; memset(e.h, 0, sizeof(e.h)); for (int i = 0; i < n; i++) e.h[i] = 1;
    vector<unordered_set<Key, KH>> prev(T), nxt(T);
    prev[shardOf(e) % T].insert(e);
    long long total = 1; int len = 0; long long reported = 0;
    while (true) {
        vector<vector<vector<Key>>> buckets(T, vector<vector<Key>>(T));
        vector<thread> th;
        for (int t = 0; t < T; t++) th.emplace_back([&, t]() {
            for (const Key& u : prev[t]) for (int s = 0; s < n; s++) {
                if (u.h[s] < 0) continue;
                Key v = u; rmul(v, s);
                int R[20], nr = 0; for (int j = 0; j < n; j++) if (v.h[j] < 0) R[nr++] = j;
                bool ok = true;
                for (int a = 0; a < nr && ok; a++) for (int b = a + 1; b < nr; b++) if (A[R[a]][R[b]]) { ok = false; break; }
                if (!ok) continue;
                for (int a = 0; a < nr; a++) { if (R[a] == s) continue; Key x = v; rmul(x, R[a]); if (!prev[shardOf(x) % T].count(x)) { ok = false; break; } }
                if (!ok) continue;
                buckets[t][shardOf(v) % T].push_back(v);
            }
        });
        for (auto& x : th) x.join(); th.clear();
        for (int t = 0; t < T; t++) th.emplace_back([&, t]() {
            nxt[t].clear(); size_t tot = 0; for (int s = 0; s < T; s++) tot += buckets[s][t].size(); nxt[t].reserve(tot);
            for (int s = 0; s < T; s++) { for (const Key& k : buckets[s][t]) nxt[t].insert(k); vector<Key>().swap(buckets[s][t]); }
            // filter
            for (const Key& k : nxt[t]) {
                int rc = 0; for (int i : I) if (k.h[i] < 0) rc++;
                if (rc < (int)I.size() - 1) continue;
                vector<int> word = wordOf(k);
                Key inv = e; for (int i = (int)word.size() - 1; i >= 0; i--) rmul(inv, word[i]);
                int lc = 0; for (int i : I) if (inv.h[i] < 0) lc++;
                if (lc < (int)I.size() - 1) continue;
                if (rc < (int)I.size() && lc < (int)I.size()) continue;
                lock_guard<mutex> g(outm);
                reported++;
                printf("len %zu  L={", word.size()); for (int j = 0; j < n; j++) if (inv.h[j] < 0) printf("%d,", j);
                printf("} R={"); for (int j = 0; j < n; j++) if (k.h[j] < 0) printf("%d,", j);
                printf("} word="); for (int x : word) printf("%d ", x); printf("\n");
            }
        });
        for (auto& x : th) x.join();
        size_t cnt = 0; for (int t = 0; t < T; t++) cnt += nxt[t].size();
        if (cnt == 0) break;
        len++; total += cnt; prev.swap(nxt);
    }
    printf("E%d FC count %lld max length %d; reported %lld elements\n", n, total, len, reported);
}
