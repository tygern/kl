// Fresh FC enumerator for simply-laced Coxeter graphs.
// Each FC element is visited exactly once via its lexicographic (Anisimov-Knuth)
// normal form; FC-ness of w*s (w FC, s not a right descent) is decided by
// Stembridge's criterion: the number of N(s)-letters after the last s must be >= 2.
// Usage: fc_enum E n   |  fc_enum D n   |  fc_enum A n
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
using namespace std;

int n;
bool adj[64][64];
vector<int> nb[64];
long long count_ = 0;
int maxlen = 0;
int word[4096];
int cnt[64];       // number of N(s) letters since last s; -1 if s never occurred
int bestword[4096];
int bestlen = 0;
int ndistinct_best = 0;

void dfs(int len) {
    count_++;
    if (len > maxlen) { maxlen = len; bestlen = len; memcpy(bestword, word, len*sizeof(int)); }
    for (int s = 0; s < n; s++) {
        // FC + length-increasing criterion
        if (cnt[s] >= 0 && cnt[s] < 2) continue;
        // lexicographic normal form check: scanning back while letters commute with s,
        // none may be > s
        bool ok = true;
        for (int j = len - 1; j >= 0; j--) {
            int t = word[j];
            if (t == s || adj[s][t]) break;
            if (t > s) { ok = false; break; }
        }
        if (!ok) continue;
        // apply
        int saved[64]; int ns = (int)nb[s].size();
        for (int i = 0; i < ns; i++) { int t = nb[s][i]; saved[i] = cnt[t]; if (cnt[t] >= 0) cnt[t]++; }
        int savedS = cnt[s]; cnt[s] = 0;
        word[len] = s;
        dfs(len + 1);
        cnt[s] = savedS;
        for (int i = 0; i < ns; i++) cnt[nb[s][i]] = saved[i];
    }
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: fc_enum {A|D|E} n\n"); return 1; }
    char type = argv[1][0]; n = atoi(argv[2]);
    memset(adj, 0, sizeof adj);
    auto edge = [&](int a, int b){ adj[a][b] = adj[b][a] = true; };
    if (type == 'E') {           // chain 0-1-...-(n-2), node n-1 attached to 2
        for (int i = 0; i + 1 <= n - 2; i++) edge(i, i + 1);
        edge(n - 1, 2);
    } else if (type == 'D') {    // chain 0-...-(n-2), node n-1 attached to n-3
        for (int i = 0; i + 1 <= n - 2; i++) edge(i, i + 1);
        edge(n - 1, n - 3);
    } else {                     // A_n chain
        for (int i = 0; i + 1 <= n - 1; i++) edge(i, i + 1);
    }
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (adj[i][j]) nb[i].push_back(j);
    for (int i = 0; i < n; i++) cnt[i] = -1;
    dfs(0);
    printf("%c%d: FC count = %lld, max FC length = %d\n", type, n, count_, maxlen);
    printf("one max-length element (lex normal form): ");
    for (int i = 0; i < bestlen; i++) printf("%d ", bestword[i]);
    printf("\nmultiplicities:");
    int mult[64]; memset(mult, 0, sizeof mult);
    for (int i = 0; i < bestlen; i++) mult[bestword[i]]++;
    for (int i = 0; i < n; i++) printf(" %d", mult[i]);
    printf("\n");
    return 0;
}
