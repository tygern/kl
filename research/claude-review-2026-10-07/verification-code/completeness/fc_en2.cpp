// FC catalogue of generalized E_n (paper numbering) by level-wise closure under right
// ascents with the recurrence  w in FC <=> R(w) commuting and w s in FC for all s in R(w).
// Only two consecutive levels are kept.  Exact int32 arithmetic in the geometric representation.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <string>
#include <unordered_set>
#include <cstdint>
using namespace std;
static int n;
static vector<vector<int>> adj;
typedef vector<int32_t> Elt; // n*n, column j at [j*n .. j*n+n)
struct H { size_t operator()(const Elt& e) const { size_t h=1469598103934665603ULL; for (auto v: e){ h^=(size_t)(uint32_t)v; h*=1099511628211ULL;} return h; } };
static inline void rmul(Elt& w, int i) { // w <- w s_i
    for (int j: adj[i]) for (int k=0;k<n;k++) w[j*n+k]+=w[i*n+k];
    for (int k=0;k<n;k++) w[i*n+k]=-w[i*n+k];
}
static inline bool neg(const Elt& w, int j){ bool any=false; for(int k=0;k<n;k++){ if(w[j*n+k]>0) return false; if(w[j*n+k]<0) any=true;} return any; }
static int lengthOf(Elt w){ int l=0; while(true){ int d=-1; for(int j=0;j<n;j++) if(neg(w,j)){d=j;break;} if(d<0) return l; rmul(w,d); l++; } }
static Elt inverseOf(Elt w){ vector<int> word; while(true){ int d=-1; for(int j=0;j<n;j++) if(neg(w,j)){d=j;break;} if(d<0) break; rmul(w,d); word.push_back(d);} Elt v(n*n,0); for(int j=0;j<n;j++) v[j*n+j]=1; for(int i=(int)word.size()-1;i>=0;i--) rmul(v,word[i]); /* w = s_{word[k-1]}...s_{word[0]} ; w^{-1} = s_{word[0]}...s_{word[k-1]} */ Elt u(n*n,0); for(int j=0;j<n;j++) u[j*n+j]=1; for(int i=0;i<(int)word.size();i++) rmul(u,word[i]); return u; }
int main(int argc,char**argv){
    n=atoi(argv[1]);
    vector<int> I; if(argc>2){ char* p=strtok(argv[2],","); while(p){ I.push_back(atoi(p)); p=strtok(NULL,","); } }
    long long hits=0; int hitlen=-1; adj.assign(n,{});
    for(int i=0;i<n-2;i++){adj[i].push_back(i+1);adj[i+1].push_back(i);} adj[n-1].push_back(2);adj[2].push_back(n-1);
    vector<vector<char>> A(n,vector<char>(n,0)); for(int i=0;i<n;i++)for(int j:adj[i])A[i][j]=1;
    Elt e(n*n,0); for(int j=0;j<n;j++) e[j*n+j]=1;
    unordered_set<Elt,H> prev; prev.insert(e);
    long long total=1; int len=0;
    while(true){
        unordered_set<Elt,H> nxt;
        for(const Elt& u: prev){
            for(int s=0;s<n;s++){
                if(neg(u,s)) continue;
                Elt w=u; rmul(w,s);
                if(nxt.count(w)) continue;
                vector<int> R; for(int j=0;j<n;j++) if(neg(w,j)) R.push_back(j);
                bool ok=true;
                for(size_t a=0;a<R.size()&&ok;a++) for(size_t b=a+1;b<R.size();b++) if(A[R[a]][R[b]]){ok=false;break;}
                if(!ok) continue;
                for(int t: R){ Elt v=w; rmul(v,t); if(!prev.count(v)){ok=false;break;} }
                if(ok){ if(!I.empty()){ bool sub=true; for(int i: I) if(!neg(w,i)){sub=false;break;} if(sub){ Elt wi=inverseOf(w); bool subL=true; for(int i: I) if(!neg(wi,i)){subL=false;break;} if(subL){ hits++; hitlen=lengthOf(w); fprintf(stderr,"  hit: length %d\n",hitlen);} } } nxt.insert(std::move(w)); }
            }
        }
        if(nxt.empty()) break;
        len++; total+=nxt.size();
        fprintf(stderr,"E%d level %d: %zu\n",n,len,nxt.size());
        prev.swap(nxt);
    }
    printf("E%d FC count %lld max length %d\n",n,total,len); if(!I.empty()) printf("  FC elements with I in L and R: %lld (last hit length %d)\n",hits,hitlen);
}
