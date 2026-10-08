// Independent enumeration of fully commutative elements of E_n (chain 0..n-2, node n-1 at 2).
// Height-vector representation p_i = ht(w(alpha_i)); FC test via parabolic strip:
// for FC w and s not a right descent, ws is non-FC iff w^J (J = nodes commuting with s,
// J excludes s) has a reduced word ending in s t with t adjacent to s.
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <algorithm>
#include <cstring>
#include <cstdint>
typedef std::array<int32_t,16> Vec;
int n; int A[16][16]; std::vector<int> nb[16], J[16];
inline void mul(Vec& p, int s){ int32_t ps=p[s]; for(int i=0;i<n;i++) if(A[i][s]) p[i]-=A[i][s]*ps; }
bool extend_fc(const Vec& p, int s){
    Vec q=p; bool changed=true;
    while(changed){ changed=false; for(int j: J[s]) if(q[j]<0){ mul(q,j); changed=true; } }
    for(int t: nb[s]) if(q[t]<0){ Vec qt=q; mul(qt,t); if(qt[s]<0) return false; }
    return true;
}
int main(int argc,char**argv){
    n=atoi(argv[1]);
    memset(A,0,sizeof A);
    for(int i=0;i<n;i++) A[i][i]=2;
    for(int i=0;i+1<=n-2;i++){ A[i][i+1]=A[i+1][i]=-1; }
    A[2][n-1]=A[n-1][2]=-1;
    for(int s=0;s<n;s++) for(int j=0;j<n;j++){ if(j==s) continue; if(A[s][j]) nb[s].push_back(j); else J[s].push_back(j); }
    Vec id; id.fill(0); for(int i=0;i<n;i++) id[i]=1;
    std::vector<Vec> layer{id}; long long total=0; int len=0; Vec lastmax=id;
    while(!layer.empty()){
        total+=layer.size();
        std::vector<Vec> nxt; nxt.reserve(layer.size()*2+16);
        for(const Vec& p: layer) for(int s=0;s<n;s++) if(p[s]>0 && extend_fc(p,s)){ Vec q=p; mul(q,s); nxt.push_back(q); }
        std::sort(nxt.begin(),nxt.end()); nxt.erase(std::unique(nxt.begin(),nxt.end()),nxt.end());
        if(!nxt.empty()){ len++; lastmax=nxt[0]; }
        else break;
        layer.swap(nxt);
        if(len%10==0){ fprintf(stderr,"E%d length %d layer %zu total %lld\n",n,len,layer.size(),total); }
    }
    printf("E%d FC count %lld max length %d\n",n,total,len);
    // a reduced word for one maximal element, by stripping right descents
    Vec p=lastmax; std::vector<int> word; std::vector<int> mult(n,0);
    for(;;){ int s=-1; for(int i=0;i<n;i++) if(p[i]<0){ s=i; break; } if(s<0) break; mul(p,s); word.push_back(s); mult[s]++; }
    std::reverse(word.begin(),word.end());
    printf("  max element word:"); for(int s: word) printf(" %d",s); printf("\n  multiplicities:"); for(int i=0;i<n;i++) printf(" %d",mult[i]); printf("\n");
    return 0;
}
