// klcheck.cpp -- verifier's own full-group KL / mu computation (exact integers).
// Numbering: chain 0-1-...-(n-2); D: node n-1 attached to node n-3; E: node n-1 attached to node 2.
// Usage: klcheck <A|D|E> <rank> [xword wword]...
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <map>
using namespace std;

static int n;
static int Cart[8][8];
struct Mat { int8_t m[8][8]; bool operator==(const Mat&o) const { return memcmp(m,o.m,64)==0; } };
struct MH { size_t operator()(const Mat&a) const { uint64_t h=14695981039346656037ULL; for(int i=0;i<64;i++){h^=((const unsigned char*)a.m)[i];h*=1099511628211ULL;} return h; } };

static Mat rmul(const Mat&w,int i){ // w * s_i : row j <- row j - Cart[i][j] row i ; row i <- -row i
    Mat r=w;
    for(int j=0;j<n;j++) if(j!=i && Cart[i][j]) for(int k=0;k<n;k++) r.m[j][k]=(int8_t)(w.m[j][k]-Cart[i][j]*w.m[i][k]);
    for(int k=0;k<n;k++) r.m[i][k]=(int8_t)(-w.m[i][k]);
    return r;
}
static Mat lmul(int i,const Mat&w){ // s_i * w : row j <- row j - <alpha_i^vee, row j> alpha_i
    Mat r=w;
    for(int j=0;j<n;j++){ int p=0; for(int k=0;k<n;k++) p+=Cart[i][k]*w.m[j][k]; r.m[j][i]=(int8_t)(w.m[j][i]-p); }
    return r;
}
static bool rowneg(const Mat&w,int i){ for(int k=0;k<n;k++) if(w.m[i][k]) return w.m[i][k]<0; return false; }

static int N;
static vector<Mat> els; static unordered_map<Mat,int,MH> idx;
static vector<int> len, par, parlet; static vector<int> R[8], L[8]; // R[s][w]=index of w s ; L[s][w]=index of s w
static vector<unsigned> Rdes, Ldes;
static size_t WPL; static vector<uint64_t> low; // bitsets
static inline bool leq(int x,int w){ return (low[(size_t)w*WPL+(x>>6)]>>(x&63))&1ULL; }
static inline void setb(int w,int x){ low[(size_t)w*WPL+(x>>6)]|=1ULL<<(x&63); }

#define MAXD 20
struct Poly { int32_t c[MAXD]; };
static vector<vector<int>> ext; static vector<vector<Poly>> pol;
static vector<vector<pair<int,int>>> muL;
static Poly ZERO;

static const Poly* lookup(int x,int w){
    if(!leq(x,w)) return &ZERO;
    for(;;){
        bool mv=false;
        for(int s=0;s<n&&!mv;s++) if((Ldes[w]>>s)&1){ int y=L[s][x]; if(len[y]>len[x]){x=y;mv=true;} }
        if(mv) continue;
        for(int s=0;s<n&&!mv;s++) if((Rdes[w]>>s)&1){ int y=R[s][x]; if(len[y]>len[x]){x=y;mv=true;} }
        if(!mv) break;
    }
    auto&E=ext[w]; auto it=lower_bound(E.begin(),E.end(),x);
    if(it==E.end()||*it!=x){ fprintf(stderr,"lookup failure x=%d w=%d\n",x,w); exit(1); }
    return &pol[w][it-E.begin()];
}
static int mu(int z,int v){ int d=len[v]-len[z]; if(d<=0||d%2==0) return 0; return lookup(z,v)->c[(d-1)/2]; }
static string word(int w){ string s; while(w){ s.push_back('0'+parlet[w]); w=par[w]; } reverse(s.begin(),s.end()); return s; }
static int elem(const string&s){ int w=0; for(char c:s) w=R[c-'0'][w]; return w; }
static void pprint(const Poly&p,int d){ printf("["); for(int i=0;i<=d/2;i++) printf("%s%d",i?", ":"",p.c[i]); printf("]"); }

int main(int argc,char**argv){
    char t=argv[1][0]; n=atoi(argv[2]);
    for(int i=0;i<n;i++) for(int j=0;j<n;j++) Cart[i][j]=(i==j)?2:0;
    auto edge=[&](int i,int j){Cart[i][j]=Cart[j][i]=-1;};
    for(int i=0;i+1<=n-2;i++) edge(i,i+1);
    if(t=='D') edge(n-3,n-1); else if(t=='E') edge(2,n-1);
    Mat e; memset(&e,0,sizeof e); for(int i=0;i<n;i++) e.m[i][i]=1;
    els.push_back(e); idx[e]=0; len.push_back(0); par.push_back(-1); parlet.push_back(-1);
    for(int s=0;s<n;s++) R[s].push_back(-1);
    for(size_t w=0; w<els.size(); w++){
        Mat cur=els[w];
        for(int s=0;s<n;s++){
            Mat v=rmul(cur,s); auto it=idx.find(v); int vi;
            if(it==idx.end()){ vi=els.size(); idx[v]=vi; els.push_back(v); len.push_back(len[w]+1); par.push_back(w); parlet.push_back(s); for(int q=0;q<n;q++) R[q].push_back(-1); }
            else vi=it->second;
            R[s][w]=vi;
        }
    }
    N=els.size(); fprintf(stderr,"|W|=%d longest=%d\n",N,len[N-1]);
    for(int s=0;s<n;s++) L[s].assign(N,-1);
    Rdes.assign(N,0); Ldes.assign(N,0);
    for(int w=0;w<N;w++){
        for(int s=0;s<n;s++){ int v=idx[lmul(s,els[w])]; L[s][w]=v; if(len[v]<len[w]) Ldes[w]|=1u<<s; if(rowneg(els[w],s)) Rdes[w]|=1u<<s; if((len[R[s][w]]<len[w]) != rowneg(els[w],s)){fprintf(stderr,"descent mismatch\n");return 1;} }
    }
    // Bruhat ideals
    WPL=(N+63)/64; low.assign((size_t)N*WPL,0); setb(0,0);
    for(int w=1;w<N;w++){
        int s=__builtin_ctz(Rdes[w]); int v=R[s][w];
        memcpy(&low[(size_t)w*WPL],&low[(size_t)v*WPL],WPL*8);
        for(size_t q=0;q<WPL;q++){ uint64_t b=low[(size_t)v*WPL+q]; while(b){ int x=q*64+__builtin_ctzll(b); b&=b-1; setb(w,R[s][x]); } }
    }
    fprintf(stderr,"Bruhat done\n");
    memset(&ZERO,0,sizeof ZERO);
    ext.resize(N); pol.resize(N); muL.resize(N);
    long long nExt=0;
    for(int w=0;w<N;w++){
        // extremal x
        for(size_t q=0;q<WPL;q++){ uint64_t b=low[(size_t)w*WPL+q]; while(b){ int x=q*64+__builtin_ctzll(b); b&=b-1; if((Ldes[w]&~Ldes[x])==0 && (Rdes[w]&~Rdes[x])==0) ext[w].push_back(x);} }
        pol[w].resize(ext[w].size()); nExt+=ext[w].size();
        if(w==0){ memset(&pol[0][0],0,sizeof(Poly)); pol[0][0].c[0]=1; continue; }
        int s=__builtin_ctz(Ldes[w]); int v=L[s][w];
        for(size_t i=0;i<ext[w].size();i++){
            int x=ext[w][i]; Poly p; memset(&p,0,sizeof p);
            if(x==w){ p.c[0]=1; pol[w][i]=p; continue; }
            int sx=L[s][x]; // sx<x since s in L(w) subset L(x)
            const Poly*a=lookup(sx,v); for(int k=0;k<MAXD;k++) p.c[k]+=a->c[k];
            const Poly*b=lookup(x,v); for(int k=0;k+1<MAXD;k++) p.c[k+1]+=b->c[k];
            for(auto&zm:muL[v]){ int z=zm.first; if(!((Ldes[z]>>s)&1)) continue; if(!leq(x,z)) continue; int sh=(len[w]-len[z])/2; const Poly*c=lookup(x,z); for(int k=0;k+sh<MAXD;k++) p.c[k+sh]-=zm.second*c->c[k]; }
            int d=len[w]-len[x]; for(int k=(d-1)/2+1;k<MAXD;k++) if(p.c[k]){fprintf(stderr,"degree violation\n");return 1;}
            pol[w][i]=p;
        }
        // mu list of w: covers plus extremal odd-gap nonzero
        for(size_t q=0;q<WPL;q++){ uint64_t b=low[(size_t)w*WPL+q]; while(b){ int z=q*64+__builtin_ctzll(b); b&=b-1; if(len[z]+1==len[w]) muL[w].push_back({z,1}); } }
        for(size_t i=0;i<ext[w].size();i++){ int z=ext[w][i]; int d=len[w]-len[z]; if(d>1&&d%2==1){ int m=pol[w][i].c[(d-1)/2]; if(m) muL[w].push_back({z,m}); } }
        if(w%5000==0) fprintf(stderr,"w=%d len=%d ext so far=%lld\n",w,len[w],nExt);
    }
    // tally
    long long covers=0; map<int,long long> hist; int maxmu=0; int maxcoef=0;
    for(int w=0;w<N;w++){ for(auto&zm:muL[w]){ if(len[w]-len[zm.first]==1) covers++; else hist[zm.second]++; if(zm.second>maxmu) maxmu=zm.second; }
        for(auto&p:pol[w]) for(int k=0;k<MAXD;k++) if(p.c[k]>maxcoef) maxcoef=p.c[k]; }
    printf("RESULT type=%c%d |W|=%d extremal_pairs=%lld max_coefficient=%d covers=%lld max_mu=%d\n",t,n,N,nExt,maxcoef,covers,maxmu);
    printf("noncover nonzero mu histogram:"); for(auto&h:hist) printf(" mu=%d: %lld;",h.first,h.second); printf("\n");
    // print all pairs with mu>=5 (or the maximum in D)
    int thr = (t=='E')?5:maxmu;
    for(int w=0;w<N;w++) for(auto&zm:muL[w]) if(len[w]-len[zm.first]>1 && zm.second>=thr){ int x=zm.first; printf("mu=%d x=%s (len %d) w=%s (len %d) P=",zm.second,word(x).c_str(),len[x],word(w).c_str(),len[w]); pprint(*lookup(x,w),len[w]-len[x]); printf("\n"); }
    for(int a=3;a+1<argc;a+=2){ int x=elem(argv[a]), w=elem(argv[a+1]); int d=len[w]-len[x]; printf("PAIR x=%s(len %d) w=%s(len %d) leq=%d P=",argv[a],len[x],argv[a+1],len[w],(int)leq(x,w)); pprint(*lookup(x,w),d); printf(" mu=%d\n",mu(x,w)); }
    return 0;
}
