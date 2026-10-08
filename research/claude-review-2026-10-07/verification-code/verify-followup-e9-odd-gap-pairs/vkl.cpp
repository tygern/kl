// Independent verifier: ordinary KL polynomials P_{x,W} for a target W in a
// simply laced Coxeter group, computed via the left-coset quotient W/W_J with
// J = R(W), using the identity P_{x w_J, w w_J} = Q_{x,w} (x,w in W^J) and the
// recursion derived from the ordinary left-multiplication recursion:
//   Q_{x,z} = [sx<x](Q_{sx,v}+q Q_{x,v}) + [sx>x, sx in W^J](q Q_{sx,v}+Q_{x,v})
//           + [sx notin W^J](1+q) Q_{x,v}
//           - sum_{y in W^J, sy<y, y<v, mu(y,v)!=0} mu(y,v) q^{(l(z)-l(y))/2} Q_{x,y}
// where v = sz, s in L(z).
// Usage: ./vkl <type> <word letters...>   type in {E6,E7,E8,E9,E10,D4,A3,A4}
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <tuple>
#include <cctype>
using namespace std;
typedef long long ll;
typedef vector<ll> Poly;

int n; int A[16][16];
struct Elt { int16_t M[16][16]; int16_t Mi[16][16]; };
// columns: M[i][j] = coeff of alpha_i in x(alpha_j)

void setup(const string& t){
    memset(A,0,sizeof A);
    vector<pair<int,int>> E;
    if(t[0]=='E'){ n=stoi(t.substr(1)); for(int i=0;i+1<n-1;i++)E.push_back({i,i+1}); E.push_back({2,n-1}); }
    else if(t[0]=='D'){ n=stoi(t.substr(1)); for(int i=0;i+1<n-1;i++)E.push_back({i,i+1}); E.push_back({n-3,n-1}); }
    else if(t[0]=='A'){ n=stoi(t.substr(1)); for(int i=0;i+1<n;i++)E.push_back({i,i+1}); }
    else { fprintf(stderr,"bad type\n"); exit(1);}
    for(int i=0;i<n;i++)A[i][i]=2;
    for(auto [a,b]:E){A[a][b]=-1;A[b][a]=-1;}
}
Elt identity(){ Elt e; memset(&e,0,sizeof e); for(int i=0;i<n;i++){e.M[i][i]=1;e.Mi[i][i]=1;} return e; }
int pairing(int s, const int16_t* col /*col[i]*/){ int p=0; for(int i=0;i<n;i++)p+=A[s][i]*col[i]; return p; }
// left multiply: (s x)(alpha_j) = x(alpha_j) - <alpha_s^vee, x(alpha_j)> alpha_s ; inverse: x^{-1} s
Elt lmul(const Elt& x,int s){
    Elt y=x;
    for(int j=0;j<n;j++){ int p=0; for(int i=0;i<n;i++)p+=A[s][i]*x.M[i][j]; y.M[s][j]-=p; }
    for(int j=0;j<n;j++){ if(A[s][j]!=0) for(int i=0;i<n;i++) y.Mi[i][j]-=A[s][j]*x.Mi[i][s]; }
    return y;
}
// right multiply: (x s)(alpha_j) = x(alpha_j) - a_{sj} x(alpha_s); inverse: s x^{-1}
Elt rmul(const Elt& x,int s){
    Elt y=x;
    for(int j=0;j<n;j++){ if(A[s][j]!=0) for(int i=0;i<n;i++) y.M[i][j]-=A[s][j]*x.M[i][s]; }
    for(int j=0;j<n;j++){ int p=0; for(int i=0;i<n;i++)p+=A[s][i]*x.Mi[i][j]; y.Mi[s][j]-=p; }
    return y;
}
bool colneg(const int16_t (*M)[16], int j){ bool nz=false; for(int i=0;i<n;i++){ if(M[i][j]>0)return false; if(M[i][j]<0)nz=true;} return nz; }
bool colpos(const int16_t (*M)[16], int j){ bool nz=false; for(int i=0;i<n;i++){ if(M[i][j]<0)return false; if(M[i][j]>0)nz=true;} return nz; }
bool rdesc(const Elt& x,int s){ return colneg(x.M,s); }
bool ldesc(const Elt& x,int s){ return colneg(x.Mi,s); }
bool isid(const Elt& x){ for(int i=0;i<n;i++)for(int j=0;j<n;j++) if(x.M[i][j]!=(i==j))return false; return true; }
vector<int> redword(Elt x){ vector<int> w; while(!isid(x)){ int s=-1; for(int t=0;t<n;t++)if(rdesc(x,t)){s=t;break;} if(s<0){fprintf(stderr,"no descent?\n");exit(1);} w.push_back(s); x=rmul(x,s);} reverse(w.begin(),w.end()); return w; }
int len(const Elt& x){ return (int)redword(x).size(); }
struct Key { array<int16_t,256> a; bool operator==(const Key&o)const{return a==o.a;} };
Key key(const Elt& x){ Key k; k.a.fill(0); for(int i=0;i<n;i++)for(int j=0;j<n;j++)k.a[i*16+j]=x.M[i][j]; return k; }
struct KH { size_t operator()(const Key& k)const{ size_t h=1469598103934665603ULL; for(int i=0;i<256;i++){ h^=(uint16_t)k.a[i]; h*=1099511628211ULL;} return h; } };
string polystr(const Poly& p){ string s; bool first=true; for(size_t i=0;i<p.size();i++){ if(!p[i])continue; if(!first)s+=" + "; first=false; if(i==0)s+=to_string(p[i]); else s+=to_string(p[i])+"q^"+to_string(i);} if(first)s="0"; return s; }
void trim(Poly& p){ while(!p.empty()&&p.back()==0)p.pop_back(); }

int main(int argc,char**argv){
    if(argc<3){fprintf(stderr,"usage\n");return 1;}
    setup(argv[1]);
    vector<int> word; for(int i=2;i<argc;i++)word.push_back(atoi(argv[i]));
    Elt W=identity(); for(int s:word)W=rmul(W,s);
    int lW=len(W);
    printf("word:"); for(int s:word)printf(" %d",s); printf("\nlength %d (word length %zu)%s\n",lW,word.size(), lW==(int)word.size()?"":"  [NOT REDUCED]");
    vector<int> J; for(int s=0;s<n;s++)if(rdesc(W,s))J.push_back(s);
    printf("R(W)={"); for(int s:J)printf("%d,",s); printf("}  L(W)={"); for(int s=0;s<n;s++)if(ldesc(W,s))printf("%d,",s); printf("}\n");
    // is W a reflection? W^2 = 1 and W = 1 - beta beta^T A (rank-one difference): check involution and -1 eigenvector
    { Elt W2=W; vector<int> rw=redword(W); for(int s:rw)W2=rmul(W2,s); printf("involution: %s\n", isid(W2)?"yes":"no"); }
    // w_J: longest element of W_J, by BFS over W_J (must be finite)
    Elt wJ=identity();
    { unordered_map<Key,int,KH> seen; vector<Elt> q{identity()}; seen[key(q[0])]=0; size_t head=0; Elt best=q[0]; int bl=0;
      while(head<q.size()){ Elt x=q[head++]; for(int s:J){ Elt y=rmul(x,s); Key k=key(y); if(!seen.count(k)){ seen[k]=1; q.push_back(y);} } if(q.size()>100000){fprintf(stderr,"W_J too big\n");return 1;} }
      for(auto& x:q){ int l=len(x); if(l>bl){bl=l;best=x;} }
      wJ=best; printf("|W_J| = %zu, l(w_J)=%d\n", q.size(), bl); }
    Elt Wp=W; { vector<int> rw=redword(wJ); for(int s:rw)Wp=rmul(Wp,s); }  // W' = W w_J
    int lWp=len(Wp);
    for(int t:J) if(!colpos(Wp.M,t)){fprintf(stderr,"W' not min rep\n");return 1;}
    printf("W' = W w_J has length %d (expect %d)\n", lWp, lW-len(wJ));
    // quotient lower interval of W' via suffix DP with left multiplication and projection
    vector<int> rwp=redword(Wp);
    vector<Elt> elts; unordered_map<Key,int,KH> idx;
    auto addE=[&](const Elt& x)->int{ Key k=key(x); auto it=idx.find(k); if(it!=idx.end())return it->second; int id=elts.size(); elts.push_back(x); idx[k]=id; return id; };
    addE(identity());
    for(int i=(int)rwp.size()-1;i>=0;i--){ int s=rwp[i]; size_t cur=elts.size();
        for(size_t a=0;a<cur;a++){ Elt y=lmul(elts[a],s); bool inWJ=true; for(int t:J) if(!colpos(y.M,t)){inWJ=false;break;} if(inWJ) addE(y); }
    }
    int N=elts.size();
    printf("quotient ideal size %d\n", N);
    vector<int> L(N); for(int i=0;i<N;i++)L[i]=len(elts[i]);
    ll full=0; { // full ideal size = sum over quotient elements of |W_J| ... only if all cosets fully included: [e,W] is a union of cosets since J subset R(W). so full = N*|W_J|
    }
    int idWp=idx.count(key(Wp))?idx[key(Wp)]:-1; if(idWp<0){fprintf(stderr,"W' missing\n");return 1;}
    if(L[idWp]!=lWp){fprintf(stderr,"length mismatch\n");return 1;}
    // left mult table: lm[x][s] = index of sx if in W^J (and in ideal), -1 if sx notin W^J, -2 if not in ideal
    vector<array<int,16>> lm(N);
    for(int a=0;a<N;a++) for(int s=0;s<n;s++){ Elt y=lmul(elts[a],s); bool inWJ=true; for(int t:J) if(!colpos(y.M,t)){inWJ=false;break;} if(!inWJ){lm[a][s]=-1;continue;} auto it=idx.find(key(y)); lm[a][s]= it==idx.end()? -2 : it->second; }
    // ldesc table
    vector<uint16_t> LD(N,0); for(int a=0;a<N;a++) for(int s=0;s<n;s++) if(ldesc(elts[a],s)) LD[a]|=(1<<s);
    // sanity: left descent s of x in W^J implies sx in W^J and in ideal
    for(int a=0;a<N;a++) for(int s=0;s<n;s++) if(LD[a]>>s&1){ if(lm[a][s]<0){fprintf(stderr,"sanity fail\n");return 1;} if(L[lm[a][s]]!=L[a]-1){fprintf(stderr,"len sanity fail\n");return 1;} }
    // KL recursion with memo
    unordered_map<int, vector<pair<int,Poly>>> memo;
    vector<Poly> acc(N); vector<char> touched(N,0); vector<int> tl;
    ll pairs=0; int computed=0;
    function<const vector<pair<int,Poly>>&(int)> Q = [&](int z)->const vector<pair<int,Poly>>&{
        auto it=memo.find(z); if(it!=memo.end())return it->second;
        vector<pair<int,Poly>> res;
        if(L[z]==0){ res.push_back({z,Poly{1}}); return memo[z]=res; }
        int s=-1; for(int t=0;t<n;t++) if(LD[z]>>t&1){s=t;break;}
        int v=lm[z][s];
        // need Q(v) and the mu-terms' Q(y) computed first (recursion), collect them
        vector<pair<int,Poly>> Qv = Q(v); // copy (memo may rehash)
        vector<tuple<int,int,ll>> corr; // (y, shift, mu)
        for(auto& [y,p]: Qv){ if(y==v)continue; int d=L[v]-L[y]; if(d<=0||d%2==0)continue; if(!((LD[y]>>s&1) || lm[y][s]==-1))continue; /* sZ<Z for Z=y w_J iff sy<y or sy notin W^J */ int dg=(d-1)/2; if((int)p.size()>dg && p[dg]!=0) corr.push_back({y,(d+1)/2,p[dg]}); }
        for(auto& [y,sh,mu]:corr) Q(y);
        // accumulate
        auto addto=[&](int x,const Poly& p,int shift,ll scale){ if(!touched[x]){touched[x]=1;tl.push_back(x);} Poly& a=acc[x]; if(a.size()<p.size()+shift)a.resize(p.size()+shift,0); for(size_t i=0;i<p.size();i++)a[i+shift]+=scale*p[i]; };
        for(auto& [x,p]:Qv){
            int sx=lm[x][s];
            if(LD[x]>>s&1){ addto(x,p,1,1); addto(sx,p,1,1); }
            else if(sx>=0){ addto(x,p,0,1); addto(sx,p,0,1); }
            else if(sx==-1){ addto(x,p,0,1); addto(x,p,1,1); }
            else { fprintf(stderr,"sx not in ideal: lifting violated\n"); exit(1);} }
        for(auto& [y,sh,mu]:corr){ const auto& Qy=memo[y]; for(auto& [x,p]:Qy) addto(x,p,sh,-mu); }
        for(int x:tl){ trim(acc[x]); if(!acc[x].empty()){ res.push_back({x,acc[x]}); } acc[x].clear(); touched[x]=0; }
        tl.clear();
        sort(res.begin(),res.end(),[](auto&a,auto&b){return a.first<b.first;});
        // sanity: constant terms 1, nonneg coefficients, degree bound
        for(auto& [x,p]:res){ if(p[0]!=1){fprintf(stderr,"const term !=1 at z=%d x=%d: %s\n",z,x,polystr(p).c_str());exit(1);} for(ll c:p) if(c<0){fprintf(stderr,"negative coeff\n");exit(1);} int d=L[z]-L[x]; if(x!=z && 2*((int)p.size()-1)>d-1){ auto zw=redword(elts[z]); auto xw=redword(elts[x]); fprintf(stderr,"degree bound violated: z="); for(int t:zw)fprintf(stderr,"%d",t); fprintf(stderr," x="); for(int t:xw)fprintf(stderr,"%d",t); fprintf(stderr," P=%s  v=",polystr(p).c_str()); for(int t:redword(elts[v]))fprintf(stderr,"%d",t); fprintf(stderr," s=%d corr:",s); for(auto&[y,sh,mu]:corr){fprintf(stderr," (y="); for(int t:redword(elts[y]))fprintf(stderr,"%d",t); fprintf(stderr,",sh=%d,mu=%lld)",sh,mu);} fprintf(stderr,"\n"); exit(1);} }
        pairs+=res.size(); computed++;
        return memo[z]=res;
    };
    auto t0=chrono::steady_clock::now();
    const auto& QW = Q(idWp);
    double secs=chrono::duration<double>(chrono::steady_clock::now()-t0).count();
    printf("KL done: %d quotient elements computed, %lld pairs, %.1f s; |supp Q_{.,W'}| = %zu (should equal quotient size %d)\n", computed,pairs,secs,QW.size(),N);
    // Report P_{x,W} for x = u w_J with u in W^J of small length (all u with l(u)<=3), plus requested ones
    // Also report the FC bottoms: here we print all u with l(u)<=2 and the specific u's given via env var XWORDS (semicolon separated words)
    map<int,Poly> QWmap; for(auto& [x,p]:QW)QWmap[x]=p;
    auto report=[&](const vector<int>& uword){ Elt u=identity(); for(int s:uword)u=rmul(u,s); auto it=idx.find(key(u)); if(it==idx.end()){printf("u="); for(int s:uword)printf("%d",s); printf(" : not in quotient ideal / not min rep\n"); return;} int ui=it->second; Poly p=QWmap.count(ui)?QWmap[ui]:Poly{}; int gap=lW-(L[ui]+len(wJ)); int dg=(gap-1)/2; ll mu= (gap%2==1 && (int)p.size()>dg)? p[dg]:0; printf("x = u*w_J with u="); for(int s:uword)printf("%d",s); printf(" l(x)=%d gap %d : P = %s   mu = %lld%s\n", L[ui]+len(wJ), gap, polystr(p).c_str(), mu, gap%2==0?" (even gap)":""); };
    const char* xw=getenv("UWORDS");
    if(xw){ string s(xw); stringstream ss(s); string tok; while(getline(ss,tok,';')){ vector<int> u; for(char c:tok) if(isdigit(c)) u.push_back(c-'0'); report(u);} }
    // histogram of mu over all x in quotient with odd gap
    map<ll,int> hist; for(auto& [x,p]:QW){ int gap=lW-(L[x]+len(wJ)); if(gap%2==0||gap<=0)continue; int dg=(gap-1)/2; ll mu=(int)p.size()>dg?p[dg]:0; if(mu)hist[mu]++; }
    printf("mu histogram over quotient reps with odd gap:"); for(auto&[m,c]:hist)printf(" mu=%lld:%d",m,c); printf("\n");
    if(getenv("DUMPIDEAL")){ FILE* f=fopen(getenv("DUMPIDEAL"),"w"); for(int a=0;a<N;a++){ for(int i=0;i<n;i++)for(int j=0;j<n;j++)fprintf(f,"%d%c",(int)elts[a].M[i][j], (i==n-1&&j==n-1)?'\n':','); } fclose(f); }
    if(getenv("DUMP")){ // dump all P_{u w_J, W} as u-word : poly
        for(auto& [x,p]:QW){ vector<int> uw=redword(elts[x]); printf("DUMP u="); for(int s:uw)printf("%d",s); printf(" : %s\n",polystr(p).c_str()); }
    }
    return 0;
}
