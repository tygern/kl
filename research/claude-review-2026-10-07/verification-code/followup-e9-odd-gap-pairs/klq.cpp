// Independent parabolic (maximal-coset) Kazhdan-Lusztig computation for
// simply laced Coxeter groups, elements represented by their integer matrix
// on the span of the simple roots.  Written for the referee follow-up on the
// affine E8 reflection family.  Exact integer arithmetic throughout (int64
// with overflow checks).
//
// Conventions: P_{x,w} in q, standard KL polynomials; right multiplication by
// s acts on columns, left multiplication on rows.  For I subset R(w),
// w = w' w_I with w' in W^I (minimal right-coset reps).  We compute
// Pq(x,y) := P_{x w_I, y w_I} for x,y in W^I using the recursion derived from
// the left ideal H C'_{w_I} (see notes in the report).
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <climits>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <tuple>
#include <random>
#include <chrono>
#include <cctype>
using namespace std;
typedef long long ll;
typedef array<int8_t,81> Mat;
static int N; // rank
static vector<vector<int>> ADJ;
struct MatHash{ size_t operator()(const Mat&m)const{ uint64_t h=1469598103934665603ULL; for(int i=0;i<N*9;i++){h^=(uint8_t)m[i];h*=1099511628211ULL;} return h;} };
static Mat identity(){Mat m; m.fill(0); for(int j=0;j<N;j++)m[j*9+j]=1; return m;}
static inline int8_t chk(int v){ if(v>120||v<-120){fprintf(stderr,"int8 overflow\n");exit(3);} return (int8_t)v; }
static Mat rmul(const Mat&w,int s){ Mat r=w; for(int i=0;i<N;i++) r[s*9+i]=chk(-w[s*9+i]); for(int t:ADJ[s]) for(int i=0;i<N;i++) r[t*9+i]=chk(w[t*9+i]+w[s*9+i]); return r; }
static Mat lmul(const Mat&w,int s){ Mat r=w; for(int j=0;j<N;j++){ int p=2*w[j*9+s]; for(int t:ADJ[s]) p-=w[j*9+t]; r[j*9+s]=chk(w[j*9+s]-p);} return r; }
static bool negcol(const Mat&w,int j){ bool nz=false; for(int i=0;i<N;i++){ if(w[j*9+i]>0) return false; if(w[j*9+i]<0) nz=true;} if(!nz){fprintf(stderr,"zero column\n");exit(3);} return true; }
static bool poscol(const Mat&w,int j){ return !negcol(w,j); }
static vector<int> rword(Mat w){ vector<int> out; Mat e=identity(); while(w!=e){ int s=-1; for(int j=0;j<N;j++) if(negcol(w,j)){s=j;break;} if(s<0){fprintf(stderr,"no descent\n");exit(3);} w=rmul(w,s); out.push_back(s);} reverse(out.begin(),out.end()); return out; }
static Mat fromword(const vector<int>&wd){ Mat w=identity(); for(int s:wd) w=rmul(w,s); return w; }
static Mat inverse(const Mat&w){ vector<int> wd=rword(w); reverse(wd.begin(),wd.end()); return fromword(wd); }
static int Rmask(const Mat&w){int m=0; for(int j=0;j<N;j++) if(negcol(w,j)) m|=1<<j; return m;}
static int Lmask(const Mat&w){ return Rmask(inverse(w)); }
static int pairing(const vector<int>&a,const vector<int>&b){ // (a,b) with Cartan matrix
  int p=0; for(int i=0;i<N;i++){ p+=2*a[i]*b[i]; for(int t:ADJ[i]) p-=a[i]*b[t]; } return p; }
static Mat reflection(const vector<int>&beta){ // r_beta(alpha_j) = alpha_j - (alpha_j,beta) beta
  Mat m=identity(); for(int j=0;j<N;j++){ vector<int> aj(N,0); aj[j]=1; int p=pairing(aj,beta); for(int i=0;i<N;i++) m[j*9+i]=chk((i==j)-p*beta[i]); } return m; }
// FC test (Stembridge): w is FC iff no reduced word contains sts; recursive: all right weak predecessors FC and no two adjacent right descents... Actually
// the standard criterion used (as in the authors' code): w FC iff for all s in R(w): FC(ws) and no t~s with both s,t in R(w) -- this is the test that
// the right descent set is commuting at every stage of the right weak order; equivalent to FC (Stembridge Prop 2.? : w FC iff every element of
// the right weak lower interval has commuting right descents is FALSE in general?)  We instead use the braid-avoidance definition directly via
// heap: w is FC iff for every element u <=_R w (right weak order), R(u) is a commuting set.  (Stembridge, FC elements of Coxeter groups, Prop 2.1/ Thm 3.1:
// w is FC iff no reduced word contains a factor sts, iff no u in the right weak lower interval has a non-commuting pair of right descents --
// if u has s,t non-commuting in R(u) then u = u' sts reduced so w has a reduced word with sts.)  Conversely a reduced word of w with factor sts gives
// such a u.  So the criterion is exact.
static unordered_map<Mat,bool,MatHash> FCMEMO;
static bool isFC(const Mat&w){ auto it=FCMEMO.find(w); if(it!=FCMEMO.end()) return it->second; int R=Rmask(w); bool ok=true; for(int s=0;s<N&&ok;s++) if(R>>s&1){ for(int t:ADJ[s]) if(R>>t&1) ok=false; if(ok&&!isFC(rmul(w,s))) ok=false; } FCMEMO[w]=ok; return ok; }

// ---------------- polynomial helpers --------------
typedef vector<ll> Poly;
static void trim(Poly&p){ while(!p.empty()&&p.back()==0) p.pop_back(); }
static inline ll addchk(ll a,ll b){ ll c; if(__builtin_add_overflow(a,b,&c)){fprintf(stderr,"overflow\n");exit(3);} return c; }
static void addto(Poly&a,const Poly&b,int shift,ll scale){ if(a.size()<b.size()+shift) a.resize(b.size()+shift,0); for(size_t i=0;i<b.size();i++){ ll t; if(__builtin_mul_overflow(b[i],scale,&t)){fprintf(stderr,"overflow\n");exit(3);} a[i+shift]=addchk(a[i+shift],t);} }
static string pstr(const Poly&p){ if(p.empty()) return "0"; string s; for(size_t i=0;i<p.size();i++){ if(p[i]==0) continue; if(!s.empty()) s+=" + "; s+=to_string(p[i]); if(i>0) s+="q^"+to_string(i);} return s; }

// ---------------- quotient structure --------------
struct Quot {
  int I; // mask
  Mat wI; int lwI;
  vector<Mat> el; vector<int> len; vector<int> Lq; // Lq mask
  vector<array<int,9>> left; // index, -1 NOTMIN (sx = xt), -2 OUT (not in Q, sx>x in W^I)
  unordered_map<Mat,int,MatHash> idx;
  // buckets by Lq mask: list of indices sorted by length
  vector<vector<int>> bucket;
  bool inWI(const Mat&x){ for(int t=0;t<N;t++) if(I>>t&1) if(negcol(x,t)) return false; return true; }
  void build(const Mat&w, int Imask, bool verbose){
    I=Imask;
    // longest element of W_I by greedy ascent
    wI=identity(); while(true){ bool moved=false; for(int t=0;t<N;t++) if(I>>t&1){ if(!negcol(wI,t)){ wI=rmul(wI,t); moved=true; } } if(!moved) break; }
    lwI=rword(wI).size();
    Mat wp=w; for(int t:rword(wI)) wp=rmul(wp,t); // w' = w wI^{-1} = w wI
    if(!inWI(wp)){fprintf(stderr,"w' not minimal: I not subset of R(w)?\n");exit(3);}
    vector<int> wd=rword(wp);
    unordered_map<Mat,int,MatHash> seen; vector<Mat> cur; Mat e=identity(); seen[e]=0; cur.push_back(e);
    for(int k=(int)wd.size()-1;k>=0;k--){ int s=wd[k]; size_t sz=cur.size(); for(size_t i=0;i<sz;i++){ Mat y=lmul(cur[i],s); if(!inWI(y)) continue; if(seen.count(y)) continue; seen[y]=cur.size(); cur.push_back(y);} if(verbose) fprintf(stderr,"suffix %d size %zu\n",(int)wd.size()-k,cur.size()); }
    // lengths and sort
    vector<pair<int,int>> order; order.reserve(cur.size());
    // length: compute via reduced word of x (cheap enough) -- use rword length
    vector<int> L(cur.size()); for(size_t i=0;i<cur.size();i++){ L[i]=rword(cur[i]).size(); order.push_back({L[i],(int)i}); }
    sort(order.begin(),order.end());
    el.resize(cur.size()); len.resize(cur.size());
    for(size_t i=0;i<order.size();i++){ el[i]=cur[order[i].second]; len[i]=order[i].first; idx[el[i]]=i; }
    left.resize(el.size()); Lq.assign(el.size(),0);
    for(size_t i=0;i<el.size();i++) for(int s=0;s<N;s++){ Mat y=lmul(el[i],s); if(!inWI(y)){ left[i][s]=-1; Lq[i]|=1<<s; continue;} auto it=idx.find(y); if(it==idx.end()){ left[i][s]=-2; continue;} left[i][s]=it->second; if(len[it->second]<len[i]) Lq[i]|=1<<s; else if(len[it->second]!=len[i]+1){fprintf(stderr,"length mismatch\n");exit(3);} }
    bucket.assign(1<<N,{}); for(size_t i=0;i<el.size();i++) bucket[Lq[i]].push_back(i); // already length-sorted
  }
};

// ---------------- KL in quotient --------------
struct KLQ {
  Quot&Q; std::mt19937 rng; bool randomize=false;
  vector<Poly> polys; unordered_map<string,int> pid; // hash-consing via byte string
  vector<vector<pair<int,int>>> lst; vector<char> done; // per y: sorted (x, polyid)
  size_t pairs=0; size_t computed=0;
  KLQ(Quot&q):Q(q){ lst.resize(Q.el.size()); done.assign(Q.el.size(),0); }
  int id(const Poly&p){ string key((const char*)p.data(),p.size()*sizeof(ll)); auto it=pid.find(key); if(it!=pid.end()) return it->second; int i=polys.size(); polys.push_back(p); pid[key]=i; return i; }
  // P_{x,z} lookup with descent reduction
  const Poly* get(int x,int z){ static Poly zero; if(Q.len[x]>Q.len[z]) return &zero; int Lz=Q.Lq[z];
    while(true){ int d=Lz&~Q.Lq[x]; if(!d) break; int t=__builtin_ctz(d); int nx=Q.left[x][t]; if(nx<0){ if(nx==-2) return &zero; fprintf(stderr,"reduction NOTMIN impossible\n"); exit(3);} x=nx; if(Q.len[x]>Q.len[z]) return &zero; }
    auto&v=lst[z]; auto it=lower_bound(v.begin(),v.end(),make_pair(x,INT_MIN)); if(it==v.end()||it->first!=x) return &zero; return &polys[it->second]; }
  void compute(int y){ if(done[y]) return; if(Q.len[y]==0){ lst[y]={{y,id(Poly{1})}}; done[y]=1; return; }
    int s=-1; if(randomize){ vector<int> opts; for(int c=0;c<N;c++) if((Q.Lq[y]>>c&1)&&Q.left[y][c]>=0) opts.push_back(c); s=opts[rng()%opts.size()]; }
    else for(int c=0;c<N;c++) if((Q.Lq[y]>>c&1)&&Q.left[y][c]>=0){ if(s<0) s=c; if(done[Q.left[y][c]]){s=c;break;} }
    if(s<0){fprintf(stderr,"no quotient descent\n");exit(3);} int u=Q.left[y][s]; compute(u);
    // corrections
    vector<tuple<int,ll,int>> corr; // z, mu, shift
    for(auto&pr:lst[u]){ int z=pr.first; if(z==u) continue; if(!(Q.Lq[z]>>s&1)) continue; int gap=Q.len[u]-Q.len[z]; if(gap%2==0) continue; const Poly&p=polys[pr.second]; size_t top=(gap-1)/2; if(p.size()>top+1){fprintf(stderr,"degree violation\n");exit(3);} if(p.size()==top+1&&p[top]!=0) corr.emplace_back(z,p[top],(Q.len[y]-Q.len[z])/2); }
    for(int t=0;t<N;t++) if((Q.Lq[u]>>t&1)&&Q.left[u][t]>=0){ int z=Q.left[u][t]; if(Q.Lq[z]>>s&1) corr.emplace_back(z,1,1); }
    for(auto&c:corr) compute(std::get<0>(c));
    vector<pair<int,int>> out;
    int Ly=Q.Lq[y]; int full=(1<<N)-1; int freebits=full&~Ly;
    // iterate supersets D of Ly
    for(int sub=freebits;;sub=(sub-1)&freebits){ int D=Ly|sub; auto&b=Q.bucket[D];
      for(int x:b){ if(Q.len[x]>=Q.len[y]) break; Poly p; int sx=Q.left[x][s];
        if(sx==-1){ const Poly*a=get(x,u); addto(p,*a,0,1); addto(p,*a,1,1); }
        else { if(sx<0||Q.len[sx]>=Q.len[x]){fprintf(stderr,"candidate without s-descent\n");exit(3);} addto(p,*get(sx,u),0,1); addto(p,*get(x,u),1,1); }
        for(auto&c:corr){ int z=std::get<0>(c); if(Q.len[x]>Q.len[z]) continue; const Poly*pz=get(x,z); if(!pz->empty()) addto(p,*pz,std::get<2>(c),-std::get<1>(c)); }
        trim(p); if(p.empty()) continue; int gap=Q.len[y]-Q.len[x]; if((int)p.size()-1>(gap-1)/2){fprintf(stderr,"degree bound violated y=%d x=%d\n",y,x);exit(3);} for(ll c:p) if(c<0){fprintf(stderr,"negative coefficient y=%d x=%d: %s\n",y,x,pstr(p).c_str());exit(3);} if(p[0]!=1){fprintf(stderr,"constant term !=1\n");exit(3);} out.push_back({x,id(p)}); }
      if(sub==0) break; }
    out.push_back({y,id(Poly{1})}); sort(out.begin(),out.end()); lst[y]=out; done[y]=1; pairs+=out.size(); computed++;
    if(computed%2000==0) fprintf(stderr,"computed %zu elements, %zu pairs, %zu polys\n",computed,pairs,polys.size()); }
};

// ---------------- ordinary naive KL (for validation) --------------
struct Ord {
  vector<Mat> el; vector<int> len; unordered_map<Mat,int,MatHash> idx; vector<array<int,9>> left; vector<int> L;
  map<pair<int,int>,Poly> memo;
  void build(const Mat&w){ vector<int> wd=rword(w); unordered_map<Mat,int,MatHash> seen; vector<Mat> cur; Mat e=identity(); seen[e]=0; cur.push_back(e);
    for(int s:wd){ size_t sz=cur.size(); for(size_t i=0;i<sz;i++){ Mat y=rmul(cur[i],s); if(seen.count(y)) continue; seen[y]=cur.size(); cur.push_back(y);} }
    vector<pair<int,int>> order; for(size_t i=0;i<cur.size();i++) order.push_back({(int)rword(cur[i]).size(),(int)i}); sort(order.begin(),order.end());
    el.resize(cur.size()); len.resize(cur.size()); for(size_t i=0;i<order.size();i++){ el[i]=cur[order[i].second]; len[i]=order[i].first; idx[el[i]]=i; }
    left.resize(el.size()); L.assign(el.size(),0); for(size_t i=0;i<el.size();i++) for(int s=0;s<N;s++){ Mat y=lmul(el[i],s); auto it=idx.find(y); if(it==idx.end()){left[i][s]=-2; continue;} left[i][s]=it->second; if(len[it->second]<len[i]) L[i]|=1<<s; } }
  bool leq(int x,int y){ while(x!=y){ if(len[x]>=len[y]) return false; int s=__builtin_ctz(L[y]); if(L[x]>>s&1) x=left[x][s]; y=left[y][s]; } return true; }
  Poly P(int x,int y){ if(x==y) return Poly{1}; if(!leq(x,y)) return Poly{}; auto key=make_pair(x,y); auto it=memo.find(key); if(it!=memo.end()) return it->second;
    int s=__builtin_ctz(L[y]); int u=left[y][s]; Poly p; int c=(L[x]>>s&1); int sx=left[x][s]; if(sx==-2){ /* sx not in ideal => sx>y... x<=y and sx>x, sy<y gives sx<=y; so cannot happen */ fprintf(stderr,"ord: sx missing\n"); exit(3);}
    addto(p,P(sx,u),1-c,1); addto(p,P(x,u),c,1);
    for(int z=0;z<(int)el.size();z++){ if(len[z]>=len[u]) break; if(!(L[z]>>s&1)) continue; int gap=len[u]-len[z]; if(gap%2==0) continue; if(!leq(x,z)) continue; Poly pz=P(z,u); size_t top=(gap-1)/2; if(pz.size()==top+1&&pz[top]) addto(p,P(x,z),(len[y]-len[z])/2,-pz[top]); }
    trim(p); memo[key]=p; return p; }
};

static void setup(const string&type){ ADJ.clear(); if(type=="E9"){ N=9; ADJ.assign(9,{}); auto add=[&](int a,int b){ADJ[a].push_back(b);ADJ[b].push_back(a);}; for(int i=0;i<7;i++) add(i,i+1); add(2,8);} else if(type=="E6"){ N=6; ADJ.assign(6,{}); auto add=[&](int a,int b){ADJ[a].push_back(b);ADJ[b].push_back(a);}; for(int i=0;i<4;i++) add(i,i+1); add(2,5);} else if(type=="D4"){ N=4; ADJ.assign(4,{}); auto add=[&](int a,int b){ADJ[a].push_back(b);ADJ[b].push_back(a);}; add(0,1);add(1,2);add(1,3);} else if(type=="A4"){N=4; ADJ.assign(4,{}); auto add=[&](int a,int b){ADJ[a].push_back(b);ADJ[b].push_back(a);}; for(int i=0;i<3;i++) add(i,i+1);} else {fprintf(stderr,"unknown type\n");exit(2);} }

static long RANDSEED=-1;
static string maskstr(int m){ string s="{"; for(int i=0;i<N;i++) if(m>>i&1){ if(s.size()>1) s+=","; s+=to_string(i);} return s+"}"; }

int main(int argc,char**argv){
  if(argc<2){fprintf(stderr,"usage: klq validate <type> <seed> <count> <maxlen> | klq run <word letters...> --I a,b,c --x words | klq info ...\n");return 2;}
  string mode=argv[1];
  if(mode=="validate"){ setup(argv[2]); unsigned seed=atoi(argv[3]); int count=atoi(argv[4]); int maxlen=atoi(argv[5]); mt19937 rng(seed); int tested=0, pairs=0;
    for(int it=0;it<count;it++){ int L=1+rng()%maxlen; Mat w=identity(); for(int i=0;i<L;i++){ int s=rng()%N; Mat y=rmul(w,s); if(negcol(w,s)) continue; w=y; } int R=Rmask(w); if(!R) continue; // random subset I of R generating a finite parabolic: restrict to commuting subsets or small
      int I=0; for(int s=0;s<N;s++) if((R>>s&1)&&(rng()%2)) I|=1<<s; if(!I) I=1<<__builtin_ctz(R); // W_I finite? in E9 any proper subset is finite; the only infinite is all 9 letters
      if(I==(1<<N)-1) continue;
      Ord o; o.build(w); if(o.el.size()>6000) continue; Quot q; q.build(w,I,false); KLQ k(q); int yq=q.idx[ fromword( rword(fromword([&]{ vector<int> wd=rword(w); return wd;}())) ) ]; // w' index
      // w' = w wI
      Mat wp=w; for(int t:rword(q.wI)) wp=rmul(wp,t); yq=q.idx[wp]; k.compute(yq);
      // compare for all x in quotient
      int ow=o.idx[w]; for(size_t xi=0;xi<q.el.size();xi++){ Mat xm=q.el[xi]; for(int t:rword(q.wI)) xm=rmul(xm,t); auto it=o.idx.find(xm); Poly ref= (it==o.idx.end())?Poly{}:o.P(it->second,ow); const Poly*mine=k.get(xi,yq); if(ref!=*mine){ printf("MISMATCH w len %d I=%s x len %d: ord %s quot %s\n",o.len[ow],maskstr(I).c_str(),q.len[xi],pstr(ref).c_str(),pstr(*mine).c_str()); return 1;} pairs++; }
      tested++; }
    printf("validated %d elements, %d quotient pairs, all equal\n",tested,pairs); return 0; }
  if(mode=="run"){ setup("E9"); vector<int> beta, wordv, ext; int I=0; int i=2; string cur="";
    while(i<argc){ string a=argv[i]; if(a.rfind("--",0)==0){ cur=a; if(cur=="--I"){ i++; string Is=argv[i]; for(char c:Is) if(isdigit(c)) I|=1<<(c-'0'); } else if(cur=="--rand"){ i++; RANDSEED=atoi(argv[i]); } }
      else if(cur=="--beta") beta.push_back(atoi(a.c_str())); else if(cur=="--word") wordv.push_back(atoi(a.c_str())); else if(cur=="--ext") ext.push_back(atoi(a.c_str())); else {fprintf(stderr,"bad arg %s\n",a.c_str()); return 2;} i++; }
    Mat w; if(!beta.empty()){ if((int)beta.size()!=N){fprintf(stderr,"beta needs %d coords\n",N);return 2;} w=reflection(beta);} else w=fromword(wordv);
    for(int s:ext){ if(negcol(w,s)){fprintf(stderr,"extension letter is a descent\n");return 2;} w=rmul(w,s);}
    vector<int> wd=rword(w); printf("word:"); for(int s:wd) printf(" %d",s); printf("\nlength %zu\n",wd.size()); int R=Rmask(w), L=Lmask(w); printf("R=%s L=%s\n",maskstr(R).c_str(),maskstr(L).c_str()); if(!I) I=R; printf("I=%s\n",maskstr(I).c_str());
    Mat winv=inverse(w); printf("involution: %s\n", winv==w?"yes":"no");
    // terminal test: for s in R(w), t~s: ws t > ws
    bool term=true; for(auto v:{w,winv}) for(int s=0;s<N;s++) if(negcol(v,s)){ Mat vs=rmul(v,s); for(int t:ADJ[s]) if(negcol(vs,t)) term=false; } printf("terminal: %s\n",term?"yes":"no");
    Quot q; auto t0=chrono::steady_clock::now(); q.build(w,I,true); printf("quotient size %zu (|W_I| = %d, full ideal size %zu)\n",q.el.size(), (1<<__builtin_popcount(I)) , q.el.size()*(1<<__builtin_popcount(I))); fflush(stdout);
    if(!(__builtin_popcount(I)==(int)rword(q.wI).size())){ printf("I not commuting; full ideal size claim invalid\n"); }
    // eligible FC bottoms: x = x' wI with I subset L(x), FC
    vector<int> elig; for(size_t xi=0;xi<q.el.size();xi++){ Mat xm=q.el[xi]; for(int t:rword(q.wI)) xm=rmul(xm,t); if((Lmask(xm)&L)!=L) continue; if((Rmask(xm)&R)!=R) continue; if(!isFC(xm)) continue; elig.push_back(xi); }
    printf("eligible FC bottoms (descent-containing, FC, <= w): %zu\n",elig.size()); for(int xi:elig){ Mat xm=q.el[xi]; for(int t:rword(q.wI)) xm=rmul(xm,t); auto xw=rword(xm); printf("  x = "); for(int s:xw) printf("%d",s); printf("  length %zu gap %zu\n",xw.size(),wd.size()-xw.size()); }
    fflush(stdout);
    Mat wp=w; for(int t:rword(q.wI)) wp=rmul(wp,t); int yq=q.idx[wp]; KLQ k(q); if(RANDSEED>=0){k.randomize=true; k.rng.seed(RANDSEED);} k.compute(yq);
    auto t1=chrono::steady_clock::now(); printf("KL done: %zu quotient elements computed, %zu stored pairs, %zu distinct polys, %.1f s\n",k.computed,k.pairs,k.polys.size(),chrono::duration<double>(t1-t0).count());
    for(int xi:elig){ Mat xm=q.el[xi]; for(int t:rword(q.wI)) xm=rmul(xm,t); auto xw=rword(xm); const Poly*p=k.get(xi,yq); int gap=wd.size()-xw.size(); ll mu=0; if(gap%2==1){ size_t top=(gap-1)/2; if(p->size()==top+1) mu=(*p)[top]; } printf("P(x,w) for x="); for(int s:xw) printf("%d",s); printf(" : %s   gap %d  mu = %lld\n",pstr(*p).c_str(),gap,mu); }
    // also report full statistics: max mu over all extremal x in list of w'
    ll maxmu=0; int cnt_nonzero_mu=0; map<ll,int> muhist; for(auto&pr:k.lst[yq]){ int x=pr.first; if(x==yq) continue; int gap=q.len[yq]-q.len[x]; if(gap%2==0) continue; const Poly&p=k.polys[pr.second]; size_t top=(gap-1)/2; if(p.size()==top+1&&p[top]){ muhist[p[top]]++; maxmu=max(maxmu,p[top]); cnt_nonzero_mu++; } }
    printf("nonzero mu(x,w) over extremal x in quotient: %d; histogram:",cnt_nonzero_mu); for(auto&h:muhist) printf(" mu=%lld:%d",h.first,h.second); printf("\n");
    if(winv==w){ size_t checked=0, bad=0, missing=0; for(auto&pr:k.lst[yq]){ Mat X=q.el[pr.first]; for(int t:rword(q.wI)) X=rmul(X,t); Mat Xi=inverse(X); for(int t:rword(q.wI)) Xi=rmul(Xi,t); auto it=q.idx.find(Xi); if(it==q.idx.end()){missing++; continue;} const Poly*p2=k.get(it->second,yq); if(*p2!=k.polys[pr.second]) bad++; checked++; }
      printf("inversion symmetry check P(x,w)=P(x^-1,w): %zu pairs checked, %zu mismatches, %zu inverses missing from quotient\n",checked,bad,missing); }
    // max coefficient anywhere
    ll maxc=0; size_t maxdeg=0; for(auto&p:k.polys){ for(ll c:p) maxc=max(maxc,c); maxdeg=max(maxdeg,p.size()); } printf("max coefficient seen %lld, max degree %zu\n",maxc,maxdeg?maxdeg-1:0);
    return 0; }
  if(mode=="dump"){ setup(argv[2]); vector<int> wd; for(int i=3;i<argc;i++) wd.push_back(atoi(argv[i])); Mat w=fromword(wd); Ord o; o.build(w); int ow=o.idx[w];
    for(size_t x=0;x<o.el.size();x++){ Poly p=o.P(x,ow); if(p.empty()) continue; for(int j=0;j<N;j++){ for(int i=0;i<N;i++) printf("%d%c",(int)o.el[x][j*9+i],i+1<N?',':';'); } printf(" %s\n",pstr(p).c_str()); }
    return 0; }
  fprintf(stderr,"unknown mode\n"); return 2; }
