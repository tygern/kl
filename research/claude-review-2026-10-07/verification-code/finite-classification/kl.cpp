// From-scratch KL polynomial computation on the lower Bruhat interval of a given word in E_n.
// Geometric representation, exact 64-bit integer coefficients.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <map>
#include <unordered_map>
using namespace std;
static int N; static int C[8][8]; static bool ADJ[8][8];
struct Mat { int8_t a[8][8]; bool operator<(const Mat&o)const{return memcmp(a,o.a,64)<0;} bool operator==(const Mat&o)const{return memcmp(a,o.a,64)==0;} };
struct MatHash{ size_t operator()(const Mat&m)const{ uint64_t h=1469598103934665603ULL; const uint8_t*p=(const uint8_t*)m.a; for(int i=0;i<64;i++){h^=p[i];h*=1099511628211ULL;} return h;} };
static void setup(int n){N=n;memset(C,0,sizeof C);memset(ADJ,0,sizeof ADJ);for(int i=0;i<n;i++)C[i][i]=2;auto e=[&](int i,int j){ADJ[i][j]=ADJ[j][i]=1;C[i][j]=C[j][i]=-1;};for(int i=0;i+1<=n-2;i++)e(i,i+1);e(2,n-1);}
static Mat identity(){Mat m;memset(m.a,0,64);for(int i=0;i<N;i++)m.a[i][i]=1;return m;}
static void rightMul(Mat&m,int i){int8_t o[8];memcpy(o,m.a[i],8);for(int k=0;k<N;k++)m.a[i][k]=-o[k];for(int j=0;j<N;j++)if(ADJ[i][j])for(int k=0;k<N;k++)m.a[j][k]+=o[k];}
static void leftMul(Mat&m,int i){for(int j=0;j<N;j++){int p=0;for(int k=0;k<N;k++)p+=m.a[j][k]*C[k][i];m.a[j][i]-=p;}}
static bool neg(const int8_t*v){bool nz=false;for(int k=0;k<N;k++){if(v[k]>0)return false;if(v[k]<0)nz=true;}return nz;}
static bool rightDescent(const Mat&m,int i){return neg(m.a[i]);}
// left descent: s in L(w) iff w^{-1}(alpha_s)<0 iff (alpha_s, w(alpha_j))-pattern... use: s in L(w) iff ell(sw)<ell(w); compute via inversion count of columns? simpler: w^{-1}alpha_s <0 iff the root alpha_s is a left inversion; test: s in L(w) iff s w has smaller length. Use length via positive roots count below.
static int length(Mat m){int l=0;while(true){int s=-1;for(int i=0;i<N;i++)if(rightDescent(m,i)){s=i;break;}if(s<0)break;rightMul(m,s);l++;}return l;}
static bool leftDescent(const Mat&m,int s){ // s in L(w) iff w^{-1}(alpha_s) < 0 iff B(w(alpha_j), alpha_s) pattern: w^{-1}alpha_s = sum_j B(alpha_s, w alpha_j) omega_j^vee... use length instead (cheap enough for small interval)
    Mat t=m; leftMul(t,s); return length(t)<length(m); }
typedef vector<long long> Poly;
static void addShift(Poly&a,const Poly&b,int sh,long long c){ if(a.size()<b.size()+sh)a.resize(b.size()+sh,0); for(size_t i=0;i<b.size();i++)a[i+sh]+=c*b[i]; }
static void trim(Poly&a){while(!a.empty()&&a.back()==0)a.pop_back();}
int main(int argc,char**argv){
    if(argc<4){fprintf(stderr,"usage: kl <n> <word w> <word x>\n");return 1;}
    setup(atoi(argv[1])); string W=argv[2], X=argv[3];
    // lower interval via subwords
    vector<Mat> elems; { unordered_map<Mat,int,MatHash> seen; int L=W.size();
      for(long long mask=0;mask<(1LL<<L);mask++){ Mat m=identity(); for(int i=0;i<L;i++) if(mask>>i&1) rightMul(m,W[i]-'0'); if(!seen.count(m)){seen[m]=elems.size();elems.push_back(m);} } }
    int M=elems.size(); vector<int> len(M); for(int i=0;i<M;i++)len[i]=length(elems[i]);
    vector<int> order(M); for(int i=0;i<M;i++)order[i]=i; sort(order.begin(),order.end(),[&](int a,int b){return len[a]<len[b];});
    vector<Mat> E(M); vector<int> LEN(M); for(int i=0;i<M;i++){E[i]=elems[order[i]];LEN[i]=len[order[i]];}
    unordered_map<Mat,int,MatHash> idx; for(int i=0;i<M;i++)idx[E[i]]=i;
    printf("interval size %d, top length %d\n",M,LEN[M-1]);
    // left descents and s-multiples
    vector<array<int,8>> leftMulIdx(M); vector<array<bool,8>> LD(M);
    for(int i=0;i<M;i++) for(int s=0;s<N;s++){ Mat t=E[i]; leftMul(t,s); auto it=idx.find(t); leftMulIdx[i][s]= it==idx.end()?-1:it->second; LD[i][s]= (it!=idx.end() && LEN[it->second]<LEN[i]); }
    // P[y] : sparse map x -> Poly (only x<=y)
    vector<unordered_map<int,Poly>> P(M);
    P[0][0]=Poly{1};
    for(int y=1;y<M;y++){
        int s=-1; for(int t=0;t<N;t++) if(LD[y][t]){s=t;break;}
        int v=leftMulIdx[y][s];
        // z candidates: z<=v with s in L(z), mu(z,v)!=0
        vector<pair<int,long long>> zs;
        for(auto&kv:P[v]){ int z=kv.first; if(!LD[z][s]) continue; int d=LEN[v]-LEN[z]; if(d%2==0) continue; int k=(d-1)/2; const Poly&p=kv.second; if((int)p.size()>k && p[k]!=0) zs.push_back({z,p[k]}); }
        for(int x=0;x<M;x++){ if(LEN[x]>LEN[y]) break;
            int sx=leftMulIdx[x][s]; int c = LD[x][s]?1:0; Poly r;
            if(sx>=0){ auto it=P[v].find(sx); if(it!=P[v].end()) addShift(r,it->second,1-c,1); }
            { auto it=P[v].find(x); if(it!=P[v].end()) addShift(r,it->second,c,1); }
            for(auto&zm:zs){ auto it=P[zm.first].find(x); if(it!=P[zm.first].end()) addShift(r,it->second,(LEN[y]-LEN[zm.first])/2,-zm.second); }
            trim(r); if(!r.empty()){ if(r[0]!=1){fprintf(stderr,"constant term !=1 at y=%d x=%d\n",y,x);exit(1);} P[y][x]=r; }
        }
    }
    Mat w=identity(); for(char ch:W) rightMul(w,ch-'0'); Mat x=identity(); for(char ch:X) rightMul(x,ch-'0');
    int wi=idx[w], xi=idx.count(x)?idx[x]:-1;
    printf("|[e,w]| = %zu ; x in interval: %d\n", P[wi].size(), xi>=0);
    if(xi>=0){ Poly p=P[wi][xi]; printf("P_{x,w} ="); for(size_t i=0;i<p.size();i++) printf(" %lld q^%zu", p[i], i); int d=LEN[wi]-LEN[xi]; printf("\nlength gap %d, mu(x,w) = %lld\n", d, (d%2==1 && (int)p.size()>(d-1)/2)? p[(d-1)/2]:0); }
    // all x with mu(x,w)!=0 and both descent inclusions, FC or not
    int cnt=0; for(auto&kv:P[wi]){ int d=LEN[wi]-LEN[kv.first]; if(d%2==1 && d>1){ int k=(d-1)/2; if((int)kv.second.size()>k && kv.second[k]!=0) cnt++; } }
    printf("number of non-cover x<w with mu(x,w)!=0: %d\n",cnt);
    // sanity: P_{e,w}
    Poly pe=P[wi][0]; printf("P_{e,w} ="); for(size_t i=0;i<pe.size();i++) printf(" %lld q^%zu", pe[i], i); printf("\n");
}
