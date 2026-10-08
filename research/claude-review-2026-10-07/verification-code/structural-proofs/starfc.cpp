#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <array>
#include <vector>
#include <set>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <unordered_map>
using namespace std;
typedef array<int8_t,64> Mat; int n; vector<vector<int>> adj;
struct H { size_t operator()(const Mat&m) const { size_t h=1469598103934665603ULL; for(int k=0;k<64;k++){h^=(uint8_t)m[k]; h*=1099511628211ULL;} return h; } };
Mat ident(){ Mat m{}; for(int i=0;i<n;i++) m[i*8+i]=1; return m; }
Mat rmul(const Mat& w,int i){ Mat r=w; for(int k=0;k<n;k++){ r[k*8+i]=-w[k*8+i]; } for(int j:adj[i]) for(int k=0;k<n;k++) r[k*8+j]=w[k*8+j]+w[k*8+i]; return r; }
Mat lmul(int i,const Mat& w){ Mat r=w; for(int c=0;c<n;c++){ int p=2*w[i*8+c]; for(int j:adj[i]) p-=w[j*8+c]; r[i*8+c]=w[i*8+c]-p; } return r; }
int colsign(const Mat&w,int j){ for(int i=0;i<n;i++){ if(w[i*8+j]>0) return 1; if(w[i*8+j]<0) return -1;} return 0; }
bool rdesc(const Mat&w,int j){ return colsign(w,j)<0; }
void setType(int N){ n=N; adj.assign(n,{}); for(int i=0;i+1<=n-2;i++){ adj[i].push_back(i+1); adj[i+1].push_back(i);} adj[n-1].push_back(2); adj[2].push_back(n-1); }
int main(int argc,char**argv){ setType(atoi(argv[1]));
  // BFS enumerating group with lengths and FC flag: w is FC iff for every generator s in R(w), ws is FC, and R(w) commutative, and ... Use Stembridge: w FC iff w has no reduced word containing sts. Equivalent recursive: w is FC iff (R(w) is commutative) and (ws is FC for all s in R(w))? Not exactly (braid could be in the middle). Correct: w not FC iff exists reduced w = u sts v. Equivalently exists right factor y (w = u y reduced) with s,t in L(y) adjacent. So define bad(w) = exists s in R(w): bad(ws), or L(w) non-commutative. Then FC(w) = !bad(w). Need L(w): compute via inverse map. We store inverse index via BFS: inv of ws is s*inv(w).
  unordered_map<Mat,int,H> idx; vector<Mat> el, inv; vector<int> len; vector<char> fc;
  el.push_back(ident()); inv.push_back(ident()); len.push_back(0); idx[ident()]=0; fc.push_back(1);
  for(size_t i=0;i<el.size();i++){ for(int s=0;s<n;s++){ if(rdesc(el[i],s)) continue; Mat v=rmul(el[i],s); if(!idx.count(v)){ idx[v]=el.size(); el.push_back(v); inv.push_back(lmul(s,inv[i])); len.push_back(len[i]+1); fc.push_back(0);} } }
  // compute FC by increasing length: elements appended in BFS order = nondecreasing length
  for(size_t i=1;i<el.size();i++){ bool bad=false; // L noncommutative?
    for(int s=0;s<n&&!bad;s++) if(rdesc(inv[i],s)) for(int t:adj[s]) if(rdesc(inv[i],t)) {bad=true;break;}
    for(int s=0;s<n&&!bad;s++) if(rdesc(inv[i],s)){ int j=idx[lmul(s,el[i])]; if(!fc[j]) bad=true; }
    fc[i]=!bad; }
  long nfc=0; for(auto c:fc) nfc+=c; printf("order %zu, FC count %ld\n",el.size(),nfc);
  // star test: for all FC x, all adjacent (s,t) with x in D_L(s,t): *x = unique of sx,tx in D_L; check FC. Also right.
  long tested=0,fail=0,longer=0;
  for(size_t i=0;i<el.size();i++){ if(!fc[i]) continue; for(int s=0;s<n;s++) for(int t:adj[s]){ bool ls=rdesc(inv[i],s), lt=rdesc(inv[i],t); if(ls==lt) continue; // x in D_L(s,t)
      Mat a=lmul(s,el[i]), b=lmul(t,el[i]); int ia=idx[a], ib=idx[b]; auto inDL=[&](int k){ bool p=rdesc(inv[k],s), q=rdesc(inv[k],t); return p!=q; };
      int star=-1; if(inDL(ia)) star=ia; if(inDL(ib)){ if(star>=0) printf("two candidates!\n"); star=ib; } if(star<0){printf("none\n"); continue;}
      tested++; if(len[star]>len[i]) longer++; if(!fc[star]) { fail++; if(fail<5) printf("FAIL: x idx %zu len %d, star len %d\n",i,len[i],len[star]); } } }
  printf("left star tests %ld (longer %ld), FC failures %ld\n",tested,longer,fail); return 0; }
