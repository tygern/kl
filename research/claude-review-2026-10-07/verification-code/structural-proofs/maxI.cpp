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
bool indep(int mask){ for(int a=0;a<n;a++) if(mask>>a&1) for(int b:adj[a]) if(mask>>b&1) return false; return true; }
int main(int argc,char**argv){ setType(atoi(argv[1]));
  unordered_map<Mat,int,H> idx; vector<Mat> el, inv; vector<int> len, supp; vector<char> fc;
  el.push_back(ident()); inv.push_back(ident()); len.push_back(0); idx[ident()]=0; fc.push_back(1); supp.push_back(0);
  for(size_t i=0;i<el.size();i++){ for(int s=0;s<n;s++){ if(rdesc(el[i],s)) continue; Mat v=rmul(el[i],s); if(!idx.count(v)){ idx[v]=el.size(); el.push_back(v); inv.push_back(lmul(s,inv[i])); len.push_back(len[i]+1); fc.push_back(0); supp.push_back(supp[i]|(1<<s));} } }
  for(size_t i=1;i<el.size();i++){ bool bad=false; for(int s=0;s<n&&!bad;s++) if(rdesc(inv[i],s)) for(int t:adj[s]) if(rdesc(inv[i],t)) {bad=true;break;}
    for(int s=0;s<n&&!bad;s++) if(rdesc(inv[i],s)){ int j=idx[lmul(s,el[i])]; if(!fc[j]) bad=true; } fc[i]=!bad; }
  // alpha(K) for all K
  vector<int> alpha(1<<n,0); for(int K=0;K<(1<<n);K++){ int best=0; for(int I=K;;I=(I-1)&K){ if(indep(I)) best=max(best,__builtin_popcount(I)); if(I==0)break;} alpha[K]=best; }
  // commuting product element for independent I
  auto iI=[&](int I){ Mat m=ident(); for(int s=0;s<n;s++) if(I>>s&1) m=rmul(m,s); return m; };
  long checks=0, fails=0;
  for(size_t i=0;i<el.size();i++){ if(!fc[i]) continue; int Lm=0,Rm=0; for(int s=0;s<n;s++){ if(rdesc(inv[i],s)) Lm|=1<<s; if(rdesc(el[i],s)) Rm|=1<<s; } int LR=Lm&Rm;
    for(int K=0;K<(1<<n);K++){ if((supp[i]&~K)) continue; // K contains supp
      for(int I=K;;I=(I-1)&K){ if(indep(I) && __builtin_popcount(I)==alpha[K] && (I&~LR)==0){ checks++; if(!(el[i]==iI(I))){ fails++; if(fails<5) printf("FAIL x idx %zu len %d K=%x I=%x\n",i,len[i],K,I);} } if(I==0)break; } } }
  printf("E%d: FC %ld, hypothesis instances %ld, failures %ld\n",n,(long)count(fc.begin(),fc.end(),1),checks,fails); return 0; }
