#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <array>
#include <vector>
#include <unordered_map>
#include <algorithm>
using namespace std;
typedef array<int8_t,64> Mat; int n; vector<vector<int>> adj;
struct H { size_t operator()(const Mat&m) const { size_t h=1469598103934665603ULL; for(int k=0;k<64;k++){h^=(uint8_t)m[k]; h*=1099511628211ULL;} return h; } };
Mat ident(){ Mat m{}; for(int i=0;i<n;i++) m[i*8+i]=1; return m; }
Mat rmul(const Mat& w,int i){ Mat r=w; for(int k=0;k<n;k++){ r[k*8+i]=-w[k*8+i]; } for(int j:adj[i]) for(int k=0;k<n;k++) r[k*8+j]=w[k*8+j]+w[k*8+i]; return r; }
int colsign(const Mat&w,int j){ for(int i=0;i<n;i++){ if(w[i*8+j]>0) return 1; if(w[i*8+j]<0) return -1;} return 0; }
bool rdesc(const Mat&w,int j){ return colsign(w,j)<0; }
void setType(int N){ n=N; adj.assign(n,{}); for(int i=0;i+1<=n-2;i++){ adj[i].push_back(i+1); adj[i+1].push_back(i);} adj[n-1].push_back(2); adj[2].push_back(n-1); }
int main(int argc,char**argv){ setType(atoi(argv[1]));
  // FC elements form a lower order ideal in right weak order. BFS: from FC w, for s not in R(w), v=ws; v is FC iff R(v) commuting and v t is FC for all t in R(v) (eq. FCrecurrence of the manuscript).
  unordered_map<Mat,int,H> fcidx; vector<Mat> fc; vector<int> len; fc.push_back(ident()); fcidx[ident()]=0; len.push_back(0);
  for(size_t i=0;i<fc.size();i++){ for(int s=0;s<n;s++){ if(rdesc(fc[i],s)) continue; Mat v=rmul(fc[i],s); if(fcidx.count(v)) continue; bool ok=true; vector<int> R; for(int j=0;j<n;j++) if(rdesc(v,j)) R.push_back(j);
      for(int a:R) for(int b:adj[a]) if(rdesc(v,b)) ok=false; if(ok) for(int t:R){ if(!fcidx.count(rmul(v,t))) { ok=false; break; } } // all shorter FC elements of length len+... already present? BFS order by length ensures elements of length len(v)-1 all processed before? Not necessarily; handle by processing in length layers.
      if(ok){ fcidx[v]=fc.size(); fc.push_back(v); len.push_back(len[i]+1);} } }
  int mx=0; for(int l:len) mx=max(mx,l); printf("E%d FC count %zu, max FC length %d\n",n,fc.size(),mx); }
