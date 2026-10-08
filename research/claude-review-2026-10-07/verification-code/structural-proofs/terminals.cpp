// Independent enumeration of terminal elements in E6, E7 (full BFS) and E8 (parabolic pruning).
// Elements are stored as integer matrices in the simple-root basis: column j = w(alpha_j).
// Numbering (manuscript): chain 0-1-...-(n-2), node n-1 attached to node 2.
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <array>
#include <vector>
#include <set>
#include <string>
#include <unordered_set>
#include <algorithm>
using namespace std;
typedef array<int8_t,64> Mat; // n<=8, entry (i,j) at i*8+j : coefficient of alpha_i in w(alpha_j)
int n;
vector<vector<int>> adj;
struct H { size_t operator()(const Mat&m) const { size_t h=1469598103934665603ULL; for(int k=0;k<64;k++){h^=(uint8_t)m[k]; h*=1099511628211ULL;} return h; } };
Mat ident(){ Mat m{}; for(int i=0;i<n;i++) m[i*8+i]=1; return m; }
// right multiply by s_i: a_i -> -a_i, a_j -> a_j + a_i for j~i
Mat rmul(const Mat& w,int i){ Mat r=w; for(int k=0;k<n;k++){ r[k*8+i]=-w[k*8+i]; } for(int j:adj[i]) for(int k=0;k<n;k++) r[k*8+j]=w[k*8+j]+w[k*8+i]; return r; }
// left multiply by s_i: apply s_i to each column: v -> v - (v,alpha_i) alpha_i, (v,alpha_i)=2v_i - sum_{j~i} v_j
Mat lmul(int i,const Mat& w){ Mat r=w; for(int c=0;c<n;c++){ int p=2*w[i*8+c]; for(int j:adj[i]) p-=w[j*8+c]; r[i*8+c]=w[i*8+c]-p; } return r; }
Mat mul(const Mat&a,const Mat&b){ Mat r{}; for(int i=0;i<n;i++)for(int j=0;j<n;j++){int s=0; for(int k=0;k<n;k++) s+=a[i*8+k]*b[k*8+j]; r[i*8+j]=s;} return r; }
// sign of column j (positive root => +1)
int colsign(const Mat&w,int j){ for(int i=0;i<n;i++){ if(w[i*8+j]>0) return 1; if(w[i*8+j]<0) return -1;} return 0; }
bool rdesc(const Mat&w,int j){ return colsign(w,j)<0; }
// right-terminal test, eq (terminaltest): s in R(w) => w(alpha_s+alpha_t) > 0 for all t~s
bool rightTerminal(const Mat&w){ for(int s=0;s<n;s++) if(rdesc(w,s)) for(int t:adj[s]){ int sg=0; for(int i=0;i<n;i++){int v=w[i*8+s]+w[i*8+t]; if(v>0){sg=1;break;} if(v<0){sg=-1;break;}} if(sg<=0) return false; } return true; }
// independent brute check of right-terminality: exists s in R(w), t~s with t in R(ws)?
bool rightTerminalBrute(const Mat&w){ for(int s=0;s<n;s++) if(rdesc(w,s)){ Mat ws=rmul(w,s); for(int t:adj[s]) if(rdesc(ws,t)) return false;} return true; }
Mat inverse(const Mat&w){ // via reduced word: strip right descents
  vector<int> word; Mat cur=w; while(true){ int s=-1; for(int j=0;j<n;j++) if(rdesc(cur,j)){s=j;break;} if(s<0)break; cur=rmul(cur,s); word.push_back(s);} // w = ... ; cur*word reversed... w = e * s_{k} ... careful
  // cur is identity now; w = s_{word[k-1]} ... s_{word[0]}  (since we removed word[0] last from the right)
  Mat inv=ident(); for(int k=0;k<(int)word.size();k++) inv=rmul(inv,word[k]); // inv = s_{word[0]} s_{word[1]} ... = w^{-1}
  return inv; }
vector<int> redword(const Mat&w){ vector<int> word; Mat cur=w; while(true){ int s=-1; for(int j=0;j<n;j++) if(rdesc(cur,j)){s=j;break;} if(s<0)break; cur=rmul(cur,s); word.push_back(s);} reverse(word.begin(),word.end()); return word; }
int length(const Mat&w){ return redword(w).size(); }
bool isFC(const Mat&w){ // check: no reduced word contains sts; brute: BFS over reduced words via commutation class? simpler: element FC iff for every prefix u (w = u v) and s,t adjacent, not both s,t in L(v)... use: w FC iff no x<=w ... use Stembridge: w is FC iff for all reduced words no sts factor. Compute all reduced words set of heaps is costly; instead use: w not FC iff exists s,t adjacent and u with w = u*sts*v reduced, i.e. exists prefix u of w such that sts is a reduced prefix of u^{-1}w ... iterate over all prefixes via BFS on the set of right factors.
  // BFS over right factors v of w (all elements v with w = u v length-additive) - these are the elements v with v <= w in weak order sense from the right. Count small.
  set<vector<int8_t>> seen; vector<Mat> st; st.push_back(w);
  auto key=[&](const Mat&m){ return vector<int8_t>(m.begin(),m.end()); };
  seen.insert(key(w));
  while(!st.empty()){ Mat v=st.back(); st.pop_back(); // check if v has two adjacent left descents
    Mat vi=inverse(v); for(int s=0;s<n;s++) if(rdesc(vi,s)) for(int t:adj[s]) if(rdesc(vi,t)) return false;
    for(int s=0;s<n;s++) if(rdesc(vi,s)){ Mat nv=lmul(s,v); auto k=key(nv); if(!seen.count(k)){seen.insert(k); st.push_back(nv);} } }
  return true; }
string wstr(const vector<int>&w){ string s; for(int x:w) s+=char('0'+x); return s; }
void setType(int N){ n=N; adj.assign(n,{}); for(int i=0;i+1<=n-2;i++){ adj[i].push_back(i+1); adj[i+1].push_back(i);} adj[n-1].push_back(2); adj[2].push_back(n-1); }
set<int> descR(const Mat&w){ set<int> r; for(int j=0;j<n;j++) if(rdesc(w,j)) r.insert(j); return r; }
string sstr(const set<int>&s){ string r; for(int x:s) r+=char('0'+x); return r; }
bool isCommProduct(const Mat&w){ auto word=redword(w); set<int> s(word.begin(),word.end()); if(s.size()!=word.size()) return false; for(int a:s) for(int b:adj[a]) if(s.count(b)) return false; return true; }
void report(const Mat&w){ Mat inv=inverse(w); auto wd=redword(w); printf("  terminal len=%d word=%s L=%s R=%s comm=%d  rightTermBrute=%d leftTermBrute=%d\n",(int)wd.size(),wstr(wd).c_str(),sstr(descR(inv)).c_str(),sstr(descR(w)).c_str(),(int)isCommProduct(w),(int)rightTerminalBrute(w),(int)rightTerminalBrute(inv)); }
// Full BFS enumeration restricted to generators in gens; returns all elements
vector<Mat> enumerate(const vector<int>&gens){ unordered_set<Mat,H> seen; vector<Mat> all; vector<Mat> frontier{ident()}; seen.insert(ident()); all.push_back(ident());
  while(!frontier.empty()){ vector<Mat> nxt; for(auto&w:frontier) for(int s:gens){ Mat v=rmul(w,s); if(!seen.count(v)){ seen.insert(v); all.push_back(v); nxt.push_back(v);} } frontier.swap(nxt);} return all; }
int main(int argc,char**argv){
  int N=atoi(argv[1]); setType(N);
  if(N<=7){ vector<int> gens; for(int i=0;i<n;i++) gens.push_back(i); auto all=enumerate(gens); printf("E%d order %zu\n",N,all.size());
    int nright=0, nterm=0, ncomm=0; vector<Mat> terms;
    for(auto&w:all){ if(rightTerminal(w)!=rightTerminalBrute(w)){printf("MISMATCH test\n");return 1;} if(rightTerminal(w)){ nright++; Mat inv=inverse(w); if(rightTerminal(inv)){ nterm++; terms.push_back(w); if(isCommProduct(w)) ncomm++; } } }
    printf("right-terminals %d, two-sided terminals %d, commuting products %d, noncommuting %d\n",nright,nterm,ncomm,nterm-ncomm);
    for(auto&w:terms) if(!isCommProduct(w)) report(w);
    return 0; }
  // E8: pruning with parabolic J given by argv[2] (string of nodes), e.g. "0123456" (E7) or "1234567" (D7)
  string J=argv[2]; vector<int> gens; for(char c:J) gens.push_back(c-'0'); vector<bool> inJ(n,false); for(int g:gens) inJ[g]=true;
  auto WJ=enumerate(gens); printf("E8 pruning with J=%s, |W_J|=%zu\n",J.c_str(),WJ.size());
  vector<Mat> vt; for(auto&v:WJ) if(rightTerminal(v)) vt.push_back(v); printf("right-terminals in W_J: %zu\n",vt.size());
  // coset reps: closure under left multiplication, reduce by stripping right J-descents
  auto reduce=[&](Mat a){ while(true){ int s=-1; for(int g:gens) if(rdesc(a,g)){s=g;break;} if(s<0)break; a=rmul(a,s);} return a; };
  unordered_set<Mat,H> seen; vector<Mat> reps{ident()}; seen.insert(ident()); for(size_t i=0;i<reps.size();i++){ for(int s=0;s<n;s++){ Mat b=reduce(lmul(s,reps[i])); if(!seen.count(b)){seen.insert(b); reps.push_back(b);} } }
  printf("coset reps: %zu\n",reps.size());
  // cross-check: number of reps should equal |W|/|W_J| ; also verify each rep has no right J-descent and reps distinct
  long long order=696729600LL; printf("expected cosets %lld\n", order/(long long)WJ.size());
  vector<Mat> rt; for(auto&a:reps) for(auto&v:vt){ Mat w=mul(a,v); if(rightTerminal(w)) rt.push_back(w);} printf("right-terminals in E8: %zu\n",rt.size());
  int nterm=0,ncomm=0; vector<Mat> terms; for(auto&w:rt){ Mat inv=inverse(w); if(rightTerminal(inv)){ nterm++; terms.push_back(w); if(isCommProduct(w)) ncomm++; } }
  printf("two-sided terminals %d, commuting %d, noncommuting %d\n",nterm,ncomm,nterm-ncomm);
  for(auto&w:terms) if(!isCommProduct(w)) report(w);
  // also length check for products: verify length additivity l(av)=l(a)+l(v) for a sample
  return 0; }
