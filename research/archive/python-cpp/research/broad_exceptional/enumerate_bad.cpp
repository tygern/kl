// Exhaustive exact-root enumeration, E6/E7; E8 available but very expensive.
// Diagram 0--1--2--3--...--(n-2), node n-1 attached to 2.
#include <array>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <cstdint>
using namespace std;
using Root=array<int,8>;
int n; vector<int> adj[8];
vector<Root> roots;map<Root,int> ri; int negs[256],pos[256],adds[256][256];
int rget(uint64_t x,int s){return (x>>(8*s))&255;}
uint64_t rset(uint64_t x,int s,int v){return (x&~(uint64_t(255)<<(8*s)))|(uint64_t(v)<<(8*s));}
uint64_t mult(uint64_t w,int s){auto v=rset(w,s,negs[rget(w,s)]);for(int t:adj[s])v=rset(v,t,adds[rget(w,s)][rget(w,t)]);return v;}
vector<uint64_t> elems;unordered_map<uint64_t,uint32_t> indexer;
vector<array<uint32_t,8>> rights;vector<uint32_t> parent;vector<uint8_t> lastg,lens,ds,fc;
vector<int> word(uint32_t i){vector<int> w;while(i){w.push_back(lastg[i]);i=parent[i];}reverse(w.begin(),w.end());return w;}
uint32_t inv(uint32_t i){auto w=word(i);uint32_t v=0;for(auto s=w.rbegin();s!=w.rend();++s)v=rights[v][*s];return v;}
bool weakright(uint32_t i){for(int s=0;s<n;s++)if(ds[i]&(1<<s))for(int t:adj[s])if(ds[rights[i][s]]&(1<<t))return false;return true;}
bool leq(uint32_t x,uint32_t w){while(x!=w){if(lens[x]>=lens[w])return false;int s=__builtin_ctz(ds[w]);if(ds[x]&(1<<s))x=rights[x][s];w=rights[w][s];}return true;}
void printword(uint32_t i){auto w=word(i);cout<<'[';for(size_t j=0;j<w.size();j++)cout<<(j?",":"")<<w[j];cout<<']';}
int main(int argc,char**argv){n=argc>1?atoi(argv[1]):7;assert(n>=6&&n<=8);for(int s=0;s<n-2;s++){adj[s].push_back(s+1);adj[s+1].push_back(s);}adj[2].push_back(n-1);adj[n-1].push_back(2);
for(int s=0;s<n;s++){Root a{};a[s]=1;ri[a]=roots.size();roots.push_back(a);}for(size_t i=0;i<roots.size();i++)for(int s=0;s<n;s++){auto a=roots[i];a[s]=-a[s];for(int t:adj[s])a[s]+=roots[i][t];if(!ri.count(a)){ri[a]=roots.size();roots.push_back(a);}}
assert(roots.size()==(n==6?72:n==7?126:240));for(int i=0;i<(int)roots.size();i++){auto a=roots[i];pos[i]=true;for(int s=0;s<n;s++){pos[i]&=a[s]>=0;a[s]=-a[s];}negs[i]=ri.at(a);for(int j=0;j<(int)roots.size();j++){Root b{};for(int s=0;s<n;s++)b[s]=roots[i][s]+roots[j][s];auto it=ri.find(b);adds[i][j]=(it==ri.end()?-1:it->second);}}
size_t expected=n==6?51840:n==7?2903040:696729600;indexer.reserve(min(expected,size_t(3000000)));elems.reserve(min(expected,size_t(3000000)));rights.reserve(min(expected,size_t(3000000)));uint64_t e=0;for(int s=0;s<n;s++)e=rset(e,s,s);elems.push_back(e);indexer[e]=0;parent.push_back(0);lastg.push_back(0);lens.push_back(0);
for(size_t i=0;i<elems.size();i++){array<uint32_t,8> row{};uint8_t d=0;for(int s=0;s<n;s++){auto v=mult(elems[i],s);auto it=indexer.find(v);uint32_t k;if(it==indexer.end()){k=elems.size();indexer[v]=k;elems.push_back(v);parent.push_back(i);lastg.push_back(s);lens.push_back(lens[i]+1);}else k=it->second;row[s]=k;if(!pos[rget(elems[i],s)])d|=(1<<s);}rights.push_back(row);ds.push_back(d);bool f=true;for(int s=0;s<n;s++)if(d&(1<<s)){f&=fc[row[s]];for(int t:adj[s])f&=!(d&(1<<t));}fc.push_back(f);}
assert(elems.size()==expected);
size_t commuting_weak=0,weakfc_noncommuting=0,star_checks=0; vector<size_t> hist(lens.back()+1,0);
for(size_t i=0;i<elems.size();i++){
 hist[lens[i]]++;
 for(int s=0;s<n;s++){assert(abs(int(lens[rights[i][s]])-int(lens[i]))==1);assert(bool(ds[i]&(1<<s))==(lens[rights[i][s]]<lens[i]));assert(rights[rights[i][s]][s]==i);}
 if(fc[i]){
  if(weakright(i)&&weakright(inv(i))){auto wd=word(i);bool comm=true;for(int s:wd)for(int t:adj[s])if(find(wd.begin(),wd.end(),t)!=wd.end())comm=false;comm&=set<int>(wd.begin(),wd.end()).size()==wd.size();if(comm)commuting_weak++;else weakfc_noncommuting++;}
  for(int s=0;s<n;s++)for(int t:adj[s])if(s<t && bool(ds[i]&(1<<s))!=bool(ds[i]&(1<<t))){auto a=rights[i][s],b=rights[i][t];bool da=bool(ds[a]&(1<<s))!=bool(ds[a]&(1<<t)),db=bool(ds[b]&(1<<s))!=bool(ds[b]&(1<<t));assert(da!=db);assert(fc[da?a:b]);star_checks++;}
 }
}
assert(!weakfc_noncommuting);
vector<uint32_t> fcs,bads;vector<uint8_t> lds(elems.size());for(size_t i=0;i<elems.size();i++){if(fc[i]){fcs.push_back(i);lds[i]=ds[inv(i)];}if(!fc[i]&&weakright(i)&&weakright(inv(i)))bads.push_back(i);}
cout<<"{\"rank\":"<<n<<",\"roots\":"<<roots.size()<<",\"order\":"<<elems.size()<<",\"fc_count\":"<<fcs.size()<<",\"max_length\":"<<int(lens.back())<<",\"commuting_weak\":"<<commuting_weak<<",\"fc_right_star_checks\":"<<star_checks<<",\"length_distribution\":[";for(size_t k=0;k<hist.size();k++)cout<<(k?",":"")<<hist[k];cout<<"],\"bad\":[\n";bool first=true;
for(auto w:bads){if(!first)cout<<",\n";first=false;auto ld=ds[inv(w)];cout<<"{\"id\":"<<w<<",\"word\":";printword(w);cout<<",\"length\":"<<int(lens[w])<<",\"Rmask\":"<<int(ds[w])<<",\"Lmask\":"<<int(ld)<<",\"bottoms\":[";bool ff=true;for(auto x:fcs)if((ds[x]&ds[w])==ds[w]&&(lds[x]&ld)==ld&&leq(x,w)){if(!ff)cout<<',';ff=false;cout<<"{\"id\":"<<x<<",\"word\":";printword(x);cout<<",\"length\":"<<int(lens[x])<<",\"rank\":"<<int(lens[w]-lens[x])<<'}';}cout<<"]}";}
cout<<"\n]}\n";}
