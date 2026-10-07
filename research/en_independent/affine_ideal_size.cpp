// Exact bounded principal-lower-ideal enumeration in affine E8, using subwords.
// Input: a reduced word (space-separated digits0..8) on stdin.
// Argument: maximum number of elements retained (default3000000).
// A cap exit is explicitly incomplete; no polynomial or cardinality inferred.
#include <array>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <vector>
using Mat=std::array<int8_t,81>;
struct Hash {size_t operator()(const Mat&a)const{size_t h=1469598103934665603ULL;for(auto c:a)h=(h^uint8_t(c))*1099511628211ULL;return h;}};
std::array<std::vector<int>,9> adj;
Mat right(const Mat&a,int s){Mat b=a;for(int k=0;k<9;k++)b[s*9+k]=-a[s*9+k];
  for(int t:adj[s])for(int k=0;k<9;k++){int x=int(a[s*9+k])+a[t*9+k];assert(x>=-127&&x<=127);b[t*9+k]=x;}return b;}
bool negative(const Mat&a,int s){bool p=false,m=false;for(int k=0;k<9;k++){p|=a[s*9+k]>0;m|=a[s*9+k]<0;}assert(p!=m);return m;}
int main(int argc,char**argv){uint64_t cap=argc>1?strtoull(argv[1],nullptr,10):3000000;
  for(int s=0;s<7;s++){adj[s].push_back(s+1);adj[s+1].push_back(s);}adj[2].push_back(8);adj[8].push_back(2);
  std::vector<int> word;int s;while(std::cin>>s){assert(s>=0&&s<9);word.push_back(s);}
  Mat e{};for(int j=0;j<9;j++)e[j*9+j]=1;
  std::vector<Mat> els{e};std::vector<uint8_t> len{0};std::unordered_map<Mat,uint32_t,Hash> ids;ids[e]=0;
  std::vector<uint64_t> sizes;Mat top=e;bool complete=true;
  for(size_t k=0;k<word.size();k++){
    int a=word[k];assert(!negative(top,a));top=right(top,a);size_t count=els.size();
    for(size_t i=0;i<count;i++){
      Mat b=right(els[i],a);if(ids.count(b))continue;
      assert(!negative(els[i],a));
      if(els.size()>=cap){complete=false;break;}
      ids.emplace(b,els.size());els.push_back(b);len.push_back(len[i]+1);
    }
    sizes.push_back(els.size());std::cerr<<"prefix "<<k+1<<" lower_size "<<els.size()<<'\n';
    if(!complete)break;
  }
  std::vector<uint64_t> hist(word.size()+1,0);for(auto l:len)hist[l]++;
  std::cout<<"{\"complete\":"<<(complete?"true":"false")<<",\"cap\":"<<cap<<",\"elements_retained\":"<<els.size()<<",\"word\":[";
  for(size_t k=0;k<word.size();k++)std::cout<<(k?",":"")<<word[k];
  std::cout<<"],\"prefix_lower_sizes\":[";for(size_t k=0;k<sizes.size();k++)std::cout<<(k?",":"")<<sizes[k];
  std::cout<<"],\"length_histogram\":[";for(size_t k=0;k<hist.size();k++)std::cout<<(k?",":"")<<hist[k];std::cout<<"]}\n";
  return complete?0:3;
}
