// Direct enumeration of all fully commutative elements of generalized E_n.
// Unlike group BFS, only FC states are retained. Exact integer root matrices.
// Build: clang++ -O3 -std=c++17 fc_catalogue.cpp -o /private/tmp/en-fc
// Run: /private/tmp/en-fc 8 > e8-fc.json
#include <array>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <vector>
template<int N> struct Engine {
  using Mat=std::array<int32_t,N*N>;
  struct Hash { size_t operator()(const Mat&a) const {
    size_t h=1469598103934665603ULL;
    for(auto v:a) h=(h^uint32_t(v))*1099511628211ULL;
    return h;
  }};
  std::array<std::vector<int>,N> adj;
  std::vector<Mat> els;
  std::vector<uint32_t> parents, lengths;
  std::vector<uint8_t> last;
  std::unordered_map<Mat,uint32_t,Hash> ids;
  Engine(){
    for(int s=0;s<N-2;s++){adj[s].push_back(s+1);adj[s+1].push_back(s);}
    adj[2].push_back(N-1);adj[N-1].push_back(2);
  }
  Mat right(const Mat&a,int s) const {
    Mat b=a;
    for(int k=0;k<N;k++)b[s*N+k]=-a[s*N+k];
    for(int t:adj[s])for(int k=0;k<N;k++){
      int64_t c=int64_t(a[t*N+k])+a[s*N+k];
      assert(c>=INT32_MIN && c<=INT32_MAX);b[t*N+k]=c;
    }
    return b;
  }
  unsigned desc(const Mat&a) const {
    unsigned d=0;
    for(int s=0;s<N;s++){
      bool p=false,m=false;
      for(int k=0;k<N;k++){p|=a[s*N+k]>0;m|=a[s*N+k]<0;}
      assert(p!=m);if(m)d|=1u<<s;
    }
    return d;
  }
  bool independent(unsigned d) const {
    for(int s=0;s<N;s++)if(d&(1u<<s))for(int t:adj[s])if(d&(1u<<t))return false;
    return true;
  }
  std::vector<int> word(uint32_t i) const {
    std::vector<int>w;while(i){w.push_back(last[i]);i=parents[i];}
    std::reverse(w.begin(),w.end());return w;
  }
  void run(){
    Mat e{};for(int s=0;s<N;s++)e[s*N+s]=1;
    els.push_back(e);ids[e]=0;parents.push_back(0);lengths.push_back(0);last.push_back(0);
    uint64_t attempted=0;
    for(uint32_t i=0;i<els.size();i++){
      // Copy: push_back can reallocate the underlying matrix vector.
      const Mat a=els[i];unsigned d=desc(a);
      for(int s=0;s<N;s++)if(!(d&(1u<<s))){
        attempted++;Mat b=right(a,s);if(ids.count(b))continue;
        unsigned bd=desc(b);if(!independent(bd))continue;
        bool good=true;
        for(int t=0;t<N;t++)if(bd&(1u<<t)){
          auto it=ids.find(right(b,t));
          if(it==ids.end() || lengths[it->second]!=lengths[i]){good=false;break;}
        }
        if(!good)continue;
        uint32_t j=els.size();ids.emplace(b,j);els.push_back(b);
        parents.push_back(i);last.push_back(s);lengths.push_back(lengths[i]+1);
      }
      if(i && i%100000==0)std::cerr<<"E"<<N<<" processed "<<i<<" discovered "<<els.size()<<'\n';
    }
    std::vector<uint64_t> hist(lengths.back()+1,0);int maxcoord=0;
    for(uint32_t i=0;i<els.size();i++){
      hist[lengths[i]]++;assert(independent(desc(els[i])));
      for(auto c:els[i])maxcoord=std::max(maxcoord,int(std::abs(c)));
    }
    std::cout<<"{\"rank\":"<<N<<",\"fc_count\":"<<els.size()
      <<",\"max_fc_length\":"<<lengths.back()<<",\"max_root_coordinate\":"<<maxcoord
      <<",\"ascents_tested\":"<<attempted<<",\"length_distribution\":[";
    for(size_t k=0;k<hist.size();k++)std::cout<<(k?",":"")<<hist[k];
    std::cout<<"],\"elements\":[\n";
    for(uint32_t i=0;i<els.size();i++){
      auto w=word(i);Mat inverse=e;for(auto it=w.rbegin();it!=w.rend();++it)inverse=right(inverse,*it);
      std::cout<<(i?",\n":"")<<"{\"word\":[";
      for(size_t k=0;k<w.size();k++)std::cout<<(k?",":"")<<w[k];
      std::cout<<"],\"length\":"<<lengths[i]<<",\"Rmask\":"<<desc(els[i])<<",\"Lmask\":"<<desc(inverse)<<'}';
    }
    std::cout<<"\n],\"complete\":true,\"criterion\":\"commuting descents and all descent predecessors FC; BFS closed under FC ascents\"}\n";
  }
};
int main(int argc,char**argv){int n=argc>1?atoi(argv[1]):8;
  switch(n){case 6:Engine<6>().run();break;case 7:Engine<7>().run();break;
  case 8:Engine<8>().run();break;case 9:Engine<9>().run();break;case 10:Engine<10>().run();break;
  default:std::cerr<<"Supported ranks 6..10\n";return 2;}
}
