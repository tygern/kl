// Independent exhaustive E8 terminal enumeration through D7 (not E7) cosets.
// E8: chain0--1--2--3--4--5--6, branch7 at2. J={1,...,7} is D7.
// Build: c++ -O3 -std=c++17 e8_d7_cosets.cpp -o e8-d7
#include <array>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
using Root=std::array<int,8>;
std::array<std::vector<int>,8> adj;
std::vector<Root> roots;
std::map<Root,int> rid;
int negs[240],positive[240],adds[240][240];
int get(uint64_t w,int s){return (w>>(8*s))&255;}
uint64_t set(uint64_t w,int s,int a){return (w&~(uint64_t(255)<<(8*s)))|(uint64_t(a)<<(8*s));}
uint64_t right(uint64_t w,int s){
  uint64_t v=set(w,s,negs[get(w,s)]);
  for(int t:adj[s]){int a=adds[get(w,s)][get(w,t)];assert(a>=0);v=set(v,t,a);}return v;
}
unsigned desc(uint64_t w){unsigned d=0;for(int s=0;s<8;s++)if(!positive[get(w,s)])d|=1u<<s;return d;}
bool weak(uint64_t w,unsigned allowed=255){
  unsigned d=desc(w)&allowed;
  for(int s=0;s<8;s++)if(d&(1u<<s))for(int t:adj[s])if(allowed&(1u<<t)){
    int sum=adds[get(w,s)][get(w,t)];assert(sum>=0);if(!positive[sum])return false;
  }return true;
}
bool independent(unsigned d){for(int s=0;s<8;s++)if(d&(1u<<s))for(int t:adj[s])if(d&(1u<<t))return false;return true;}
Root image(uint64_t w,const Root&a){Root b{};
  for(int s=0;s<8;s++)for(int k=0;k<8;k++)b[k]+=a[s]*roots[get(w,s)][k];return b;
}
std::array<uint8_t,240> permutation(uint64_t w){std::array<uint8_t,240> p{};
  for(int r=0;r<240;r++)p[r]=rid.at(image(w,roots[r]));return p;
}
uint64_t compose(const std::array<uint8_t,240>&p,uint64_t v){uint64_t w=0;for(int s=0;s<8;s++)w=set(w,s,p[get(v,s)]);return w;}
uint64_t identity(){uint64_t e=0;for(int s=0;s<8;s++)e=set(e,s,s);return e;}
uint64_t element(const std::vector<int>&word){uint64_t w=identity();for(int s:word)w=right(w,s);return w;}
uint64_t inverse(const std::vector<int>&word){uint64_t w=identity();for(auto s=word.rbegin();s!=word.rend();++s)w=right(w,*s);return w;}
std::vector<int> word(uint32_t w,const std::vector<uint32_t>&parent,const std::vector<uint8_t>&last){
  std::vector<int>a;while(w){a.push_back(last[w]);w=parent[w];}std::reverse(a.begin(),a.end());return a;
}
void printword(const std::vector<int>&a){std::cout<<'[';for(size_t i=0;i<a.size();i++)std::cout<<(i?",":"")<<a[i];std::cout<<']';}
bool leq(uint64_t x,unsigned lx,uint64_t w,unsigned lw){while(x!=w){if(lx>=lw)return false;
  unsigned d=desc(w);assert(d);int s=__builtin_ctz(d);
  if(desc(x)&(1u<<s)){x=right(x,s);lx--;}w=right(w,s);lw--;}return true;
}
int main(){
  for(int s=0;s<6;s++){adj[s].push_back(s+1);adj[s+1].push_back(s);}adj[2].push_back(7);adj[7].push_back(2);
  for(int s=0;s<8;s++){Root a{};a[s]=1;rid[a]=roots.size();roots.push_back(a);}
  for(size_t i=0;i<roots.size();i++)for(int s=0;s<8;s++){
    Root a=roots[i];a[s]=-a[s];for(int t:adj[s])a[s]+=roots[i][t];
    if(!rid.count(a)){rid[a]=roots.size();roots.push_back(a);}
  }
  assert(roots.size()==240);
  for(int i=0;i<240;i++){
    Root a=roots[i];positive[i]=true;for(int k=0;k<8;k++){positive[i]&=a[k]>=0;a[k]=-a[k];}negs[i]=rid.at(a);
    for(int j=0;j<240;j++){Root b{};for(int k=0;k<8;k++)b[k]=roots[i][k]+roots[j][k];auto it=rid.find(b);adds[i][j]=it==rid.end()?-1:it->second;}
  }
  // Enumerate only D7, retaining full E8 root actions.
  std::vector<uint64_t> els{identity()};std::unordered_map<uint64_t,uint32_t> ids;ids.reserve(400000);ids[identity()]=0;
  std::vector<uint32_t> parent{0};std::vector<uint8_t> last{0},len{0};
  for(uint32_t i=0;i<els.size();i++)for(int s=1;s<8;s++){
    uint64_t v=right(els[i],s);if(!ids.count(v)){uint32_t j=els.size();ids[v]=j;els.push_back(v);parent.push_back(i);last.push_back(s);len.push_back(len[i]+1);}
  }
  assert(els.size()==322560 && len.back()==42);
  std::vector<uint32_t> rightweak;std::array<uint64_t,43> d7hist{};
  for(uint32_t i=0;i<els.size();i++){d7hist[len[i]]++;if(weak(els[i],254))rightweak.push_back(i);}
  // The fundamental weight at deleted node0 is a norm-four integral vector.
  const Root beta={4,7,10,8,6,4,2,5};
  for(int s=0;s<8;s++){int pairing=2*beta[s];for(int t:adj[s])pairing-=beta[t];assert(pairing==(s==0));}
  std::vector<Root> orbit{beta};std::map<Root,uint32_t> orbit_ids;orbit_ids[beta]=0;
  std::vector<uint32_t> cp{0};std::vector<uint8_t> cg{0},cl{0};
  for(uint32_t i=0;i<orbit.size();i++)for(int s=0;s<8;s++){
    Root b=orbit[i];b[s]=-b[s];for(int t:adj[s])b[s]+=orbit[i][t];
    if(!orbit_ids.count(b)){orbit_ids[b]=orbit.size();orbit.push_back(b);cp.push_back(i);cg.push_back(s);cl.push_back(cl[i]+1);}
  }
  assert(orbit.size()==2160 && cl.back()==78);
  std::vector<uint64_t> cosets,cosetinv;std::vector<std::vector<int>> cosetword;std::array<uint64_t,79> cosethist{};
  std::vector<std::array<uint8_t,240>> cperm;
  for(uint32_t i=0;i<orbit.size();i++){
    // Left BFS: generators are recovered latest-first, with no reversal.
    std::vector<int>a;for(uint32_t j=i;j;j=cp[j])a.push_back(cg[j]);
    uint64_t w=element(a);assert(image(w,beta)==orbit[i]);assert(!(desc(w)&254));
    cosetword.push_back(a);cosets.push_back(w);cosetinv.push_back(inverse(a));cperm.push_back(permutation(w));cosethist[cl[i]]++;
  }
  // Independent root-ID FC catalogue for checking all eligible bottoms.
  std::vector<uint64_t> fcs{identity()};std::unordered_map<uint64_t,uint32_t> fid;fid[identity()]=0;
  std::vector<uint32_t> fp{0};std::vector<uint8_t> fg{0},fl{0};
  for(uint32_t i=0;i<fcs.size();i++){
    const uint64_t a=fcs[i];unsigned ad=desc(a);
    for(int s=0;s<8;s++)if(!(ad&(1u<<s))){uint64_t b=right(a,s);if(fid.count(b))continue;
      unsigned d=desc(b);if(!independent(d))continue;bool good=true;
      for(int t=0;t<8;t++)if(d&(1u<<t)){
        auto it=fid.find(right(b,t));if(it==fid.end()||fl[it->second]!=fl[i]){good=false;break;}}
      if(good){fid[b]=fcs.size();fcs.push_back(b);fp.push_back(i);fg.push_back(s);fl.push_back(fl[i]+1);}
    }
  }
  assert(fcs.size()==10846);
  std::vector<unsigned> fleft;for(uint32_t i=0;i<fcs.size();i++)fleft.push_back(desc(inverse(word(i,fp,fg))));
  struct Terminal{uint64_t w,inv;unsigned length;std::vector<int> wd;};std::vector<Terminal> terminals;
  uint64_t tested=0,rightsurvive=0;
  for(auto vi:rightweak){
    auto vw=word(vi,parent,last);uint64_t vinv=inverse(vw);auto viperm=permutation(vinv);
    for(uint32_t ai=0;ai<cosets.size();ai++){
      tested++;uint64_t w=compose(cperm[ai],els[vi]);if(!weak(w))continue;rightsurvive++;
      uint64_t inv=compose(viperm,cosetinv[ai]);if(!weak(inv))continue;
      auto wd=cosetword[ai];wd.insert(wd.end(),vw.begin(),vw.end());
      assert(element(wd)==w);assert(wd.size()==unsigned(cl[ai])+len[vi]);
      terminals.push_back({w,inv,unsigned(wd.size()),wd});
    }
  }
  std::sort(terminals.begin(),terminals.end(),[](const Terminal&a,const Terminal&b){return a.length<b.length || (a.length==b.length&&a.wd<b.wd);});
  unsigned commuting=0;for(auto&t:terminals){unsigned support=0;for(int s:t.wd)support|=1u<<s;
    if(unsigned(__builtin_popcount(support))==t.length&&independent(support))commuting++;}
  std::cout<<"{\"rank\":8,\"parabolic\":\"D7\",\"parabolic_order\":"<<els.size()
    <<",\"cosets\":"<<cosets.size()<<",\"parabolic_right_terminals\":"<<rightweak.size()
    <<",\"candidates_tested\":"<<tested<<",\"full_right_terminals\":"<<rightsurvive
    <<",\"terminal_count\":"<<terminals.size()<<",\"commuting_terminals\":"<<commuting
    <<",\"fc_count\":"<<fcs.size()<<",\"complete\":true,\"parabolic_histogram\":[";
  for(size_t k=0;k<d7hist.size();k++)std::cout<<(k?",":"")<<d7hist[k];
  std::cout<<"],\"coset_histogram\":[";for(size_t k=0;k<cosethist.size();k++)std::cout<<(k?",":"")<<cosethist[k];
  std::cout<<"],\"bad\":[";bool first=true;
  for(auto&t:terminals){unsigned support=0;for(int s:t.wd)support|=1u<<s;
    if(unsigned(__builtin_popcount(support))==t.length&&independent(support))continue;
    std::cout<<(first?"":",")<<"\n{\"word\":";first=false;printword(t.wd);
    unsigned rd=desc(t.w),ld=desc(t.inv);std::cout<<",\"length\":"<<t.length<<",\"Rmask\":"<<rd<<",\"Lmask\":"<<ld<<",\"bottoms\":[";bool ff=true;
    for(uint32_t x=0;x<fcs.size();x++)if((desc(fcs[x])&rd)==rd && (fleft[x]&ld)==ld && leq(fcs[x],fl[x],t.w,t.length)){
      std::cout<<(ff?"":",")<<"{\"word\":";ff=false;printword(word(x,fp,fg));std::cout<<",\"length\":"<<unsigned(fl[x])<<",\"rank\":"<<t.length-fl[x]<<'}';
    }std::cout<<"]}";
  }std::cout<<"\n]}\n";
}
