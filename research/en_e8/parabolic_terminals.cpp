// Exact terminal enumeration by one-sided parabolic pruning.
// Build: c++ -std=c++17 -O3 research/en_e8/parabolic_terminals.cpp -o /private/tmp/e8-terminals
// Run: /private/tmp/e8-terminals 8 > research/en_e8/e8-terminals.json
// Optional third argument limits the number of cosets, making the result bounded.
// No enumeration or storage of the full E8 group is performed.
#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>
using State = uint64_t;
using Root = std::array<int,8>;
int n, omitted;
unsigned J;
std::vector<int> adj[8];
std::vector<std::pair<int,int>> edges;
std::vector<Root> roots;
std::map<Root,int> rid;
int adds[256][256];
uint8_t neg[256], reflect_root[8][256];
bool positive[256];
State identity;
uint64_t fc_star_checks_total=0;

inline int image(State w,int s){return (w>>(8*s))&255;}
inline State replace(State w,int s,int value){
    assert(value>=0&&value<256);
    return (w&~(State(255)<<(8*s)))|(State(value)<<(8*s));
}
inline State right(State w,int s){
    int a=image(w,s);
    State v=replace(w,s,neg[a]);
    for(int t:adj[s])v=replace(v,t,adds[a][image(w,t)]);
    return v;
}
State left(State w,int s){
    State v=0;
    for(int t=0;t<n;t++)v=replace(v,t,reflect_root[s][image(w,t)]);
    return v;
}
unsigned desc(State w){
    unsigned result=0;
    for(int s=0;s<n;s++)if(!positive[image(w,s)])result|=1u<<s;
    return result;
}
bool commuting(unsigned mask){
    for(auto [s,t]:edges)if((mask&(1u<<s))&&(mask&(1u<<t)))return false;
    return true;
}
bool weakright(State w,unsigned mask){
    unsigned d=desc(w)&mask;
    for(int s=0;s<n;s++)if(d&(1u<<s))
        for(int t:adj[s])if(mask&(1u<<t)){
            int z=adds[image(w,s)][image(w,t)];
            assert(z>=0);
            if(!positive[z])return false;
        }
    return true;
}
std::vector<int> reduced_word(State w){
    std::vector<int> reverse;
    while(w!=identity){
        auto d=desc(w);assert(d);
        int s=__builtin_ctz(d);reverse.push_back(s);w=right(w,s);
        assert(reverse.size()<=120);
    }
    return {reverse.rbegin(),reverse.rend()};
}
State inverse(const std::vector<int>& word){
    State v=identity;
    for(auto it=word.rbegin();it!=word.rend();++it)v=right(v,*it);
    return v;
}
bool leq(State x,int lx,State w,int lw){
    while(x!=w){
        if(lx>=lw)return false;
        auto d=desc(w);assert(d);int s=__builtin_ctz(d);
        if(!positive[image(x,s)]){x=right(x,s);--lx;}
        w=right(w,s);--lw;
    }
    return true;
}
void print_word(const std::vector<int>& word){
    std::cout<<'[';
    for(size_t i=0;i<word.size();i++)std::cout<<(i?",":"")<<word[i];
    std::cout<<']';
}
void initialize(){
    for(int s=0;s<n-2;s++)edges.emplace_back(s,s+1);
    edges.emplace_back(2,n-1);
    for(auto [s,t]:edges){adj[s].push_back(t);adj[t].push_back(s);}
    omitted=n-2;J=((1u<<n)-1)^(1u<<omitted);
    for(int s=0;s<n;s++){Root a{};a[s]=1;rid[a]=roots.size();roots.push_back(a);}
    for(size_t i=0;i<roots.size();i++)for(int s=0;s<n;s++){
        auto b=roots[i];b[s]=-b[s];for(int t:adj[s])b[s]+=roots[i][t];
        if(!rid.count(b)){rid[b]=roots.size();roots.push_back(b);}
    }
    assert(roots.size()==size_t(n==6?72:n==7?126:240));
    for(int i=0;i<int(roots.size());i++){
        auto b=roots[i];positive[i]=true;
        for(int s=0;s<n;s++){positive[i]&=b[s]>=0;b[s]=-b[s];}
        neg[i]=rid.at(b);
        for(int s=0;s<n;s++){
            b=roots[i];b[s]=-b[s];for(int t:adj[s])b[s]+=roots[i][t];
            reflect_root[s][i]=rid.at(b);
        }
        for(int j=0;j<int(roots.size());j++){
            Root c{};for(int s=0;s<n;s++)c[s]=roots[i][s]+roots[j][s];
            auto found=rid.find(c);adds[i][j]=found==rid.end()?-1:found->second;
        }
    }
    identity=0;for(int s=0;s<n;s++)identity=replace(identity,s,s);
}
struct Element {State w,inv;std::vector<int> word;};
std::vector<Element> parabolic_one_sided(uint64_t& parabolic_size){
    const size_t expected=n==6?1920:n==7?51840:2903040;
    std::vector<State> group{identity},inverses{identity};
    group.reserve(expected);inverses.reserve(expected);
    std::vector<uint32_t> parents{0};parents.reserve(expected);
    std::vector<uint8_t> last{0},lengths{0};last.reserve(expected);lengths.reserve(expected);
    std::unordered_map<State,uint32_t> ids;ids.reserve(expected);ids.emplace(identity,0);
    std::vector<Element> accepted;
    for(uint32_t w=0;w<group.size();w++){
        for(int s=0;s<n;s++)if(J&(1u<<s)){
            State v=right(group[w],s);
            if(ids.find(v)==ids.end()){
                auto id=group.size();ids.emplace(v,id);group.push_back(v);
                inverses.push_back(left(inverses[w],s));
                parents.push_back(w);last.push_back(s);lengths.push_back(lengths[w]+1);
            }
        }
        if(weakright(group[w],J)){
            std::vector<int> word;
            for(uint32_t p=w;p;p=parents[p])word.push_back(last[p]);
            std::reverse(word.begin(),word.end());
            assert(word.size()==lengths[w]);
            accepted.push_back({group[w],inverses[w],std::move(word)});
        }
    }
    assert(group.size()==expected);parabolic_size=group.size();
    std::cerr<<"Parabolic order "<<group.size()<<", one-sided terminals "<<accepted.size()<<'\n';
    return accepted;
}
std::vector<Element> cosets(){
    const size_t expected=n==6?27:n==7?56:240;
    std::vector<State> reps{identity};
    std::unordered_map<State,int> ids;ids.emplace(identity,0);
    for(size_t i=0;i<reps.size();i++)for(int s=0;s<n;s++){
        State a=left(reps[i],s);
        while(desc(a)&J)a=right(a,__builtin_ctz(desc(a)&J));
        if(ids.emplace(a,reps.size()).second)reps.push_back(a);
    }
    assert(reps.size()==expected);
    std::vector<Element> result;
    for(State a:reps){
        assert(!(desc(a)&J));auto word=reduced_word(a);
        result.push_back({a,inverse(word),std::move(word)});
    }
    std::cerr<<"Minimal right cosets "<<result.size()<<'\n';
    return result;
}
std::vector<Element> fully_commutative(){
    std::vector<Element> fc{{identity,identity,{}}};
    std::unordered_map<State,int> ids;ids.emplace(identity,0);
    for(size_t i=0;i<fc.size();i++)for(int s=0;s<n;s++){
        // Appending a descent cannot discover a longer FC element.
        if(!positive[image(fc[i].w,s)])continue;
        State w=right(fc[i].w,s);if(ids.count(w))continue;
        auto d=desc(w);if(!commuting(d))continue;
        bool ok=true;
        for(int t=0;t<n;t++)if(d&(1u<<t)){
            auto found=ids.find(right(w,t));
            if(found==ids.end())ok=false;
            else assert(fc[found->second].word.size()==fc[i].word.size());
        }
        if(!ok)continue;
        auto word=fc[i].word;word.push_back(s);
        auto inv=left(fc[i].inv,s);
        ids.emplace(w,fc.size());fc.push_back({w,inv,std::move(word)});
    }
    assert(fc.size()==size_t(n==6?662:n==7?2670:10846));
    uint64_t star_checks=0;
    for(const auto& x:fc)for(int side=0;side<2;side++)for(auto [s,t]:edges){
        State w=side?x.inv:x.w;
        unsigned mask=(1u<<s)|(1u<<t);
        if(__builtin_popcount(desc(w)&mask)!=1)continue;
        int valid=0;
        for(int u:{s,t}){
            State v=right(w,u);
            if(__builtin_popcount(desc(v)&mask)==1){++valid;assert(ids.count(v));}
        }
        assert(valid==1);++star_checks;
    }
    std::cerr<<"Fully commutative elements "<<fc.size()<<", star checks "<<star_checks<<'\n';
    fc_star_checks_total=star_checks;
    return fc;
}

#ifndef E8_TERMINAL_LIBRARY
int main(int argc,char**argv){
    n=argc>1?std::atoi(argv[1]):8;assert(n>=6&&n<=8);
    auto start=std::chrono::steady_clock::now();initialize();
    uint64_t parabolic_size=0;
    auto parabolic=parabolic_one_sided(parabolic_size);
    auto representatives=cosets();auto fc=fully_commutative();
    size_t limit=argc>2?std::min(representatives.size(),size_t(std::atoi(argv[2]))):representatives.size();
    const unsigned all=(1u<<n)-1;
    uint64_t candidates=0,right_terminals=0;
    std::vector<Element> terminals;
    for(size_t ci=0;ci<limit;ci++){
        const auto& a=representatives[ci];
        std::array<uint8_t,256> permutation{};
        for(int r=0;r<int(roots.size());r++)permutation[r]=r;
        for(auto it=a.word.rbegin();it!=a.word.rend();++it)
            for(int r=0;r<int(roots.size());r++)permutation[r]=reflect_root[*it][permutation[r]];
        for(int s=0;s<n;s++)assert(permutation[s]==image(a.w,s));
        for(const auto& v:parabolic){
            ++candidates;State w=0;
            for(int s=0;s<n;s++)w=replace(w,s,permutation[image(v.w,s)]);
            if(!weakright(w,all))continue;
            ++right_terminals;
            State wi=v.inv;
            for(auto it=a.word.rbegin();it!=a.word.rend();++it)wi=right(wi,*it);
            if(!weakright(wi,all))continue;
            auto word=a.word;word.insert(word.end(),v.word.begin(),v.word.end());
            assert(reduced_word(w).size()==word.size());
            assert(inverse(word)==wi);
            terminals.push_back({w,wi,std::move(word)});
        }
        if(ci%20==19||ci+1==limit)
            std::cerr<<"Cosets "<<ci+1<<'/'<<limit<<", candidates "<<candidates
                     <<", right terminals "<<right_terminals<<", terminals "<<terminals.size()<<'\n';
    }
    std::sort(terminals.begin(),terminals.end(),[](const auto&a,const auto&b){
        return a.word.size()!=b.word.size()?a.word.size()<b.word.size():a.word<b.word;
    });
    std::unordered_map<State,int> fc_ids;
    for(size_t i=0;i<fc.size();i++)fc_ids.emplace(fc[i].w,i);
    unsigned commuting_count=0;std::vector<Element> bads;
    for(auto& b:terminals){
        if(fc_ids.count(b.w)){
            unsigned support=0;for(int s:b.word)support|=1u<<s;
            assert(__builtin_popcount(support)==int(b.word.size())&&commuting(support));
            ++commuting_count;
        }else bads.push_back(std::move(b));
    }
    auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"{\"type\":\"E"<<n<<"\",\"complete\":"<<(limit==representatives.size()?"true":"false")
             <<",\"method\":\"minimal right parabolic cosets times one-sided terminal parabolic elements\""
             <<",\"parabolic_omitted_generator\":"<<omitted
             <<",\"parabolic_order\":"<<parabolic_size
             <<",\"parabolic_right_terminal_count\":"<<parabolic.size()
             <<",\"coset_count\":"<<representatives.size()<<",\"cosets_processed\":"<<limit
             <<",\"candidates_tested\":"<<candidates<<",\"right_terminal_count\":"<<right_terminals
             <<",\"fc_count\":"<<fc.size()<<",\"fc_star_checks\":"<<fc_star_checks_total
             <<",\"commuting_terminals\":"<<commuting_count
             <<",\"noncommuting_terminal_count\":"<<bads.size()<<",\"seconds\":"<<seconds<<",\"bad\":[";
    for(size_t i=0;i<bads.size();i++){
        const auto& b=bads[i];unsigned rd=desc(b.w),ld=desc(b.inv);
        std::cout<<(i?",":"")<<"{\"word\":";print_word(b.word);
        std::cout<<",\"length\":"<<b.word.size()<<",\"Rmask\":"<<rd<<",\"Lmask\":"<<ld<<",\"bottoms\":[";
        bool first=true;
        for(const auto& x:fc)if((desc(x.w)&rd)==rd&&(desc(x.inv)&ld)==ld&&leq(x.w,x.word.size(),b.w,b.word.size())){
            std::cout<<(first?"":",")<<"{\"word\":";first=false;print_word(x.word);
            std::cout<<",\"length\":"<<x.word.size()<<",\"rank\":"<<b.word.size()-x.word.size()<<'}';
        }
        std::cout<<"]}";
    }
    std::cout<<"]}\n";
}
#endif
