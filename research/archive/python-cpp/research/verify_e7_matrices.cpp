// Independent E7 terminal verification using integral matrices, not root IDs.
// Build: c++ -std=c++17 -O3 research/verify_e7_matrices.cpp -o verify-e7
// Run: ./verify-e7 > results/e7-independent-certificate.json
#include <array>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <cstdint>
using Matrix=std::array<int8_t,49>;
struct Hash {size_t operator()(const Matrix&a)const {
    size_t h=1469598103934665603ULL;
    for(auto c:a)h=(h^uint8_t(c))*1099511628211ULL;
    return h;
}};
const std::array<std::pair<int,int>,6> edges={{{0,1},{1,2},{2,3},{3,4},{4,5},{2,6}}};
std::array<std::vector<int>,7> adj;
Matrix multiply(const Matrix&a,int s){
    Matrix b=a;
    for(int k=0;k<7;k++)b[s*7+k]=-a[s*7+k];
    for(int t:adj[s])for(int k=0;k<7;k++)b[t*7+k]=a[t*7+k]+a[s*7+k];
    return b;
}
unsigned descents(const Matrix&a){
    unsigned mask=0;
    for(int s=0;s<7;s++){
        bool pos=false,neg=false;
        for(int k=0;k<7;k++){pos|=a[s*7+k]>0;neg|=a[s*7+k]<0;}
        assert(pos!=neg);
        if(neg)mask|=1<<s;
    }
    return mask;
}
int main(){
    for(auto [s,t]:edges){adj[s].push_back(t);adj[t].push_back(s);}
    constexpr uint32_t expected=2903040;
    Matrix identity{};for(int s=0;s<7;s++)identity[s*7+s]=1;
    std::vector<Matrix> el;el.reserve(expected);el.push_back(identity);
    std::unordered_map<Matrix,uint32_t,Hash> ids;ids.reserve(expected);ids[identity]=0;
    std::vector<std::array<uint32_t,7>> action;action.reserve(expected);
    std::vector<uint32_t> parent{0};parent.reserve(expected);
    std::vector<uint8_t> last{0},len{0},rd,fc;
    last.reserve(expected);len.reserve(expected);rd.reserve(expected);fc.reserve(expected);
    for(uint32_t w=0;w<el.size();w++){
        std::array<uint32_t,7> row{};
        for(int s=0;s<7;s++){
            Matrix b=multiply(el[w],s);
            auto it=ids.find(b);
            if(it==ids.end()){
                uint32_t v=el.size();ids.emplace(b,v);el.push_back(b);
                parent.push_back(w);last.push_back(s);len.push_back(len[w]+1);row[s]=v;
            }else row[s]=it->second;
        }
        action.push_back(row);unsigned d=descents(el[w]);rd.push_back(d);
        bool good=true;
        for(int s=0;s<7;s++)if(d&(1<<s)){assert(len[row[s]]+1==len[w]);good&=fc[row[s]];}
        for(auto [s,t]:edges)good&=!((d&(1<<s))&&(d&(1<<t)));
        fc.push_back(good);
    }
    assert(el.size()==expected);
    ids.clear();ids.rehash(0);el.clear();el.shrink_to_fit();
    std::vector<uint32_t> inverse(expected);
    for(uint32_t w=0;w<expected;w++){
        uint32_t z=0;
        for(uint32_t v=w;v;v=parent[v])z=action[z][last[v]];
        inverse[w]=z;
    }
    auto left=[&](uint32_t w,int s){return inverse[action[inverse[w]][s]];};
    auto weak=[&](uint32_t w){for(int s=0;s<7;s++)if(rd[w]&(1<<s))
        for(int t:adj[s])if(rd[action[w][s]]&(1<<t))return false;return true;};
    auto leq=[&](uint32_t x,uint32_t w){while(x!=w){
        if(len[x]>=len[w])return false;
        int s=__builtin_ctz(unsigned(rd[w]));
        if(rd[x]&(1<<s))x=action[x][s];w=action[w][s];}return true;};
    auto printword=[&](uint32_t w){std::vector<int> word;for(;w;w=parent[w])word.push_back(last[w]);
        std::reverse(word.begin(),word.end());std::cout<<'[';
        for(size_t k=0;k<word.size();k++)std::cout<<(k?",":"")<<word[k];std::cout<<']';};
    std::vector<uint32_t> fcs,bads;uint64_t star_checks=0;std::array<unsigned,64> levels{};
    unsigned commuting_terminals=0;
    for(uint32_t w=0;w<expected;w++){
        levels[len[w]]++;
        if(fc[w]){
            fcs.push_back(w);
            for(int side=0;side<2;side++)for(auto [s,t]:edges){
                auto ds=[&](uint32_t v){return unsigned(rd[side?inverse[v]:v]);};
                if(__builtin_popcount(ds(w)&((1<<s)|(1<<t)))!=1)continue;
                int valid=0;
                for(int u:{s,t}){auto v=side?left(w,u):action[w][u];
                    if(__builtin_popcount(ds(v)&((1<<s)|(1<<t)))==1){valid++;assert(fc[v]);}}
                assert(valid==1);star_checks++;
            }
        }
        if(weak(w)&&weak(inverse[w])){
            if(!fc[w])bads.push_back(w);
            else {unsigned support=0;for(uint32_t v=w;v;v=parent[v])support|=1<<last[v];
                assert(__builtin_popcount(support)==len[w]);
                for(auto [s,t]:edges)assert(!((support&(1<<s))&&(support&(1<<t))));
                commuting_terminals++;}
        }
    }
    assert(fcs.size()==2670);assert(bads.size()==4);
    std::cout<<"{\"representation\":\"independent integer matrices\",\"order\":"<<expected
        <<",\"fc_count\":"<<fcs.size()<<",\"fc_star_checks\":"<<star_checks
        <<",\"commuting_terminals\":"<<commuting_terminals<<",\"length_distribution\":[";
    for(int k=0;k<64;k++)std::cout<<(k?",":"")<<levels[k];std::cout<<"],\"bad\":[";
    for(size_t k=0;k<bads.size();k++){
        uint32_t w=bads[k];std::cout<<(k?",":"")<<"{\"word\":";printword(w);
        std::cout<<",\"length\":"<<int(len[w])<<",\"Rmask\":"<<int(rd[w])
                 <<",\"Lmask\":"<<int(rd[inverse[w]])<<",\"bottoms\":[";
        bool first=true;
        for(auto x:fcs)if((rd[x]&rd[w])==rd[w]&&(rd[inverse[x]]&rd[inverse[w]])==rd[inverse[w]]&&leq(x,w)){
            std::cout<<(first?"":",")<<"{\"word\":";first=false;printword(x);
            std::cout<<",\"length\":"<<int(len[x])<<",\"rank\":"<<int(len[w]-len[x])<<'}';
        }
        std::cout<<"]}";
    }
    std::cout<<"]}\n";
}
