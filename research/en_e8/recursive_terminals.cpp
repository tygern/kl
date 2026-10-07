// Recursively prune right-terminal parabolic factors along a subgroup chain.
// This reuses exact root arithmetic, not the full-parabolic enumeration.
// Build: c++ -std=c++17 -O3 research/en_e8/recursive_terminals.cpp -o /private/tmp/e8-recursive
// Run: /private/tmp/e8-recursive 8 > research/en_e8/e8-recursive.json
#define E8_TERMINAL_LIBRARY
#include "parabolic_terminals.cpp"
#include <sys/resource.h>

struct Stage {unsigned K,J;int added;uint64_t order,cosets,previous,candidates,right_count;};

std::vector<Element> minimal_cosets(unsigned K,unsigned H,size_t expected){
    std::vector<State> reps{identity};
    std::unordered_map<State,int> ids;ids.emplace(identity,0);
    for(size_t i=0;i<reps.size();i++)for(int s=0;s<n;s++)if(K&(1u<<s)){
        State a=left(reps[i],s);
        while(desc(a)&H)a=right(a,__builtin_ctz(desc(a)&H));
        if(ids.emplace(a,reps.size()).second)reps.push_back(a);
    }
    assert(reps.size()==expected);
    std::vector<Element> result;
    for(State a:reps){
        assert(!(desc(a)&H));auto word=reduced_word(a);
        for(int s:word)assert(K&(1u<<s));
        result.push_back({a,inverse(word),std::move(word)});
    }
    return result;
}

std::vector<Element> extend_right_terminals(const std::vector<Element>& prior,
                                           unsigned K,unsigned H,size_t index){
    auto reps=minimal_cosets(K,H,index);
    std::vector<Element> result;
    std::unordered_map<State,int> seen;
    for(const auto& a:reps){
        std::array<uint8_t,256> permutation{};
        for(int r=0;r<int(roots.size());r++)permutation[r]=r;
        for(auto it=a.word.rbegin();it!=a.word.rend();++it)
            for(int r=0;r<int(roots.size());r++)permutation[r]=reflect_root[*it][permutation[r]];
        for(int s=0;s<n;s++)assert(permutation[s]==image(a.w,s));
        // Computational check of the sign-preservation used in the proof.
        for(int r=0;r<int(roots.size());r++)if(positive[r]){
            bool in_subsystem=true;
            for(int s=0;s<n;s++)if(!(H&(1u<<s))&&roots[r][s])in_subsystem=false;
            if(in_subsystem)assert(positive[permutation[r]]);
        }
        for(const auto& v:prior){
            State w=0;
            for(int s=0;s<n;s++)w=replace(w,s,permutation[image(v.w,s)]);
            // A minimal representative preserves all signs in the H subsystem.
            assert(weakright(w,H));
            if(!weakright(w,K))continue;
            assert(seen.emplace(w,result.size()).second);
            State wi=v.inv;
            for(auto it=a.word.rbegin();it!=a.word.rend();++it)wi=right(wi,*it);
            auto word=a.word;word.insert(word.end(),v.word.begin(),v.word.end());
            assert(reduced_word(w).size()==word.size());
            assert(inverse(word)==wi);
            result.push_back({w,wi,std::move(word)});
        }
    }
    return result;
}

int main(int argc,char**argv){
    n=argc>1?std::atoi(argv[1]):8;assert(n>=6&&n<=8);
    auto start=std::chrono::steady_clock::now();initialize();
    std::vector<int> additions{n-1,2,3,1,0,4};
    for(int s=5;s<=n-2;s++)additions.push_back(s);
    const std::vector<uint64_t> indices{2,3,4,8,10,27,56,240};
    assert(additions.size()==size_t(n));
    std::vector<Element> right_terminals{{identity,identity,{}}};
    std::vector<Stage> stages;
    unsigned K=0;uint64_t order=1,total_candidates=0;
    for(size_t level=0;level<additions.size();level++){
        unsigned H=K;K|=1u<<additions[level];
        auto previous=right_terminals.size();
        auto next=extend_right_terminals(right_terminals,K,H,indices[level]);
        order*=indices[level];total_candidates+=indices[level]*previous;
        stages.push_back({K,H,additions[level],order,indices[level],previous,
                          indices[level]*previous,next.size()});
        std::cerr<<"Rank "<<level+1<<", order "<<order<<", cosets "<<indices[level]
                 <<", candidates "<<indices[level]*previous<<", right terminals "<<next.size()<<'\n';
        right_terminals=std::move(next);
    }
    assert(order==uint64_t(n==6?51840:n==7?2903040:696729600));
    auto fc=fully_commutative();std::unordered_map<State,int> fc_ids;
    for(size_t i=0;i<fc.size();i++)fc_ids.emplace(fc[i].w,i);
    unsigned commuting_count=0;std::vector<Element> bads;
    for(auto& b:right_terminals)if(weakright(b.inv,K)){
        if(fc_ids.count(b.w)){
            unsigned support=0;for(int s:b.word)support|=1u<<s;
            assert(__builtin_popcount(support)==int(b.word.size())&&commuting(support));
            ++commuting_count;
        }else bads.push_back(std::move(b));
    }
    std::sort(bads.begin(),bads.end(),[](const auto&a,const auto&b){
        return a.word.size()!=b.word.size()?a.word.size()<b.word.size():a.word<b.word;
    });
    auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    rusage usage{};int status=getrusage(RUSAGE_SELF,&usage);assert(status==0);
    uint64_t peak_bytes=usage.ru_maxrss;
#ifndef __APPLE__
    peak_bytes*=1024;
#endif
    std::cout<<"{\"type\":\"E"<<n<<"\",\"complete\":true,\"method\":\"recursive one-sided parabolic pruning\""
             <<",\"group_order\":"<<order<<",\"total_candidates_tested\":"<<total_candidates
             <<",\"right_terminal_count\":"<<right_terminals.size()<<",\"fc_count\":"<<fc.size()
             <<",\"fc_star_checks\":"<<fc_star_checks_total
             <<",\"commuting_terminals\":"<<commuting_count<<",\"noncommuting_terminal_count\":"<<bads.size()
             <<",\"seconds\":"<<seconds<<",\"peak_resident_bytes\":"<<peak_bytes<<",\"stages\":[";
    for(size_t i=0;i<stages.size();i++){
        const auto& s=stages[i];
        std::cout<<(i?",":"")<<"{\"subgroup_mask\":"<<s.K<<",\"smaller_subgroup_mask\":"<<s.J
                 <<",\"added_generator\":"<<s.added<<",\"order\":"<<s.order<<",\"coset_count\":"<<s.cosets
                 <<",\"previous_right_terminal_count\":"<<s.previous<<",\"candidates_tested\":"<<s.candidates
                 <<",\"right_terminal_count\":"<<s.right_count<<'}';
    }
    std::cout<<"],\"bad\":[";
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
