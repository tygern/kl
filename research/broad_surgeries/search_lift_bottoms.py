"""Check every rank-four bottom below every lifted bad D6 top."""
from coxeter_surgery import Coxeter,d_edges,bad_word
from functools import lru_cache
import json


def run(n,edge):
    source=Coxeter(n,d_edges(n))
    edges=d_edges(n)
    i,j,_=edges[edge]
    edges[edge]=(i,j,2)
    target=Coxeter(n,edges)
    top=source.element(bad_word(n))
    botword=(0,1)+tuple(range(3,n,2))
    src_ranks=source.ranks(source.interval(source.element(botword),top))
    @lru_cache(None)
    def lifts(v):
        if v==source.e:
            return {target.e:()}
        result={}
        for s in range(n):
            if not source.positive(v,s):
                for z,word in lifts(source.right(v,s)).items():
                    result[target.right(z,s)]=word+(s,)
        return result
    for v,word in lifts(top).items():
        vertices=sorted(target.lower(v),key=target.length)
        bots=[u for u in vertices if target.length(u)==len(botword)]
        ids={u:i for i,u in enumerate(bots)}
        ranks=[[0]*len(src_ranks) for u in bots]
        descendants={}
        for u in vertices:
            length=target.length(u)
            if length<len(botword): continue
            if length==len(botword):
                down={ids[u]}
            else:
                expression=target.word(u)
                down=set()
                for k in range(length):
                    z=target.element(expression[:k]+expression[k+1:])
                    if target.length(z)==length-1:
                        down.update(descendants[z])
            descendants[u]=down
            for idx in down:
                ranks[idx][length-len(botword)]+=1
        for idx,count in enumerate(ranks):
            if count==src_ranks:
                print(json.dumps({'edge':edge,'top':word,'bottom':target.word(bots[idx]),'ranks':count}),flush=True)
        print(json.dumps({'edge':edge,'top':word,'bottoms_checked':len(bots)}),flush=True)


if __name__=='__main__':
    for edge in range(5):
        run(6,edge)
