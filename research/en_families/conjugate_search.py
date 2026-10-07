"""Bounded terminal search by involution conjugation and descent conjugation.
A heuristic candidate search, never a terminal classification.
"""
from targeted import en
from itertools import combinations
import random,json,time,argparse
from pathlib import Path


def left(g,w,s):
    out=[]
    for col in w:
        a=list(col);a[s]=-col[s]+sum(col[t] for t in g.adj[s]);out.append(tuple(a))
    return tuple(out)

def conj(g,w,s):return left(g,g.right(w,s),s)

def indep(g,I):return not any(t in I for s in I for t in g.adj[s])

def reduce_involution(g,w,rng):
    # Prefer decreasing conjugations from noncommuting suffixes. A terminal
    # with commuting descents passes immediately. Others can get stuck.
    while True:
        ds=g.desc(w)
        if not indep(g,ds):return None
        moves=[s for s in ds if any(t in g.desc(g.right(w,s)) for t in g.adj[s])]
        if not moves:return w
        rng.shuffle(moves)
        good=False
        for s in moves:
            z=conj(g,w,s)
            if z!=w:
                # Both s descents imply drop2 unless the first removal loses
                # the other descent. Root-sign test is enough.
                if s in g.desc(left(g,w,s)):
                    w=z;good=True;break
        if not good:return None

def search(n,trials,seed,maxsteps):
    g=en(n);rng=random.Random(seed);start=time.monotonic();seen={};attempted=0
    indeps=[I for k in range(2,n) for I in combinations(range(n),k) if indep(g,I)]
    # All finite-E7 seed tops embedded in E_n, then conjugated.
    base=json.loads(Path('research/broad_exceptional/e7_bad.json').read_text())['bad']
    seeds=[g.elt([n-1 if s==6 else s for s in row['word']]) for row in base]
    for trial in range(trials):
        w=rng.choice(seeds) if trial%2 else g.elt(rng.choice(indeps))
        for step in range(rng.randint(0,maxsteps)):
            w=conj(g,w,rng.randrange(n))
        v=reduce_involution(g,w,rng);attempted+=1
        if v is not None and v not in seen:
            word=g.word(v)
            if len(word)>len(set(word)):
                row=dict(rank=n,word=word,length=len(word),support=sorted(set(word)),R=g.desc(v),full_support=len(set(word))==n,trial=trial)
                seen[v]=row
                print(json.dumps(row),flush=True)
        if trial%1000==999:
            # Caches outside stored exact candidate matrices are disposable.
            g.right.cache_clear();g.desc.cache_clear();g.word.cache_clear()
    return dict(rank=n,trials=trials,seed=seed,max_conjugation_steps=maxsteps,seconds=time.monotonic()-start,terminals=list(seen.values()))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('rank',type=int);p.add_argument('--trials',type=int,default=10000);p.add_argument('--steps',type=int,default=50);p.add_argument('--seed',type=int,default=934);a=p.parse_args()
    out=search(a.rank,a.trials,a.seed,a.steps)
    Path(__file__).with_name(f'conjugates_E{a.rank}_{a.trials}_{a.seed}.json').write_text(json.dumps(out,indent=2)+'\n')
    print('DONE',out['seconds'],len(out['terminals']),flush=True)
