"""Independent support-sensitive check of the maximum-I lemma on supplied FC words.
This is a finite countercheck, not a proof or a certification of the input catalogues.
"""
import json
from pathlib import Path
from functools import lru_cache
ROOT=Path(__file__).resolve().parents[2]
for n in range(6,10):
    edges=[(i,i+1) for i in range(n-2)]+[(2,n-1)]
    adj=[[] for _ in range(n)]
    for s,t in edges: adj[s].append(t);adj[t].append(s)
    E=tuple(tuple(int(i==j) for i in range(n)) for j in range(n))
    def right(a,s):
        b=list(a);b[s]=tuple(-v for v in a[s])
        for t in adj[s]:b[t]=tuple(u+v for u,v in zip(a[t],a[s]))
        return tuple(b)
    def word(w):
        a=E
        for s in w:a=right(a,s)
        return a
    def desc(a):
        d=0
        for i,c in enumerate(a):
            assert any(c) and (min(c)>=0 or max(c)<=0)
            if max(c)<=0:d|=1<<i
        return d
    independents=[m for m in range(1<<n) if not any(m>>s&1 and m>>t&1 for s,t in edges)]
    @lru_cache(None)
    def alpha(support):return max(m.bit_count() for m in independents if m&support==m)
    rows=json.loads((ROOT/f'research/en_independent/e{n}-fc.json').read_text())['elements']
    checked=0
    for row in rows:
        w=row['word'];a=word(w);ai=word(w[::-1]);d=desc(a)&desc(ai)
        support=sum(1<<s for s in set(w))
        if d.bit_count()==alpha(support):
            assert len(w)==d.bit_count() and set(w)=={s for s in range(n) if d>>s&1}
            checked+=1
    print(f'E{n}: {len(rows)} words checked against support-sensitive maximum-I conclusion; {checked} equality cases, no counterexample')
