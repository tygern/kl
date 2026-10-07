from targeted import en,affine_d,bipartition
import json
from pathlib import Path
out=[]
for typ,make,ranks in [('E',en,range(6,17)),('affineD',affine_d,range(4,13))]:
    for n in ranks:
        g=make(n);colors=bipartition(g)
        for initial in (0,1):
            word=()
            for layers in range(1,22):
                word+=colors[(initial+layers-1)%2]
                w=g.elt(word)
                if g.terminal(w) and g.length(w)>len(set(word)):
                    out.append(dict(type=typ,rank=n,initial=initial,layers=layers,word=word,length=g.length(w),R=g.desc(w),reduced=g.length(w)==len(word)))
        g.right.cache_clear();g.desc.cache_clear();g.word.cache_clear();g.inv.cache_clear()
Path(__file__).with_name('bipartite_terminal_scan.json').write_text(json.dumps(out,indent=2)+'\n')
for row in out:print(row['type'],row['rank'],row['initial'],row['layers'],row['length'],row['reduced'])
