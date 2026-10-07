from targeted import star
from collections import deque
from itertools import combinations
import json
from pathlib import Path

def factor(g,w,u):
    todo=deque([(w,())]);seen={w};l=g.length(u);ui=g.inv(u)
    while todo:
        p,suffix=todo.popleft()
        a=p
        for s in g.word(ui):a=g.right(a,s)
        if g.length(a)+l==g.length(p):
            return dict(prefix=g.word(a),factor=g.word(u),suffix=suffix)
        for s in g.desc(p):
            q=g.right(p,s)
            if q not in seen:seen.add(q);todo.append((q,(s,)+suffix))
    return None

if __name__=='__main__':
    g=star(4);A=(1,2,3,4);out=[]
    u=g.elt((1,0,2,1,0,1))
    assert g.length(u)==6 and set(g.desc(u))=={0,1,2}
    for k in range(1,6):
        word=(A+(0,))*k+A;w=g.elt(word)
        out.append(dict(k=k,length=g.length(w),factor=factor(g,w,u)))
    Path(__file__).with_name('parabolic_factor_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
