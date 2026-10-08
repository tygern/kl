import sys, time
sys.path.insert(0,'/Users/ic/workspace/kl/computations')
sys.path.insert(0,'/Users/ic/workspace/kl/research/en_families')
from sparse_kl import SparseCoxeter
import targeted
g=SparseCoxeter('D',6)
x=(-1,-2,4,3,6,5); w=(-1,-6,3,-4,5,-2)
def word(z):
    out=[]
    while g.length(z):
        s=g.descents(z)[0]; out.append(s); z=g.right(z,s)
    return tuple(reversed(out))
wx,ww=word(x),word(w)
print('words',wx,ww,len(wx),len(ww))
# check commutation structure of signed-permutation generators to build graph
import itertools
edges=[]
for s,t in itertools.combinations(range(6),2):
    e=(1,2,3,4,5,6)
    st=g.right(g.right(e,s),t); ts=g.right(g.right(e,t),s)
    if st!=ts: edges.append((s,t))
print('edges',edges)
m=targeted.Coxeter(6,edges)
X=m.elt(wx); W=m.elt(ww)
print('matrix lengths',m.length(X),m.length(W),'leq',m.leq(X,W))
t=time.perf_counter(); p=m.kl(X,W); print('matrix-engine P =',p,'seconds',round(time.perf_counter()-t,2))
print('sparse P =',g.kl(x,w))
print('kl_values cached', len(m.kl_values))
L=m.lower(W); print('lower ideal', len(L), 'interval', sum(1 for z in L if m.leq(X,z)))
from collections import Counter
print('rank vector', [Counter(m.length(z)-m.length(X) for z in L if m.leq(X,z))[i] for i in range(12)])
# every right descent of W gives same value
for s in m.desc(W):
    v=m.right(W,s); xs=m.right(X,s); c=int(m.length(xs)<m.length(X))
    q=targeted.add((),m.kl(xs,v),1-c); q=targeted.add(q,m.kl(X,v),c)
    for z,d,mu in m.corrections(W,s):
        if m.leq(X,z): q=targeted.add(q,m.kl(X,z),d,-mu)
    print('descent',s,q)
