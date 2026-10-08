# Independent enumeration of terminal elements in E6 (paper numbering: chain 0-1-2-3-4, node 5 attached to 2)
# Representation: permutation of the 72 roots (root-index storage). Terminal: for every right descent s and t~s, ws has no descent at t (i.e. w(alpha_s+alpha_t)>0); same for inverse.
import sys
from collections import deque
n=6
edges={(0,1),(1,2),(2,3),(3,4),(2,5)}
adj={i:set() for i in range(n)}
for a,b in edges: adj[a].add(b); adj[b].add(a)
C=[[2 if i==j else (-1 if j in adj[i] else 0) for j in range(n)] for i in range(n)]
def refl(i,v):
    val=sum(C[i][j]*v[j] for j in range(n))
    w=list(v); w[i]-=val; return tuple(w)
# generate roots
simple=[tuple(1 if i==j else 0 for i in range(n)) for j in range(n)]
roots=set(simple); fr=list(simple)
while fr:
    nf=[]
    for v in fr:
        for i in range(n):
            u=refl(i,v)
            if u not in roots: roots.add(u); nf.append(u)
    fr=nf
roots=sorted(roots); idx={r:k for k,r in enumerate(roots)}
N=len(roots); assert N==72, N
pos=[all(c>=0 for c in r) for r in roots]
gen=[tuple(idx[refl(i,r)] for r in roots) for i in range(n)]   # permutation: s_i(root k)
simp_idx=[idx[s] for s in simple]
sumidx={}
for a,b in edges:
    v=tuple(simple[a][j]+simple[b][j] for j in range(n)); sumidx[(a,b)]=sumidx[(b,a)]=idx[v]
ident=tuple(range(N))
def rmul(w,i):  # w*s_i : root k -> w(s_i(k))
    g=gen[i]; return tuple(w[g[k]] for k in range(N))
def lmul(i,w):
    g=gen[i]; return tuple(g[w[k]] for k in range(N))
# BFS
dist={ident:0}; q=deque([ident]); words={ident:()}
while q:
    w=q.popleft()
    for i in range(n):
        u=rmul(w,i)
        if u not in dist:
            dist[u]=dist[w]+1; words[u]=words[w]+(i,); q.append(u)
print("group order",len(dist), "max length",max(dist.values()))
def rdesc(w): return [i for i in range(n) if not pos[w[simp_idx[i]]]]
def inv(w):
    r=[0]*N
    for k in range(N): r[w[k]]=k
    return tuple(r)
def right_terminal(w):
    for s in rdesc(w):
        for t in adj[s]:
            if not pos[w[sumidx[(s,t)]]]: return False
    return True
term=[w for w in dist if right_terminal(w) and right_terminal(inv(w))]
def commuting_product(w):
    d=rdesc(w)
    return dist[w]==len(d) and all((a,b) not in edges for a in d for b in d)
nonc=[w for w in term if not commuting_product(w)]
print("terminals",len(term),"commuting",len(term)-len(nonc),"noncommuting",len(nonc))
for w in nonc: print("  word",''.join(map(str,words[w])),"len",dist[w],"R",rdesc(w),"L",rdesc(inv(w)))
