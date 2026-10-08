# Independent FC catalogue of generalized E_n (chain 0..n-2, node n-1 attached to 2)
# and exhaustive test of Lemma maxI: for FC x not a commuting product, |L(x) cap R(x)| < alpha(supp x).
import sys
from itertools import combinations
n = int(sys.argv[1])
adj = [set() for _ in range(n)]
def edge(a,b): adj[a].add(b); adj[b].add(a)
for i in range(n-2): edge(i,i+1)
edge(n-1,2)
nb = [tuple(sorted(adj[i])) for i in range(n)]
Id = tuple(tuple(1 if a==b else 0 for a in range(n)) for b in range(n))
def rmul(M,i):
    cols=list(M); ai=cols[i]
    cols[i]=tuple(-a for a in ai)
    for j in nb[i]:
        cols[j]=tuple(a+b for a,b in zip(cols[j],ai))
    return tuple(cols)
def lmul(M,i):
    out=[]
    for v in M:
        vi=-v[i]+sum(v[j] for j in nb[i])
        out.append(v[:i]+(vi,)+v[i+1:])
    return tuple(out)
def neg(v):
    for a in v:
        if a<0: return True
        if a>0: return False
    return False
def rdesc(M): return frozenset(j for j in range(n) if neg(M[j]))
def independent(D): return all(b not in adj[a] for a,b in combinations(sorted(D),2))
acache={}
def alpha(mask):
    if mask in acache: return acache[mask]
    K=[j for j in range(n) if mask>>j&1]
    best=0
    for r in range(len(K),0,-1):
        if r<=best: break
        for C in combinations(K,r):
            if independent(C): best=r; break
    acache[mask]=best; return best
level={Id:(Id,0)}
total=1; l=0; viol=[]; comm=1; maxl=0
while level:
    nxt={}
    for M,(Minv,supp) in level.items():
        R=rdesc(M)
        for s in range(n):
            if s in R: continue
            M2=rmul(M,s)
            if M2 in nxt: continue
            R2=rdesc(M2)
            if not independent(R2): continue
            ok=True
            for t in R2:
                if t!=s and rmul(M2,t) not in level: ok=False; break
            if not ok: continue
            Minv2=lmul(Minv,s); supp2=supp|(1<<s)
            nxt[M2]=(Minv2,supp2)
            L2=rdesc(Minv2); D=L2&R2
            sz=bin(supp2).count('1'); suppset=frozenset(j for j in range(n) if supp2>>j&1)
            iscomm = (l+1==sz and independent(suppset))
            a=alpha(supp2)
            if iscomm:
                comm+=1
                if D!=suppset: viol.append(('comm',l+1,M2))
            else:
                if len(D)>=a: viol.append(('maxI',l+1,sorted(suppset),sorted(D),a))
    l+=1
    total+=len(nxt); 
    if nxt: maxl=l
    level=nxt
print(f"E{n}: FC count={total}, max FC length={maxl}, commuting products={comm}, violations={len(viol)}")
for v in viol[:10]: print(v)
