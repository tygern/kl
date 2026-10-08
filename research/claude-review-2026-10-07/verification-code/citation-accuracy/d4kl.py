# Exact KL polynomials in D4 via standard recursion. Generators as signed perms on {1..4}.
# Gern labeling: s1=(1,-2)(-1,2), s_i swaps i-1,i. Leaves: s1,s2,s4; center s3.
import itertools
n=4
def gen(i):
    w=list(range(1,n+1))
    if i==1: w[0]=-2; w[1]=-1
    else: w[i-2],w[i-1]=w[i-1],w[i-2]
    return tuple(w)
def ap(w,i): return w[i-1] if i>0 else -w[-i-1]
def mul(u,v): return tuple(ap(u,v[i]) for i in range(n))
def length(w):
    L=0
    for i in range(n):
        for j in range(i+1,n):
            if w[i]>w[j]: L+=1
            if -w[i]>w[j]: L+=1
    return L
e=tuple(range(1,n+1)); S=[gen(i) for i in range(1,n+1)]
# BFS the group
G={e}; frontier=[e]
while frontier:
    nf=[]
    for w in frontier:
        for s in S:
            u=mul(w,s)
            if u not in G: G.add(u); nf.append(u)
    frontier=nf
G=sorted(G,key=length); assert len(G)==192
idx={w:i for i,w in enumerate(G)}
ell={w:length(w) for w in G}
def inv(w):
    r=[0]*n
    for i in range(n):
        v=w[i]
        if v>0: r[v-1]=i+1
        else: r[-v-1]=-(i+1)
    return tuple(r)
# Bruhat order via reflections: x<w iff exists chain; use standard: compute via reduced words subword? simpler: x<=w iff using descent recursion:
# x<=w: if w=e: x==e. pick s with ws<w: x<=w iff (xs<x ? xs<=ws : x<=ws) -- standard lifting property version: if s in R(w): x<=w iff min(x,xs)<=ws
from functools import lru_cache
def rdesc(w): return [s for s in S if ell[mul(w,s)]<ell[w]]
@lru_cache(None)
def leq(x,w):
    if ell[x]>ell[w]: return False
    if w==e: return x==e
    s=rdesc(w)[0]
    ws=mul(w,s); xs=mul(x,s)
    xm = xs if ell[xs]<ell[x] else x
    return leq(xm,ws)
# polynomials as tuple of coefficients
def padd(a,b):
    m=max(len(a),len(b)); return tuple((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(m))
def pshift(a,k): return (0,)*k+tuple(a)
def pscale(a,c): return tuple(c*t for t in a)
def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0: a.pop()
    return tuple(a)
P={}
def mu(z,w):
    p=KL(z,w); d=ell[w]-ell[z]
    if d%2==0 or d<=0: return 0
    k=(d-1)//2
    return p[k] if k<len(p) else 0
def KL(x,w):
    if (x,w) in P: return P[(x,w)]
    if not leq(x,w): r=(0,)
    elif x==w: r=(1,)
    else:
        s=rdesc(w)[0]; ws=mul(w,s); xs=mul(x,s)
        c=1 if ell[xs]<ell[x] else 0
        r=padd(pshift(KL(xs,ws),1-c), pshift(KL(x,ws),c))
        for z in G:
            if ell[z]<ell[ws] and ell[mul(z,s)]<ell[z] and leq(x,z) and leq(z,ws):
                m=mu(z,ws)
                if m: r=padd(r, pscale(pshift(KL(x,z),(ell[w]-ell[z])//2), -m))
        r=trim(r)
    P[(x,w)]=r; return r
s1,s2,s3,s4=S
x=mul(mul(s1,s2),s4)
b=mul(mul(x,s3),x)
print("l(x)=",ell[x],"l(b)=",ell[b],"x<=b",leq(x,b))
print("P_{x,b} =",KL(x,b))
print("P_{e,b} =",KL(e,b))
# check Gern Example 2.1.2 element w = s1 s2 s4 s3 s1 s2 s4 equals b
