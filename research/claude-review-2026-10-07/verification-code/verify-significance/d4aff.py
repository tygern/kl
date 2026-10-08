# Exact KL computation in the star Coxeter group W_4 = affine D4: c=0, leaves 1..4.
from functools import lru_cache
from itertools import combinations
n=5
B=[[0]*n for _ in range(n)]
for i in range(n): B[i][i]=2
for i in range(1,5): B[0][i]=B[i][0]=-1
def refl(i):
    M=[[1 if r==c else 0 for c in range(n)] for r in range(n)]
    for c in range(n): M[i][c]-=B[c][i]  # s_i(e_c)=e_c - B(e_c,a_i) a_i
    return tuple(tuple(r) for r in M)
S=[refl(i) for i in range(n)]
def mul(A,Bm): return tuple(tuple(sum(A[r][k]*Bm[k][c] for k in range(n)) for c in range(n)) for r in range(n))
I=tuple(tuple(1 if r==c else 0 for c in range(n)) for r in range(n))
def col(M,j): return tuple(M[r][j] for r in range(n))
def neg(v): return all(x<=0 for x in v) and any(x<0 for x in v)
def rdesc(M): return [j for j in range(n) if neg(col(M,j))]
def ldesc(M): return rdesc(inv(M))
def length(M):
    L=0; M2=M
    while True:
        d=rdesc(M2)
        if not d: return L
        M2=mul(M2,S[d[0]]); L+=1
def inv(M):
    # elements are products of involutive integer matrices; compute by reducing
    word=[]; M2=M
    while True:
        d=rdesc(M2)
        if not d: break
        word.append(d[0]); M2=mul(M2,S[d[0]])
    R=I
    for i in word: R=mul(R,S[i])   # M = S[word reversed...]; careful
    # M * S[w0] * S[w1] ... = I  => M = S[w_last]...S[w0]; M^{-1} = S[w0]...S[w_last]
    return R
def word(M):
    w=[]; M2=M
    while True:
        d=rdesc(M2)
        if not d: break
        w.append(d[0]); M2=mul(M2,S[d[0]])
    return w[::-1]  # M = S[w[0]]...S[w[-1]]
def prod(w):
    M=I
    for i in w: M=mul(M,S[i])
    return M
x=prod([1,2,3,4]); w=prod([1,2,3,4,0,1,2,3,4])
assert length(x)==4 and length(w)==9
# lower ideal via subwords
def lower(M):
    wd=word(M); res=set()
    for mask in range(1<<len(wd)):
        res.add(prod([wd[i] for i in range(len(wd)) if mask>>i&1]))
    return res
low={}
for z in lower(w): low[z]=lower(z)
print("lower ideal size",len(low))
def leq(u,z): return u in low[z]
L={z:length(z) for z in low}
from functools import lru_cache
@lru_cache(None)
def P(u,z):
    if not leq(u,z): return ()
    if u==z: return (1,)
    s=ldesc(z)[0]; v=mul(S[s],z)
    su=mul(S[s],u); c=1 if L[su]<L[u] else 0
    def add(a,b):
        m=max(len(a),len(b)); return tuple((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(m))
    def shift(a,k): return (0,)*k+tuple(a)
    def scal(a,k): return tuple(k*t for t in a)
    res=add(shift(P(su,v),1-c), shift(P(u,v),c))
    for zz in low[v]:
        if zz==v: continue
        d=L[v]-L[zz]

    # mu(zz,v) with s zz < zz
    for zz in low[v]:
        if zz==v or not leq(u,zz): continue
        d=L[v]-L[zz]
        if d%2==0: continue
        if L[mul(S[s],zz)]>L[zz]: continue
        pz=P(zz,v); mu=pz[(d-1)//2] if len(pz)>(d-1)//2 else 0
        if mu: res=add(res, scal(shift(P(u,zz),(L[z]-L[zz])//2),-mu))
    while res and res[-1]==0: res=res[:-1]
    return res

p=P(x,w); print("P_{x,w} =",p, " ranks:", sorted(__import__('collections').Counter(L[z]-4 for z in low[w] if leq(x,z)).items()))
print("mu(x,w) =", p[(9-4-1)//2])
