# Exact KL computation of P_{x,w} for the star group (center c, m leaves), x = product of leaves, w = x c x.
# Uses the geometric representation over Fractions, enumerates the Bruhat lower interval of w via subwords.
from fractions import Fraction
from itertools import combinations
import sys
sys.setrecursionlimit(10000)

def run(m):
    n = m+1  # generators: 0 = c, 1..m = leaves
    # Coxeter matrix: m(c,a_i)=3, m(a_i,a_j)=2
    B = [[Fraction(0)]*n for _ in range(n)]
    for i in range(n):
        B[i][i] = Fraction(1)
    for i in range(1,n):
        B[0][i] = B[i][0] = Fraction(-1,2)   # -cos(pi/3)
    def apply(s, v):  # s_i acting on vector v (coords in simple-root basis): v - 2B(v,alpha_s) alpha_s
        val = sum(v[j]*B[s][j] for j in range(n))
        w = list(v)
        w[s] -= 2*val
        return tuple(w)
    # represent elements by image of simple roots (tuple of tuples) -- faithful
    def mat_apply(s, M):  # left-multiply by s: columns s(M alpha_j)
        return tuple(apply(s, col) for col in M)
    I = tuple(tuple(Fraction(1) if i==j else Fraction(0) for i in range(n)) for j in range(n))
    def right_mult(M, s):  # columns of M*s: (M s)(alpha_j) = M(s alpha_j)
        cols = []
        for j in range(n):
            e = [Fraction(0)]*n; e[j]=Fraction(1)
            v = apply(s, tuple(e))
            # M v = sum v_i M(alpha_i)
            out = [sum(v[i]*M[i][k] for i in range(n)) for k in range(n)]
            cols.append(tuple(out))
        return tuple(cols)
    def is_right_descent(M, s):  # M(alpha_s) negative
        col = M[s]
        return any(c < 0 for c in col)
    word = list(range(1,n)) + [0] + list(range(1,n))  # x c x
    # enumerate subword products with lengths
    elems = {}
    def length(M):
        # count via descents: reduce
        l=0; N=M
        while True:
            d=[s for s in range(n) if is_right_descent(N,s)]
            if not d: return l
            N = right_mult(N, d[0]); l+=1
    frontier = {I}
    for s in word:
        new=set(frontier)
        for M in frontier:
            new.add(right_mult(M,s))
        frontier=new
    L = {M: length(M) for M in frontier}
    W_el = sorted(frontier, key=lambda M: L[M])
    wM = I
    for s in word: wM = right_mult(wM, s)
    assert L[wM]==len(word), (L[wM], len(word))
    # Bruhat order via subword property on fixed reduced word of each element; simpler: use criterion
    # x<=w iff ... we compute via lifting property recursively using descents: x<=w, s in R(w): x<=w iff min(x,xs)<=ws
    from functools import lru_cache
    @lru_cache(None)
    def leq(x, w):
        if x==w: return True
        if L[x] >= L[w]: return False
        if L[w]==0: return False
        s = next(t for t in range(n) if is_right_descent(w,t))
        ws = right_mult(w,s)
        if is_right_descent(x,s):
            return leq(right_mult(x,s), ws)
        return leq(x, ws)
    # KL polynomials as dict degree->int
    def padd(p,q,sign=1):
        r=dict(p)
        for k,v in q.items(): r[k]=r.get(k,0)+sign*v
        return {k:v for k,v in r.items() if v}
    def pshift(p,k): return {d+k:v for d,v in p.items()}
    @lru_cache(None)
    def P(x,w):
        if not leq(x,w): return ()
        if x==w: return ((0,1),)
        s = next(t for t in range(n) if is_right_descent(w,t))
        v = right_mult(w,s)
        xs = right_mult(x,s)
        c = 1 if is_right_descent(x,s) else 0
        # P_{x,w} = q^{1-c} P_{xs,v} + q^c P_{x,v} - sum_{z: s in R(z), z<v} mu(z,v) q^{(l(w)-l(z))/2} P_{x,z}
        res = padd(pshift(dict(P(xs,v)),1-c), pshift(dict(P(x,v)),c))
        for z in W_el:
            if L[z] < L[v] and is_right_descent(z,s) and leq(x,z) and leq(z,v):
                d = L[v]-L[z]
                if d%2==1:
                    pz = dict(P(z,v)); mu = pz.get((d-1)//2,0)
                    if mu:
                        res = padd(res, pshift({k:mu*val for k,val in dict(P(x,z)).items()}, (L[w]-L[z])//2), -1)
        return tuple(sorted(res.items()))
    xM=I
    for s in range(1,n): xM=right_mult(xM,s)
    p = dict(P(xM,wM))
    d = L[wM]-L[xM]
    mu = p.get((d-1)//2,0) if d%2 else 0
    print(f"m={m}: l(x)={L[xM]} l(w)={L[wM]} P_x,w = {sorted(p.items())} mu={mu} (interval size {sum(1 for z in W_el if leq(xM,z) and leq(z,wM))})")
for m in (2,3,4,5,6):
    run(m)
