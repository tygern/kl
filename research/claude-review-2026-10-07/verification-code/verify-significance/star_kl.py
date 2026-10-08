# Exact KL polynomials P_{x,w} for star Coxeter group W_m (center c, m leaves), x = product of leaves, w = x c x.
# Enumerate the Bruhat interval [e, w] via the subword property on the fixed reduced word, represent elements
# via the geometric representation (exact integers), compute P via standard recursion restricted to the interval.
import sys
from fractions import Fraction
from itertools import combinations

def run(m):
    n = m+1  # generators: 0 = c, 1..m = leaves
    # Coxeter matrix: m(c,a_i)=3, m(a_i,a_j)=2. Bilinear form B(e_s,e_t) = -cos(pi/m_st): -1/2 for m=3, 0 for m=2, 1 for s=t.
    # Use 2*B to keep integers: 2B(s,s)=2, 2B(c,a_i)=-1, 2B(a_i,a_j)=0.
    def B2(s,t):
        if s==t: return 2
        if s==0 or t==0: return -1
        return 0
    # Represent element by its matrix acting on root coordinates (columns = images of simple roots), integers.
    # s(v) = v - 2B(v,e_s) e_s = v - (sum_t v_t * B2(t,s)) e_s
    def apply_s(s, v):
        coef = sum(v[t]*B2(t,s) for t in range(n))
        w = list(v); w[s] -= coef
        return tuple(w)
    def id_mat():
        return tuple(tuple(1 if i==j else 0 for i in range(n)) for j in range(n))  # columns: images of e_j
    def mult_right(M, s):
        # M*s: image of e_j under M s = M(s(e_j)); s(e_j) = e_j - B2(j,s) e_s
        cols = []
        for j in range(n):
            v = list(M[j])
            c = B2(j,s)
            if c != 0:
                v = [v[i] - c*M[s][i] for i in range(n)]
            cols.append(tuple(v))
        return tuple(cols)
    def is_pos(v):
        return all(a>=0 for a in v) and any(a>0 for a in v)
    def right_descent(M, s):
        # s is a right descent iff M(e_s) is negative
        return not is_pos(M[s])
    # reduced word for w = x c x: leaves 1..m, then 0, then leaves 1..m
    word = list(range(1,m+1)) + [0] + list(range(1,m+1))
    L = len(word)
    # Enumerate all subwords, compute element matrices, dedupe. Keep the shortest-length as length (via reduce by descents).
    elems = {}
    from itertools import product
    for mask in range(1<<L):
        M = id_mat()
        for i in range(L):
            if mask>>i & 1:
                M = mult_right(M, word[i])
        if M not in elems:
            # compute length: number of positive roots sent negative = count via descent-reduction
            N = M; ln = 0
            while True:
                found = False
                for s in range(n):
                    if right_descent(N, s):
                        N = mult_right(N, s); ln += 1; found = True; break
                if not found: break
            elems[M] = ln
    E = list(elems.keys())
    idx = {M:i for i,M in enumerate(E)}
    length = [elems[M] for M in E]
    # Bruhat order via subword: y <= z iff y is a subword of a reduced word of z. Compute via standard recursion:
    # Use the characterization: for s in D_R(z), y <= z iff min(y, ys) <= zs. Base: y<=e iff y=e.
    import functools
    sys.setrecursionlimit(10000)
    rmul = {}
    def rm(i, s):
        key=(i,s)
        if key not in rmul:
            M = mult_right(E[i], s)
            rmul[key] = idx[M] if M in idx else None
        return rmul[key]
    @functools.lru_cache(None)
    def leq(i, j):
        if length[i] > length[j]: return False
        if i == j: return True
        if length[j] == 0: return False
        z = E[j]
        for s in range(n):
            if right_descent(z, s):
                js = rm(j, s)
                if right_descent(E[i], s):
                    is_ = rm(i, s)
                    if is_ is None: return False  # ys outside interval -> cannot be <= zs (zs in interval; interval is lower set)
                    return leq(is_, js)
                else:
                    return leq(i, js)
        return False
    # KL polynomials via recursion: P_{y,z} for z with right descent s, zs = v (v < z):
    # P_{y,z} = q^{1-c} P_{ys,v} + q^c P_{y,v} - sum_{u: v>=u, us<u} mu(u,v) q^{(l(z)-l(u))/2} P_{y,u}, c = 1 if ys<y else 0
    # Polynomials as dict deg->int.
    def padd(a,b):
        r=dict(a)
        for k,v in b.items(): r[k]=r.get(k,0)+v
        return {k:v for k,v in r.items() if v}
    def pscale(a,c,shift):
        return {k+shift:v*c for k,v in a.items() if v*c}
    @functools.lru_cache(None)
    def P(i, j):
        if not leq(i,j): return ()
        if i == j: return ((0,1),)
        z = E[j]
        s = next(t for t in range(n) if right_descent(z,t))
        v = rm(j, s)
        ys = rm(i, s)
        y_desc = right_descent(E[i], s)
        res = {}
        if y_desc:
            # c = 1: q^0 P_{ys,v} + q P_{y,v}
            res = padd(res, dict(P(ys, v)))
            res = padd(res, pscale(dict(P(i, v)), 1, 1))
        else:
            # c = 0: q P_{ys,v} + P_{y,v}
            if ys is not None:
                res = padd(res, pscale(dict(P(ys, v)), 1, 1))
            res = padd(res, dict(P(i, v)))
        for u in range(len(E)):
            if leq(u, v) and right_descent(E[u], s) and leq(i, u):
                d = length[v]-length[u]
                if d % 2 == 1:
                    pu = dict(P(u, v))
                    mu = pu.get((d-1)//2, 0)
                    if mu:
                        res = padd(res, pscale(dict(P(i,u)), -mu, (length[j]-length[u])//2))
        return tuple(sorted(res.items()))
    # x and w
    Mx = id_mat()
    for s in range(1,m+1): Mx = mult_right(Mx, s)
    Mw = id_mat()
    for s in word: Mw = mult_right(Mw, s)
    ix, iw = idx[Mx], idx[Mw]
    poly = dict(P(ix, iw))
    d = length[iw]-length[ix]
    mu = poly.get((d-1)//2,0) if d%2==1 else 0
    interval = sum(1 for i in range(len(E)) if leq(ix,i) and leq(i,iw))
    print(f"m={m}: |[e,w]|={len(E)} |[x,w]|={interval} l(x)={length[ix]} l(w)={length[iw]} P_{{x,w}}={poly} mu={mu}")

for m in [2,3,4,5,6]:
    run(m)
