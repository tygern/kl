"""Star Coxeter systems W_m (center c, m leaves), via exact integer root-lattice representation.
Generators act on Z^{m+1} (basis alpha_c, alpha_1..alpha_m) by s_i(v) = v - <v, alpha_i^vee> alpha_i,
with Cartan matrix a_ij = <alpha_j, alpha_i^vee> (symmetric here)."""
import sys
from math import comb
from functools import lru_cache
def run(m):
    N = m+1
    C = [[0]*N for _ in range(N)]
    for i in range(N): C[i][i] = 2
    for i in range(1, N): C[0][i] = C[i][0] = -1
    def refl(s, v):  # s acts on vector v (coords in simple roots)
        coef = sum(C[s][j]*v[j] for j in range(N))
        return tuple(v[j] - (coef if j == s else 0) for j in range(N))
    # element = tuple of images of simple roots (matrix columns)
    E = tuple(tuple(1 if i == j else 0 for i in range(N)) for j in range(N))
    def mul_s(w, s):  # w * s : (w s)(alpha_j) = w(s(alpha_j))
        def apply(v):
            out = [0]*N
            for j in range(N):
                if v[j]:
                    for i in range(N): out[i] += v[j]*w[j][i]
            return tuple(out)
        return tuple(apply(refl(s, E[j])) for j in range(N))
    def neg(v): return all(c <= 0 for c in v) and any(c < 0 for c in v)
    def rdesc(w): return [s for s in range(N) if neg(w[s])]
    @lru_cache(None)
    def length(w):
        d = rdesc(w)
        return 0 if not d else 1 + length(mul_s(w, d[0]))
    def ev(word):
        w = E
        for s in word: w = mul_s(w, s)
        return w
    leaves = list(range(1, N))
    xw = leaves; ww = leaves + [0] + leaves
    x, w = ev(xw), ev(ww)
    assert length(x) == m and length(w) == 2*m+1
    ideal = {}
    def Ideal(v):
        if v in ideal: return ideal[v]
        d = rdesc(v)
        if not d: r = frozenset([v])
        else:
            I = Ideal(mul_s(v, d[0])); r = frozenset(I | {mul_s(y, d[0]) for y in I})
        ideal[v] = r; return r
    def leq(a, b): return a in Ideal(b)
    def padd(a, b, shift=0, scale=1):
        r = list(a) + [0]*max(0, len(b)+shift-len(a))
        for k, c in enumerate(b): r[k+shift] += scale*c
        while r and r[-1] == 0: r.pop()
        return tuple(r)
    KL = {}
    def P(a, b):
        if (a, b) in KL: return KL[(a, b)]
        if not leq(a, b): r = ()
        elif a == b: r = (1,)
        else:
            s = rdesc(b)[0]; bs = mul_s(b, s); as_ = mul_s(a, s)
            c = 1 if length(as_) < length(a) else 0
            r = padd(padd((), P(as_, bs), 1-c), P(a, bs), c)
            for z in Ideal(bs):
                if leq(a, z) and length(mul_s(z, s)) < length(z):
                    dz = length(bs) - length(z)
                    if dz % 2 == 1:
                        pz = P(z, bs); k = (dz-1)//2
                        mu = pz[k] if k < len(pz) else 0
                        if mu: r = padd(r, P(a, z), (length(b)-length(z))//2, -mu)
        KL[(a, b)] = r; return r
    Rc = {}
    def R(a, b):
        if (a, b) in Rc: return Rc[(a, b)]
        if not leq(a, b): r = ()
        elif a == b: r = (1,)
        else:
            s = rdesc(b)[0]; bs = mul_s(b, s); as_ = mul_s(a, s)
            if length(as_) < length(a): r = R(as_, bs)
            else: r = padd(padd((), R(a, bs), 0, -1), padd(R(a, bs), R(as_, bs)), 1)
        Rc[(a, b)] = r; return r
    I = [z for z in Ideal(w) if leq(x, z)]
    rv = [0]*(m+2)
    for z in I: rv[length(z)-m] += 1
    pol = P(x, w)
    pred = tuple(comb(m, j) - (comb(m, j-1) if j >= 1 else 0) for j in range(m//2+1))
    # reciprocity
    d = m+1; lhs = ()
    for z in I:
        r = R(x, z); p = P(z, w); prod = ()
        for k, c in enumerate(p): prod = padd(prod, r, k, c)
        lhs = padd(lhs, prod)
    rhs = [0]*(d+1)
    for k, c in enumerate(pol): rhs[d-k] += c
    while rhs and rhs[-1] == 0: rhs.pop()
    mu = pol[(d-1)//2] if d % 2 == 1 and (d-1)//2 < len(pol) else 0
    allP1 = all(P(z, w) == (1,) for z in I if z != x)
    print(f"m={m}: |[e,w]|={len(Ideal(w))} |[x,w]|={len(I)} (3^m+1={3**m+1}) rank vector {rv}; P={pol} predicted {pred} match={pol==pred}; mu={mu}; reciprocity={lhs==tuple(rhs)}; P(z,w)=1 for all faces z>x: {allP1}")
    # Boolean-interval check: all intervals [y,z] with x<y<=z<=w have R=(q-1)^rank and P=1
    ok = True
    for y in I:
        if y == x: continue
        for z in I:
            if leq(y, z):
                dd = length(z)-length(y)
                exp = ()
                e1 = (1,)
                for _ in range(dd): e1 = padd(padd((), e1, 0, -1), e1, 1)
                if R(y, z) != e1 or P(y, z) != (1,): ok = False
    print(f"   all face-to-face subintervals have R=(q-1)^d and P=1: {ok}")
for m in range(1, 7): run(m)
