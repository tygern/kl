#!/usr/bin/env python3
"""Independent single-pair check of a Kazhdan-Lusztig polynomial via R-polynomials.

Given reduced words x, w in a simply laced type (A/D/E, manuscript numbering for E),
compute P_{x,w} from the defining identity (Kazhdan-Lusztig 1979, (2.2.a)):

   q^{l(w)-l(x)} P_{x,w}(q^{-1}) - P_{x,w}(q) = sum_{x<y<=w} R_{x,y}(q) P_{y,w}(q),

solving top-down for P_{y,w}, y in the interval [x,w], with R-polynomials from
   R_{y,z} = R_{ys,zs}                         if s in R(z), ys < y
   R_{y,z} = (q-1) R_{y,zs} + q R_{ys,zs}       if s in R(z), ys > y.
Elements are 8x8 integer matrices (images of simple roots); no part of klmu.cpp is reused.
Exact integer arithmetic throughout.  Usage: rpoly_check.py <type> <rank> <xword> <wword>
"""
import sys
from functools import lru_cache
sys.setrecursionlimit(100000)

t, n = sys.argv[1], int(sys.argv[2])
xw, ww = [int(c) for c in sys.argv[3]], [int(c) for c in sys.argv[4]]
M = [[2]*n for _ in range(n)]
for i in range(n): M[i][i] = 1
def edge(i, j): M[i][j] = M[j][i] = 3
if t == 'A':
    for i in range(n-1): edge(i, i+1)
elif t == 'D':
    for i in range(n-2): edge(i, i+1)
    edge(n-3, n-1)
elif t == 'E':
    for i in range(n-2): edge(i, i+1)
    edge(2, n-1)
Acart = [[2 if i == j else (-1 if M[i][j] == 3 else 0) for j in range(n)] for i in range(n)]

# element = tuple of n column tuples (images of simple roots)
ident = tuple(tuple(1 if k == j else 0 for k in range(n)) for j in range(n))
def rmul(m, i):          # m * s_i
    cols = list(m)
    ci = m[i]
    for j in range(n):
        if j != i and Acart[i][j]:
            cols[j] = tuple(m[j][k] - Acart[i][j]*ci[k] for k in range(n))
    cols[i] = tuple(-c for c in ci)
    return tuple(cols)
def neg(col):
    for c in col:
        if c: return c < 0
    return False
def rdesc(m): return [i for i in range(n) if neg(m[i])]
# positive roots for length computation
simple = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
def refl(i, v):
    pair = sum(Acart[i][k]*v[k] for k in range(n)); v = list(v); v[i] -= pair; return tuple(v)
roots = set(simple); stack = list(simple)
while stack:
    r = stack.pop()
    for i in range(n):
        r2 = refl(i, r)
        if all(c >= 0 for c in r2) and r2 not in roots: roots.add(r2); stack.append(r2)
roots = sorted(roots)
def apply(m, v): return tuple(sum(v[j]*m[j][k] for j in range(n)) for k in range(n))
@lru_cache(maxsize=None)
def length(m): return sum(1 for r in roots if neg(apply(m, r)))
def inverse(m):
    # columns of m^{-1}: solve; since m is a signed permutation of roots up to basis, use linear algebra over Q via fractions
    from fractions import Fraction
    A_ = [[Fraction(m[j][k]) for j in range(n)] + [Fraction(1 if k == r else 0) for r in range(n)] for k in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A_[r][c] != 0); A_[c], A_[piv] = A_[piv], A_[c]
        pv = A_[c][c]; A_[c] = [a/pv for a in A_[c]]
        for r in range(n):
            if r != c and A_[r][c] != 0:
                f = A_[r][c]; A_[r] = [a - f*b for a, b in zip(A_[r], A_[c])]
    inv = [[A_[k][n+j] for k in range(n)] for j in range(n)]
    return tuple(tuple(int(inv[j][k]) for k in range(n)) for j in range(n))
@lru_cache(maxsize=None)
def ldesc(m): return tuple(rdesc(inverse(m)))
def lmul(m, i): return inverse(rmul(inverse(m), i))
def word_to_elt(wd):
    m = ident
    for s in wd: m = rmul(m, s)
    return m
@lru_cache(maxsize=None)
def leq(x, w):
    while True:
        if x == w: return True
        if length(x) >= length(w): return False
        s = ldesc(w)[0]
        w = lmul(w, s)
        if s in ldesc(x): x = lmul(x, s)

x, w = word_to_elt(xw), word_to_elt(ww)
assert length(x) == len(xw) and length(w) == len(ww), "words not reduced"
assert leq(x, w), "x not <= w"
# interval [x,w]: enumerate subwords of w's reduced word (elements <= w), keep those >= x
below = {ident}
for s in ww:
    below |= {rmul(m, s) for m in below}
interval = [y for y in below if leq(x, y)]
interval.sort(key=length)
print(f"l(x)={length(x)} l(w)={length(w)} |[e,w]|={len(below)} |[x,w]|={len(interval)}")

def padd(a, b, shift=0, coef=1):
    res = list(a) + [0]*max(0, len(b)+shift-len(a))
    for k, c in enumerate(b): res[k+shift] += coef*c
    while res and res[-1] == 0: res.pop()
    return res
@lru_cache(maxsize=None)
def R(y, z):
    if y == z: return (1,)
    if not leq(y, z): return ()
    s = rdesc(z)[-1]
    zs = rmul(z, s); ys = rmul(y, s)
    if length(ys) < length(y): return R(ys, zs)
    a = R(y, zs); b = R(ys, zs)
    res = padd(padd([], a, 1), a, 0, -1)       # (q-1) a
    res = padd(res, b, 1)                      # + q b
    return tuple(res)

P = {w: [1]}
for y in reversed(interval):        # decreasing length
    if y == w: continue
    d = length(w) - length(y)
    rhs = []
    for z in interval:
        if z != y and length(z) > length(y) and leq(y, z):
            r = R(y, z); pz = P[z]
            prod = [0]*(len(r)+len(pz)-1) if r and pz else []
            for i, a in enumerate(r):
                for j, b in enumerate(pz): prod[i+j] += a*b
            rhs = padd(rhs, prod)
    # q^d P(q^-1) - P(q) = rhs ; P has degree <= (d-1)/2, so coefficients of q^k for k < d/2 of -P(q) are determined by rhs:
    # [q^k] rhs = -P_k for k <= (d-1)/2  (the q^d P(1/q) part has degrees >= d - (d-1)/2 > (d-1)/2)
    Py = [-(rhs[k] if k < len(rhs) else 0) for k in range((d-1)//2 + 1)]
    while Py and Py[-1] == 0: Py.pop()
    # consistency: full identity must hold
    lhs = padd([0]*(d+1), [], 0)
    lhs = [0]*(d+1)
    for k, c in enumerate(Py): lhs[d-k] += c; lhs[k] -= c
    while lhs and lhs[-1] == 0: lhs.pop()
    assert lhs == rhs, (lhs, rhs, Py)
    P[y] = Py
d = length(w) - length(x)
print("P_{x,w} =", P[x], " mu =", P[x][(d-1)//2] if d % 2 == 1 and len(P[x]) > (d-1)//2 else 0)
