#!/usr/bin/env python3
"""Independent single-pair check of P_{x,w} via R-polynomials (right-sided only, no inverses).

   q^{l(w)-l(x)} P_{x,w}(q^{-1}) - P_{x,w}(q) = sum_{x<y<=w} R_{x,y}(q) P_{y,w}(q)      (KL79, 2.2.a)
   R_{y,z} = R_{ys,zs}  if s in R(z), ys<y ;  R_{y,z} = (q-1)R_{y,zs} + qR_{ys,zs} otherwise.
Bruhat test (right-sided lifting): for s in R(w): x<=w iff (xs<=ws if s in R(x) else x<=ws).
Usage: rpoly_check2.py <type A|D|E> <rank> <xword> <wword>   (E numbering: chain 0..n-2, node n-1 at node 2)
"""
import sys
sys.setrecursionlimit(1000000)
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
C = [[2 if i == j else (-1 if M[i][j] == 3 else 0) for j in range(n)] for i in range(n)]
ident = tuple(tuple(1 if k == j else 0 for k in range(n)) for j in range(n))
simple = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
def refl(i, v):
    pair = sum(C[i][k]*v[k] for k in range(n)); v = list(v); v[i] -= pair; return tuple(v)
roots = set(simple); st = list(simple)
while st:
    r = st.pop()
    for i in range(n):
        r2 = refl(i, r)
        if all(c >= 0 for c in r2) and r2 not in roots: roots.add(r2); st.append(r2)
roots = sorted(roots)
def neg(col):
    for c in col:
        if c: return c < 0
    return False
RM, LEN, RD = {}, {}, {}
def rmul(m, i):
    key = (m, i)
    r = RM.get(key)
    if r is None:
        cols = list(m); ci = m[i]
        for j in range(n):
            if j != i and C[i][j]: cols[j] = tuple(m[j][k] - C[i][j]*ci[k] for k in range(n))
        cols[i] = tuple(-c for c in ci); r = tuple(cols); RM[key] = r
    return r
def length(m):
    l = LEN.get(m)
    if l is None:
        l = sum(1 for r in roots if neg(tuple(sum(r[j]*m[j][k] for j in range(n)) for k in range(n)))); LEN[m] = l
    return l
def rdesc(m):
    d = RD.get(m)
    if d is None: d = tuple(i for i in range(n) if neg(m[i])); RD[m] = d
    return d
LEQ = {}
def leq(x, w):
    key = (x, w); r = LEQ.get(key)
    if r is not None: return r
    x0, w0 = x, w
    while True:
        if x == w: r = True; break
        if length(x) >= length(w): r = False; break
        s = rdesc(w)[0]
        w = rmul(w, s)
        if s in rdesc(x): x = rmul(x, s)
    LEQ[key] = r; return r
def padd(a, b, shift=0, coef=1):
    res = list(a) + [0]*max(0, len(b)+shift-len(a))
    for k, c in enumerate(b): res[k+shift] += coef*c
    while res and res[-1] == 0: res.pop()
    return tuple(res)
RP = {}
def R(y, z):
    key = (y, z); r = RP.get(key)
    if r is not None: return r
    if y == z: r = (1,)
    elif not leq(y, z): r = ()
    else:
        s = rdesc(z)[-1]; zs = rmul(z, s); ys = rmul(y, s)
        if length(ys) < length(y): r = R(ys, zs)
        else:
            a, b = R(y, zs), R(ys, zs)
            r = padd(padd(padd((), a, 1), a, 0, -1), b, 1)
    RP[key] = r; return r
x = ident
for s in xw: x = rmul(x, s)
w = ident
for s in ww: w = rmul(w, s)
assert length(x) == len(xw) and length(w) == len(ww), "words not reduced"
assert leq(x, w), "x not <= w"
below = {ident}
for s in ww: below |= {rmul(m, s) for m in below}
interval = sorted((y for y in below if leq(x, y)), key=length)
print(f"l(x)={length(x)} l(w)={length(w)} |[e,w]|={len(below)} |[x,w]|={len(interval)}", flush=True)
P = {w: (1,)}
for y in reversed(interval):
    if y == w: continue
    d = length(w) - length(y); rhs = ()
    for z in interval:
        if z != y and length(z) > length(y) and leq(y, z):
            r = R(y, z); pz = P[z]
            prod = [0]*(len(r)+len(pz)-1)
            for i, a in enumerate(r):
                for j, b in enumerate(pz): prod[i+j] += a*b
            rhs = padd(rhs, prod)
    Py = [-(rhs[k] if k < len(rhs) else 0) for k in range((d-1)//2 + 1)]
    while Py and Py[-1] == 0: Py.pop()
    lhs = [0]*(d+1)
    for k, c in enumerate(Py): lhs[d-k] += c; lhs[k] -= c
    while lhs and lhs[-1] == 0: lhs.pop()
    assert tuple(lhs) == rhs, ("identity violated", y)
    P[y] = tuple(Py)
d = length(w) - length(x)
print("P_{x,w} =", list(P[x]), " mu =", P[x][(d-1)//2] if d % 2 == 1 and len(P[x]) > (d-1)//2 else 0)
