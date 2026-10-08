#!/usr/bin/env python3
"""Naive, independent Kazhdan-Lusztig polynomial computation for a small finite
Coxeter group, for cross-checking klmu.cpp.  Group elements are realised as
permutations of the root system (vectors in the simple-root basis), the Bruhat
order by the subword/lifting recursion, and KL polynomials by the RIGHT-sided
recursion with full-table storage (no extremal-pair reductions):

  w = v s, s in R(w):  P_{x,w} = q^{1-c} P_{xs,v} + q^c P_{x,v}
                        - sum_{z < v, zs < z} mu(z,v) q^{(l(w)-l(z))/2} P_{x,z},
  c = 1 if xs < x else 0.
Writes lines "<word x> <word w> c0 c1 ..." like klmu's dump (words are canonical
reduced words obtained by stripping the smallest right descent each time).
"""
import sys
from collections import deque

def coxeter_data(t, n):
    M = [[2]*n for _ in range(n)]
    for i in range(n): M[i][i] = 1
    def edge(i, j, m): M[i][j] = M[j][i] = m
    if t == 'A':
        for i in range(n-1): edge(i, i+1, 3)
    elif t == 'B':
        for i in range(n-1): edge(i, i+1, 3)
        edge(n-2, n-1, 4)
    elif t == 'D':
        for i in range(n-2): edge(i, i+1, 3)
        edge(n-3, n-1, 3)
    elif t == 'E':
        for i in range(n-2): edge(i, i+1, 3)
        edge(2, n-1, 3)
    else:
        raise ValueError
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            A[i][j] = 2 if i == j else (0 if M[i][j] == 2 else -1)
    if t == 'B':
        A[n-2][n-1] = -2
    return M, A

def main():
    t, n = sys.argv[1], int(sys.argv[2])
    out = sys.argv[3]
    M, A = coxeter_data(t, n)
    # simple reflection on a root vector
    def refl(i, v):
        pair = sum(A[i][k]*v[k] for k in range(n))
        v = list(v); v[i] -= pair; return tuple(v)
    # positive roots
    simple = [tuple(1 if k == i else 0 for k in range(n)) for i in range(n)]
    roots = set(simple); dq = deque(simple)
    while dq:
        r = dq.popleft()
        for i in range(n):
            r2 = refl(i, r)
            if all(c >= 0 for c in r2) and r2 not in roots:
                roots.add(r2); dq.append(r2)
    roots = sorted(roots)
    allroots = roots + [tuple(-c for c in r) for r in roots]
    ridx = {r: k for k, r in enumerate(allroots)}
    # group elements as permutations of allroots (tuples); generators
    gens = [tuple(ridx[refl(i, r)] for r in allroots) for i in range(n)]
    ident = tuple(range(len(allroots)))
    def compose(p, g):   # (p * g)(r) = p(g(r)) : right multiplication w*s_i acts as w(s_i(r))
        return tuple(p[g[k]] for k in range(len(p)))
    elems = {ident: 0}; lst = [ident]; length = [0]; rm = []
    q = deque([0])
    while q:
        w = q.popleft(); rm.append([None]*n)
        for i in range(n):
            ws = compose(lst[w], gens[i])
            if ws not in elems:
                elems[ws] = len(lst); lst.append(ws); length.append(length[w]+1); q.append(elems[ws])
            rm[w][i] = elems[ws]
    N = len(lst)
    # left multiplication s_i * w : (s_i w)(r) = s_i(w(r))
    lm = [[elems[tuple(gens[i][lst[w][k]] for k in range(len(allroots)))] for i in range(n)] for w in range(N)]
    R = [set(i for i in range(n) if length[rm[w][i]] < length[w]) for w in range(N)]
    L = [set(i for i in range(n) if length[lm[w][i]] < length[w]) for w in range(N)]
    # Bruhat order via subword property on the right: [e,w] = [e,v] u [e,v]s  (w = v s)
    order = sorted(range(N), key=lambda w: length[w])
    below = [None]*N
    below[0] = {0}
    for w in order[1:]:
        s = min(R[w]); v = rm[w][s]
        below[w] = set(below[v]) | {rm[x][s] for x in below[v]}
    # KL polynomials, full table, right recursion
    P = {}
    def padd(a, b, shift=0, coef=1):
        res = list(a) + [0]*max(0, len(b)+shift-len(a))
        for k, c in enumerate(b): res[k+shift] += coef*c
        while res and res[-1] == 0: res.pop()
        return res
    mu = {}
    for w in order:
        P[(w, w)] = [1]
        if w == 0: continue
        s = min(R[w]); v = rm[w][s]
        zs = [z for z in below[v] if z != v and s in R[z] and mu.get((z, v), 0) != 0]
        for x in below[w]:
            if x == w: continue
            xs = rm[x][s]
            c = 1 if length[xs] < length[x] else 0
            res = []
            if xs in below[v]: res = padd(res, P[(xs, v)], 1-c)
            if x in below[v]: res = padd(res, P[(x, v)], c)
            for z in zs:
                if x in below[z]:
                    res = padd(res, P[(x, z)], (length[w]-length[z])//2, -mu[(z, v)])
            P[(x, w)] = res
            d = length[w]-length[x]
            if d % 2 == 1:
                m = res[(d-1)//2] if len(res) > (d-1)//2 else 0
                if m: mu[(x, w)] = m
            assert res[0] == 1 and all(c >= 0 for c in res) and len(res)-1 <= (d-1)//2
    def word(x):
        w = []
        while x != 0:
            s = min(R[x]); w.append(s); x = rm[x][s]
        return ''.join(map(str, reversed(w))) or 'e'
    with open(out, 'w') as f:
        for w in range(N):
            ww = word(w)
            for x in sorted(below[w]):
                f.write(word(x) + ' ' + ww + ' ' + ' '.join(map(str, P[(x, w)])) + '\n')
    print(f"{t}{n}: |W|={N}, pairs={sum(len(b) for b in below)}, max mu={max(mu.values())}, "
          f"max coeff={max(max(p) for p in P.values())}, distinct polys={len(set(tuple(p) for p in P.values()))}")

main()
