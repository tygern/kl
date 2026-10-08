#!/usr/bin/env python3
"""Referee's own independent cross-check of manuscript numbers (exact integer arithmetic).

Written from the manuscript's stated conventions only (columns a_j = w(alpha_j);
right multiplication by s_i: a_i -> -a_i, a_j -> a_j + a_i for j ~ i).
"""
import itertools, json, collections, sys
from fractions import Fraction

def make(n, extra_edges=None):
    edges = [(i, i+1) for i in range(n-2)] + [(2, n-1)]
    adj = [set() for _ in range(n)]
    for s, t in edges:
        adj[s].add(t); adj[t].add(s)
    return edges, adj

def identity(n):
    return tuple(tuple(int(i == j) for i in range(n)) for j in range(n))

def right(a, s, adj):
    cols = list(a)
    cols[s] = tuple(-v for v in a[s])
    for t in adj[s]:
        cols[t] = tuple(u+v for u, v in zip(a[t], a[s]))
    return tuple(cols)

def sign(v):
    assert any(v)
    if all(c >= 0 for c in v): return 1
    if all(c <= 0 for c in v): return -1
    raise ValueError("mixed sign column: not a root %r" % (v,))

def desc(a):
    return frozenset(s for s, c in enumerate(a) if sign(c) < 0)

def element(word, n, adj, check_reduced=True):
    a = identity(n)
    for s in word:
        if check_reduced:
            assert sign(a[s]) > 0, "non-reduced word"
        a = right(a, s, adj)
    return a

def inverse(a, n):
    # inverse of w from its matrix: w^{-1} columns via solving; easier: use word reversal.
    raise NotImplementedError

def reduced_word_by_stripping(a, n, adj, limit=100000):
    """Length of a by repeatedly stripping right descents (independent of inversion formulas)."""
    word = []
    e = identity(n)
    while a != e:
        d = desc(a)
        assert d, "no descent but not identity"
        s = min(d)
        word.append(s)
        a = right(a, s, adj)
        assert len(word) < limit
    return word[::-1]

def terminal_right(a, adj):
    """Manuscript eq (terminaltest): s in R(w) => w(alpha_s+alpha_t) > 0 for all t ~ s."""
    for s in desc(a):
        for t in adj[s]:
            if sign(tuple(u+v for u, v in zip(a[s], a[t]))) < 0:
                return False
    return True

def is_fc_bfs(word, n, adj, edges):
    """FC test via the manuscript's recurrence: R(w) commuting and ws FC for all s in R(w). Memoised on matrices."""
    memo = {identity(n): True}
    def fc(a):
        if a in memo: return memo[a]
        d = desc(a)
        ok = not any(s in d and t in d for s, t in edges) and all(fc(right(a, s, adj)) for s in d)
        memo[a] = ok
        return ok
    return fc(element(word, n, adj))

def independence_number(vertices, edges):
    vs = sorted(vertices)
    best = 0
    for r in range(len(vs), 0, -1):
        for sub in itertools.combinations(vs, r):
            S = set(sub)
            if not any(s in S and t in S for s, t in edges):
                return r
    return 0

report = {}

# ---------------- Table 1 ----------------
TABLE = [
    (6, '1325213', '135', 7, 4),
    (7, '1326213', '136', 7, 4),
    (7, '13256213', '1356', 8, 4),
    (7, '132543621324356', '1356', 15, 11),
    (7, '1325436210321432543621324356', '1356', 28, 24),
    (8, '1327213', '137', 7, 4),
    (8, '13257213', '1357', 8, 4),
    (8, '61327213', '1367', 8, 4),
    (8, '132543721324357', '1357', 15, 11),
    (8, '1325437210321432543721324357', '1357', 28, 24),
    (8, '7534231270123456210321432' + '5437210321432543721324357', '1357', 50, 46),
]
table_ok = True
for n, w, x, ell, gap in TABLE:
    edges, adj = make(n)
    word = list(map(int, w)); xw = list(map(int, x))
    a = element(word, n, adj)            # asserts reducedness along the way
    ai = element(word[::-1], n, adj)
    L = desc(ai); R = desc(a)
    assert len(word) == ell == len(reduced_word_by_stripping(a, n, adj))
    assert R == L == frozenset(xw), (w, sorted(R), sorted(L))
    assert terminal_right(a, adj) and terminal_right(ai, adj), ("not terminal", w)
    assert not is_fc_bfs(word, n, adj, edges), ("FC", w)
    assert ell - len(xw) == gap
    supp = set(word)
    assert independence_number(supp, edges) == len(xw)
    # x is a commuting product below b: subword check (support subset) suffices
    assert set(xw) <= supp and not any(s in xw and t in xw for s, t in edges)
    # b is an involution? (not claimed; just record)
    inv_check = (element(word, n, adj) == element(word[::-1], n, adj))
    print(f"Table1 E{n} len {ell}: OK  descents={sorted(R)} alpha(supp)={len(xw)} involution={inv_check}")
report['table1'] = 'all 11 rows verified: reduced, length, L=R=supp(x), two-sided terminal, non-FC, independence number'

# ---------------- E8 exponents / D7 weight ----------------
n = 8; edges, adj = make(8)
def cartan_pair(x, y, adj):
    n = len(x)
    cy = [2*y[j] - sum(y[t] for t in adj[j]) for j in range(n)]
    return sum(a*b for a, b in zip(x, cy))
beta = (4, 7, 10, 8, 6, 4, 2, 5)
pairs = [cartan_pair(beta, tuple(int(i == j) for i in range(8)), adj) for j in range(8)]
assert pairs == [1, 0, 0, 0, 0, 0, 0, 0], pairs
# orbit under simple reflections
orbit = {beta}; todo = [beta]
while todo:
    v = todo.pop()
    for s in range(8):
        p = cartan_pair(v, tuple(int(i == s) for i in range(8)), adj)
        w = list(v); w[s] -= p; w = tuple(w)
        if w not in orbit:
            orbit.add(w); todo.append(w)
assert len(orbit) == 2160, len(orbit)
print("D7 weight (4,7,10,8,6,4,2,5): pairings", pairs, "orbit size", len(orbit))
# Poincare: |E8|/|D7| = 2160 ; exponents
def poincare(exps):
    p = [1]
    for m in exps:
        q = [0]*(len(p)+m)
        for i, c in enumerate(p):
            for j in range(m+1): q[i+j] += c
        p = q
    return p
E8 = poincare([1,7,11,13,17,19,23,29]); D7 = poincare([1,3,5,7,9,11,6])
assert sum(E8) == 696729600 and sum(D7) == 322560 and sum(E8)//sum(D7) == 2160
assert len(E8)-1 == 120 and len(D7)-1 == 42
print("E8 exponents Poincare sum 696729600, degree 120; D7 322560 degree 42; index 2160")
E7 = poincare([1,5,7,9,11,13,17]); E6 = poincare([1,4,5,7,8,11])
assert sum(E7) == 2903040 and len(E7)-1 == 63 and sum(E6) == 51840 and len(E6)-1 == 36
assert sum(E8)//sum(E7) == 240 and sum(E7)//sum(E6) == 56 and sum(E6)//1920 == 27
print("E7: 2903040 max length 63; E6: 51840 max length 36; indices 240, 56, 27")
# number of roots
def roots(n, adj):
    basis = list(identity(n)); found = set(basis); todo = list(basis)
    while todo:
        a = todo.pop()
        for s in range(n):
            b = list(a); b[s] = -a[s] + sum(a[t] for t in adj[s]); b = tuple(b)
            if b not in found: found.add(b); todo.append(b)
    return found
for nn, expect in ((6, 72), (7, 126), (8, 240)):
    _, ad = make(nn)
    assert len(roots(nn, ad)) == expect
print("root counts 72 / 126 / 240 OK")

# ---------------- Affine E9 family ----------------
n = 9; edges, adj = make(9)
delta = (2,4,6,5,4,3,2,1,3); gamma = (1,1,1,1,1,0,0,0,0); I = {1,3,5,7,8}
b0w = list(map(int, '312875645234123012856745231'))
B0 = element(b0w, n, adj)
B0_manuscript = [
 [0,0,1,-1,1,0,0,0,-1],
 [1,0,2,-2,2,-1,1,-1,-2],
 [2,-1,3,-3,3,-2,2,-1,-2],
 [2,-1,2,-2,3,-2,2,-1,-2],
 [2,-1,1,-1,2,-1,1,-1,-1],
 [2,-1,1,-1,1,0,1,-1,-1],
 [1,-1,1,0,0,0,1,-1,-1],
 [1,-1,1,0,0,0,0,0,-1],
 [1,-1,2,-2,2,-1,1,-1,-1]]
# manuscript matrix has columns b0(alpha_j): entry (i,j) = coordinate i of column j
assert all(B0[j][i] == B0_manuscript[i][j] for i in range(9) for j in range(9)), "B0 matrix mismatch"
print("B_0 matrix matches Appendix B")
assert len(b0w) == 27 == len(reduced_word_by_stripping(B0, n, adj))
assert desc(B0) == I and desc(element(b0w[::-1], n, adj)) == I
assert terminal_right(B0, adj) and terminal_right(element(b0w[::-1], n, adj), adj)
assert not is_fc_bfs(b0w, n, adj, edges)
assert set(b0w) == set(range(9))
# involution: b0^2 = 1 -> matrix squared
def matmul_cols(A, B):
    # (A*B)(alpha_j) = A(B(alpha_j)) ; columns are images in simple-root basis
    n = len(A)
    out = []
    for j in range(n):
        v = [0]*n
        for s in range(n):
            c = B[j][s]
            if c:
                for k in range(n): v[k] += c*A[s][k]
        out.append(tuple(v))
    return tuple(out)
assert matmul_cols(B0, B0) == identity(9)
print("b_0^2 = 1, b_0 terminal both sides, L=R=I, full support, non-FC, length 27")
# delta null, gamma a finite root
assert all(cartan_pair(delta, tuple(int(i == j) for i in range(9)), adj) == 0 for j in range(9))
fin_edges, fin_adj = make(9)
# d_j = (alpha_j, gamma - b0 gamma)
def apply(A, v):
    n = len(A); out = [0]*n
    for s in range(n):
        if v[s]:
            for k in range(n): out[k] += v[s]*A[s][k]
    return tuple(out)
b0g = apply(B0, gamma)
diff = tuple(g - h for g, h in zip(gamma, b0g))
d = tuple(cartan_pair(tuple(int(i == j) for i in range(9)), diff, adj) for j in range(9))
assert d == (2,-1,1,-1,1,-1,1,-1,-1), d
print("d vector OK", d)
# reflection matrices and T = r_gamma r_{gamma+delta}
def reflection(root, adj):
    n = len(root)
    cols = []
    for j in range(n):
        e = tuple(int(i == j) for i in range(n))
        p = cartan_pair(e, root, adj)
        cols.append(tuple(e[k] - p*root[k] for k in range(n)))
    return tuple(cols)
gd = tuple(g+h for g, h in zip(gamma, delta))
T = matmul_cols(reflection(gamma, adj), reflection(gd, adj))
# T(alpha) = alpha - (alpha,gamma) delta
for j in range(9):
    e = tuple(int(i == j) for i in range(9))
    p = cartan_pair(e, gamma, adj)
    assert T[j] == tuple(e[k] - p*delta[k] for k in range(9))
Tinv = matmul_cols(reflection(gd, adj), reflection(gamma, adj))
assert matmul_cols(T, Tinv) == identity(9)
Bk = B0
lengths = []
for k in range(0, 4):
    if k > 0:
        Bk = matmul_cols(matmul_cols(T, Bk), Tinv)
    # eq affinecolumns
    assert all(Bk[j] == tuple(B0[j][i] + k*d[j]*delta[i] for i in range(9)) for j in range(9))
    assert desc(Bk) == I
    assert terminal_right(Bk, adj)
    assert matmul_cols(Bk, Bk) == identity(9)   # involution so left test equals right test
    w = reduced_word_by_stripping(Bk, n, adj)
    lengths.append(len(w))
    assert element(w, n, adj) == Bk
    assert len(w) == 27 + 92*k, (k, len(w))
    assert set(w) == set(range(9))
print("lengths by descent-stripping for k=0..3:", lengths, "= 27+92k OK")
# (A,B) multiplicity table from finite roots
finite_edges, finite_adj = make(8)
R8 = roots(8, finite_adj)  # E8 finite roots in nodes 0..7 of E9? No: finite E8 omits node 7 (affine node) -> nodes {0..6, 8}
# Build finite roots of the E8 parabolic omitting node 7, as 9-vectors with coordinate 7 = 0
fin_nodes = [0,1,2,3,4,5,6,8]
basis9 = [tuple(int(i == j) for i in range(9)) for j in fin_nodes]
found = set(basis9); todo = list(basis9)
while todo:
    a = todo.pop()
    for s in fin_nodes:
        p = cartan_pair(a, tuple(int(i == s) for i in range(9)), adj)
        b = list(a); b[s] -= p; b = tuple(b)
        if b not in found: found.add(b); todo.append(b)
assert len(found) == 240
def eps(v): return 0 if sign(v) > 0 else 1
cnt = collections.Counter()
for eta in found:
    img = apply(B0, eta)
    # img = eta' + h delta with eta' finite (coordinate 7 of eta' is zero)
    h = img[7] // delta[7]
    assert img[7] == h*delta[7]
    etap = tuple(img[i] - h*delta[i] for i in range(9))
    assert etap in found
    A = eps(etap) - h - eps(eta)
    B = -sum(eta[j]*d[j] for j in range(9))
    cnt[(A, B)] += 1
nonzero = {k: v for k, v in cnt.items() if not (k[0] <= 0 and k[1] <= 0)}
print("(A,B) multiplicities with possible positive contribution:", sorted(nonzero.items()))
assert nonzero == {(0,1):44, (0,2):8, (1,0):1, (1,1):20, (1,2):6}, nonzero
zero = sum(v for k, v in cnt.items() if k[0] <= 0 and k[1] <= 0)
assert zero == 161
for k in range(0, 4):
    tot = sum(max(0, A + k*B)*m for (A, B), m in cnt.items())
    assert tot == 27 + 92*k
print("inversion-formula sum equals 27+92k; 161 identically-zero contributions OK")
# 21 covers
covers = {}
for i in range(27):
    sub = b0w[:i] + b0w[i+1:]
    a = element(sub, n, adj, check_reduced=False)
    w = reduced_word_by_stripping(a, n, adj)
    if len(w) == 26:
        covers[a] = w
assert len(covers) == 21, len(covers)
assert all(not is_fc_bfs(w, n, adj, edges) for w in covers.values())
print("21 distinct Bruhat covers of b_0 by letter deletion, none FC: OK")
# independence number of E9 = 5
assert independence_number(set(range(9)), edges) == 5

# ---------------- E10 ----------------
n = 10; edges, adj = make(10)
delta = (2,4,6,5,4,3,2,1,0,3); gamma = (1,2,3,2,2,1,1,0,0,1)
beta0 = (3,7,10,9,7,6,4,3,1,6); I = {1,3,5,7,9}
lowering = [1,3,5,4,7,6,5,9,2,1,0,3,2,1,4,3,2, 5,6,7,8,9,2,1,0,3,2,1,4,3,2,5,4,3,6, 5,4,7,6,5,9,2,1,0,3,2,1,4,3,2]
v = beta0
for s in lowering:
    p = cartan_pair(v, tuple(int(i == s) for i in range(10)), adj)
    w = list(v); w[s] -= p; w = tuple(w)
    assert min(w) >= 0 and sum(w) < sum(v), (s, v, w)
    v = w
assert v == tuple(int(i == 9) for i in range(10)), v
print("E10 lowering word (50 steps) reduces beta_0 to alpha_9 with nonneg decreasing heights: OK")
assert cartan_pair(beta0, gamma, adj) == -3 and cartan_pair(beta0, delta, adj) == -1
assert cartan_pair(beta0, beta0, adj) == 2
ref = reflection(beta0, adj)
rw = reduced_word_by_stripping(ref, n, adj)
assert len(rw) == 101 and element(rw, n, adj) == ref and set(rw) == set(range(10))
assert desc(ref) == I and terminal_right(ref, adj)
print("r_{beta_0} length 101 by descent stripping, full support, R=I, terminal: OK")
for k in range(0, 5):
    bk = tuple(beta0[j] - k*gamma[j] + (k*k+3*k)*delta[j] for j in range(10))
    pr = [cartan_pair(tuple(int(i == j) for i in range(10)), bk, adj) for j in range(10)]
    expect = [-1,1,-2-k,1+k,-1-k,1+k,-1-k,1+k,-1-k*k-3*k,2+k]
    assert pr == expect, (k, pr)
    assert min(bk) > 0 and bk[0] == 2*k*k+5*k+3
    assert cartan_pair(bk, bk, adj) == 2
    rk = reflection(bk, adj)
    assert desc(rk) == I and terminal_right(rk, adj)
    if k <= 2:
        wk = reduced_word_by_stripping(rk, n, adj)
        assert len(wk) % 2 == 1 and set(wk) == set(range(10))
        print(f"  E10 k={k}: length {len(wk)}, pairings OK")
assert independence_number(set(range(10)), edges) == 5
print("E10 family checks OK")

# ---------------- Uniform family r=3 (E13) sanity ----------------
r = 3; n = 4*r+1; edges, adj = make(n)
beta0 = (r-1, r) + tuple(2*r - j//2 for j in range(2, n-1)) + (r,)
delta = (2,4,6,5,4,3,2,1) + (0,)*(n-9) + (3,)
gamma = (1,2,3,2,2,1,1,0) + (0,)*(n-9) + (1,)
I = {0, n-1} | set(range(3, n-1, 2))
assert cartan_pair(beta0, beta0, adj) == 2
ref = reflection(beta0, adj)
assert desc(ref) == I and terminal_right(ref, adj)
w13 = reduced_word_by_stripping(ref, n, adj)
assert set(w13) == set(range(n)) and len(w13) % 2 == 1
assert independence_number(set(range(n)), edges) == len(I) == 7
print(f"E13 (r=3) base reflection: length {len(w13)}, R=I size 7 = independence number, terminal, full support: OK")
print("ALL INDEPENDENT CHECKS PASSED")
