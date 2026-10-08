#!/usr/bin/env python3
"""Presentation-dimension sanity checks for exceptional-leading.tex (exact integer arithmetic)."""
from fractions import Fraction
from itertools import combinations
import sys

def en_edges(n):
    # generalized E_n: chain 0-1-...-(n-2), node n-1 attached to node 2
    E = [(i, i + 1) for i in range(n - 2)] + [(2, n - 1)]
    return E

def cartan(n, edges=None):
    if edges is None:
        edges = en_edges(n)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    for (i, j) in edges:
        A[i][j] = A[j][i] = -1
    return A

def adj(n, edges):
    nb = [set() for _ in range(n)]
    for i, j in edges:
        nb[i].add(j); nb[j].add(i)
    return nb

def count_independent_sets(n, edges):
    nb = adj(n, edges)
    cnt = 0
    for mask in range(1 << n):
        ok = True
        for i in range(n):
            if mask >> i & 1:
                for j in nb[i]:
                    if mask >> j & 1:
                        ok = False; break
            if not ok: break
        if ok: cnt += 1
    return cnt

def independence_number(n, edges):
    nb = adj(n, edges)
    best = 0
    for mask in range(1 << n):
        ok = True
        for i in range(n):
            if mask >> i & 1 and any(mask >> j & 1 for j in nb[i]):
                ok = False; break
        if ok: best = max(best, bin(mask).count('1'))
    return best

# geometric representation: w stored as matrix with columns w(alpha_j) in simple-root coords
def identity(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

def right_mult(M, i, A):
    """M * s_i : columns a_j -> a_j + a_i (j~i), a_i -> -a_i."""
    n = len(M)
    col_i = [M[r][i] for r in range(n)]
    N = [row[:] for row in M]
    for j in range(n):
        if j == i:
            for r in range(n): N[r][j] = -col_i[r]
        elif A[i][j] == -1:
            for r in range(n): N[r][j] = M[r][j] + col_i[r]
    return N

def word_to_matrix(word, n, A):
    M = identity(n)
    for c in word:
        M = right_mult(M, c, A)
    return M

def is_neg(col):
    return all(c <= 0 for c in col) and any(c < 0 for c in col)

def is_pos(col):
    return all(c >= 0 for c in col) and any(c > 0 for c in col)

def right_descents(M):
    n = len(M)
    return {j for j in range(n) if is_neg([M[r][j] for r in range(n)])}

def inverse(M):
    # for Weyl group matrices in root coordinates, inverse via adjugate is overkill; use Fraction Gauss.
    n = len(M)
    aug = [[Fraction(M[r][c]) for c in range(n)] + [Fraction(1 if r == c else 0) for c in range(n)] for r in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if aug[r][c] != 0)
        aug[c], aug[p] = aug[p], aug[c]
        pv = aug[c][c]
        aug[c] = [x / pv for x in aug[c]]
        for r in range(n):
            if r != c and aug[r][c] != 0:
                f = aug[r][c]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[c])]
    return [[int(aug[r][n + c]) for c in range(n)] for r in range(n)]

def length_by_descent_stripping(M, A):
    n = len(M)
    L = 0
    M = [row[:] for row in M]
    while True:
        D = right_descents(M)
        if not D: return L
        s = min(D)
        M = right_mult(M, s, A); L += 1

def terminal_right(M, A):
    n = len(M)
    D = right_descents(M)
    for s in D:
        for t in range(n):
            if A[s][t] == -1:
                col = [M[r][s] + M[r][t] for r in range(n)]
                if not is_pos(col): return False
    return True

def check_table():
    rows = [
        (6, "1325213", "135"),
        (7, "1326213", "136"),
        (7, "13256213", "1356"),
        (7, "132543621324356", "1356"),
        (7, "1325436210321432543621324356", "1356"),
        (8, "1327213", "137"),
        (8, "13257213", "1357"),
        (8, "61327213", "1367"),
        (8, "132543721324357", "1357"),
        (8, "1325437210321432543721324357", "1357"),
        (8, "7534231270123456210321432" + "5437210321432543721324357", "1357"),
    ]
    for n, w, x in rows:
        A = cartan(n)
        word = [int(c) for c in w]
        M = word_to_matrix(word, n, A)
        ln = length_by_descent_stripping(M, A)
        Minv = inverse(M)
        R = right_descents(M); L = right_descents(Minv)
        invol = (M == Minv)
        term = terminal_right(M, A) and terminal_right(Minv, A)
        X = set(int(c) for c in x)
        print(f"E{n} word={w[:12]}... len(word)={len(word)} reduced_len={ln} reduced={ln==len(word)} "
              f"R={sorted(R)} L={sorted(L)} matches_x={R==X and L==X} involution={invol} terminal={term}")

def check_affine():
    n = 9
    A = cartan(9)  # chain 0..7, node 8 attached to 2
    w = "312875645234123012856745231"
    word = [int(c) for c in w]
    M = word_to_matrix(word, n, A)
    Minv = inverse(M)
    print("affine b0: len(word)=", len(word), "involution=", M == Minv,
          "R=", sorted(right_descents(M)), "terminal=", terminal_right(M, A) and terminal_right(Minv, A))
    # reduced? compute length via stripping with affine matrices (positive roots have nonneg coords)
    print("affine b0 reduced length by stripping:", length_by_descent_stripping(M, A))
    B0 = [[0,0,1,-1,1,0,0,0,-1],[1,0,2,-2,2,-1,1,-1,-2],[2,-1,3,-3,3,-2,2,-1,-2],[2,-1,2,-2,3,-2,2,-1,-2],
          [2,-1,1,-1,2,-1,1,-1,-1],[2,-1,1,-1,1,0,1,-1,-1],[1,-1,1,0,0,0,1,-1,-1],[1,-1,1,0,0,0,0,0,-1],[1,-1,2,-2,2,-1,1,-1,-1]]
    print("B0 matches appendix:", M == B0)
    delta = [2,4,6,5,4,3,2,1,3]
    gamma = [1,1,1,1,1,0,0,0,0]
    def pair(u, v):
        return sum(u[i]*A[i][j]*v[j] for i in range(n) for j in range(n))
    b0gamma = [sum(M[r][c]*gamma[c] for c in range(n)) for r in range(n)]
    diff = [gamma[i]-b0gamma[i] for i in range(n)]
    d = [pair([1 if i==j else 0 for i in range(n)], diff) for j in range(n)]
    print("d =", d, " expected (2,-1,1,-1,1,-1,1,-1,-1)")
    print("(delta,delta)=", pair(delta,delta), "(gamma,gamma)=", pair(gamma,gamma), "(gamma,delta)=", pair(gamma,delta))
    # covers: delete one letter, keep length-26 products, count distinct
    seen = {}
    for i in range(len(word)):
        sub = word[:i] + word[i+1:]
        Ms = word_to_matrix(sub, n, A)
        if length_by_descent_stripping(Ms, A) == 26:
            seen[tuple(tuple(r) for r in Ms)] = sub
    print("distinct covers from letter deletion:", len(seen))

def lower(vec, seq, A):
    n = len(A)
    v = vec[:]
    hts = [sum(v)]
    for s in seq:
        p = sum(A[s][j]*v[j] for j in range(n))
        v[s] -= p
        if any(c < 0 for c in v):
            return None, hts
        hts.append(sum(v))
    return v, hts

def check_lowering():
    # E10 seed
    n = 10; A = cartan(10)
    beta0 = [3,7,10,9,7,6,4,3,1,6]
    seq = [1,3,5,4,7,6,5,9,2,1,0,3,2,1,4,3,2, 5,6,7,8,9,2,1,0,3,2,1,4,3,2,5,4,3,6, 5,4,7,6,5,9,2,1,0,3,2,1,4,3,2]
    v, hts = lower(beta0, seq, A)
    print("E10 lowering: steps=", len(seq), "final=", v, "heights strictly decreasing=", all(a > b for a, b in zip(hts, hts[1:])), "start ht", hts[0], "end ht", hts[-1])
    # gamma in E_{4r+1} for r=3 (n=13), branch b=12
    for r in (3, 4):
        n = 4*r+1; A = cartan(n); b = 4*r
        gamma = [0]*n; delta = [0]*n
        for j, g in enumerate([1,2,3,2,2,1,1,0]): gamma[j] = g
        for j, dd in enumerate([2,4,6,5,4,3,2,1]): delta[j] = dd
        gamma[b] = 1; delta[b] = 3
        seq1 = [2,1,0,4,3,2,1,6,5,4,3,2]
        v, hts = lower(gamma, seq1, A)
        alpha_b = [1 if i == b else 0 for i in range(n)]
        print(f"r={r}: gamma lowered to alpha_b:", v == alpha_b)
        seq2 = seq1 + [b,2,1,0,3,2,1,4,3,2,5,4,3,6,5,4,7,6,5, b,2,1,0,3,2,1,4,3,2]
        gd = [gamma[i]+delta[i] for i in range(n)]
        v, hts = lower(gd, seq2, A)
        print(f"r={r}: gamma+delta lowered to alpha_b:", v == alpha_b, "heights strictly decreasing:", all(a > bb for a, bb in zip(hts, hts[1:])))
    # E10 version with b=9
    n = 10; A = cartan(10); b = 9
    gamma = [1,2,3,2,2,1,1,0,0,1]; delta = [2,4,6,5,4,3,2,1,0,3]
    seq1 = [2,1,0,4,3,2,1,6,5,4,3,2]
    alpha_b = [1 if i == b else 0 for i in range(n)]
    print("E10 gamma lowered:", lower(gamma, seq1, A)[0] == alpha_b)
    seq2 = seq1 + [b,2,1,0,3,2,1,4,3,2,5,4,3,6,5,4,7,6,5, b,2,1,0,3,2,1,4,3,2]
    print("E10 gamma+delta lowered:", lower([gamma[i]+delta[i] for i in range(n)], seq2, A)[0] == alpha_b)

def check_seed_r():
    for r in (3, 4, 5):
        n = 4*r+1
        beta = [0]*n
        beta[0], beta[1], beta[2], beta[4*r] = r-1, r, 2*r-1, r
        for j in range(3, 4*r):
            beta[j] = 2*r - j//2
        if r == 3: print("beta_3 =", beta, " manuscript: (2,3,5,5,4,4,3,3,2,2,1,1,3)")
        A = cartan(n)
        m = [sum(A[j][i]*beta[i] for i in range(n)) for j in range(n)]
        print(f"r={r} pairings m_j:", m)

# ---- D4 KL polynomial P_{x, xcx} ----
def kl_d4():
    # D4: center 2 (index), leaves 0,1,3 ; use edges (0,2),(1,2),(2,3)
    n = 4
    edges = [(0,2),(1,2),(3,2)]
    A = cartan(n, edges)
    # enumerate group by BFS on matrices
    from collections import deque
    I = identity(n)
    key = lambda M: tuple(tuple(r) for r in M)
    elems = {key(I): (I, 0)}
    dq = deque([I])
    while dq:
        M = dq.popleft()
        for s in range(n):
            N = right_mult(M, s, A)
            k = key(N)
            if k not in elems:
                elems[k] = (N, elems[key(M)][1] + 1)
                dq.append(N)
    assert len(elems) == 192
    order = sorted(elems.values(), key=lambda t: t[1])
    idx = {key(M): i for i, (M, l) in enumerate(order)}
    mats = [M for M, l in order]
    lens = [l for M, l in order]
    N = len(mats)
    def rmul(i, s):
        return idx[key(right_mult(mats[i], s, A))]
    def rdesc(i):
        return right_descents(mats[i])
    # Bruhat order via subword / standard recursion: x <= w iff (s in R(w)): if s in R(x): xs<=ws else x<=ws
    from functools import lru_cache
    @lru_cache(None)
    def leq(x, w):
        if x == w: return True
        if lens[x] >= lens[w]: return False
        D = rdesc(w)
        s = min(D)
        ws = rmul(w, s)
        if s in rdesc(x):
            return leq(rmul(x, s), ws)
        return leq(x, ws)
    # KL polynomials via standard recursion P_{x,w} with w = v s, s in R(w)
    # P_{x,w} = q^{1-c} P_{xs,v} + q^c P_{x,v} - sum_{z: zs<z, x<=z<v} mu(z,v) q^{(l(w)-l(z))/2} P_{x,z}, c=1 if xs<x else 0
    # polynomials as tuple of ints
    def padd(p, q):
        m = max(len(p), len(q)); return tuple((p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(m))
    def pshift(p, k): return tuple([0]*k + list(p))
    def pscale(p, c): return tuple(c*a for a in p)
    def ptrim(p):
        p = list(p)
        while p and p[-1] == 0: p.pop()
        return tuple(p)
    P = {}
    def getP(x, w):
        if not leq(x, w): return ()
        if x == w: return (1,)
        if (x, w) in P: return P[(x, w)]
        s = min(rdesc(w)); v = rmul(w, s)
        xs = rmul(x, s)
        c = 1 if lens[xs] < lens[x] else 0
        res = padd(pshift(getP(xs, v), 1 - c), pshift(getP(x, v), c))
        for z in range(N):
            if lens[z] < lens[v] and s in rdesc(z) and leq(x, z) and leq(z, v):
                m = mu(z, v)
                if m:
                    res = padd(res, pscale(pshift(getP(x, z), (lens[w] - lens[z]) // 2), -m))
        res = ptrim(res)
        P[(x, w)] = res
        return res
    def mu(z, v):
        d = lens[v] - lens[z]
        if d % 2 == 0 or not leq(z, v): return 0
        p = getP(z, v)
        k = (d - 1) // 2
        return p[k] if k < len(p) else 0
    x = idx[key(word_to_matrix([0,1,3], n, A))]
    b = idx[key(word_to_matrix([0,1,3,2,0,1,3], n, A))]
    print("D4: l(x)=", lens[x], "l(b)=", lens[b], "P_{x,b}=", getP(x, b), "(expected 1+2q)")
    # also check all mu(x,w) with x FC in D4 are <=1 ? skip

if __name__ == "__main__":
    for n in (6, 7, 8):
        print(f"E{n}: #independent sets =", count_independent_sets(n, en_edges(n)), " alpha =", independence_number(n, en_edges(n)))
    check_table()
    check_affine()
    check_lowering()
    check_seed_r()
    kl_d4()
