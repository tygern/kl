#!/usr/bin/env python3
"""Independent exact-integer verification of Section 4 (rank-uniform family)
of results/exceptional-leading.tex.  Written from scratch; no repo code used.

Numbering: chain 0-1-...-(n-2), node n-1 attached to node 2.  n = 4r+1.
"""
from fractions import Fraction as Fr
import itertools, sys

def cartan(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    for i in range(n-2):          # chain 0..n-2
        A[i][i+1] = A[i+1][i] = -1
    A[n-1][2] = A[2][n-1] = -1    # node n-1 attached to node 2
    return A

def pair(A, x, y):
    n = len(A)
    return sum(x[i]*A[i][j]*y[j] for i in range(n) for j in range(n) if A[i][j] != 0)

def simple(n, j):
    v = [0]*n; v[j] = 1; return v

def refl_simple(A, j, v):
    """s_j v = v - (v, alpha_j) alpha_j"""
    n = len(A)
    c = sum(A[j][i]*v[i] for i in range(n))
    w = list(v); w[j] -= c; return w

def refl_root(A, beta, v):
    """r_beta v = v - (v,beta) beta, assuming (beta,beta)=2"""
    c = pair(A, v, beta)
    return [vi - c*bi for vi, bi in zip(v, beta)]

def det_int(M):
    """Bareiss exact determinant"""
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = 1
    for k in range(n-1):
        if M[k][k] == 0:
            sw = next((i for i in range(k+1, n) if M[i][k] != 0), None)
            if sw is None: return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k+1, n):
            for j in range(k+1, n):
                M[i][j] = (M[i][j]*M[k][k] - M[i][k]*M[k][j]) // prev
        prev = M[k][k]
    return sign*M[n-1][n-1]

def lower_to_simple(A, v, trace=False):
    """Greedy height-reduction: while v is not a simple root, find j with (v,alpha_j)>0,
    apply s_j; demand positivity throughout.  Returns (ok, sequence)."""
    n = len(A); v = list(v); seq = []
    while True:
        if any(c < 0 for c in v): return False, seq
        if sum(v) == 1 and all(c in (0, 1) for c in v):
            return True, seq
        cand = [j for j in range(n) if sum(A[j][i]*v[i] for i in range(n)) > 0]
        if not cand: return False, seq
        j = cand[0]; v = refl_simple(A, j, v); seq.append(j)
        if len(seq) > 10000: return False, seq

def apply_sequence(A, v, seq):
    """Apply s_{seq[0]}, then s_{seq[1]}, ...  Check positivity and strict height decrease."""
    v = list(v); ok = True; hts = [sum(v)]
    for j in seq:
        v = refl_simple(A, j, v)
        if any(c < 0 for c in v): ok = False
        hts.append(sum(v))
    strictly_dec = all(hts[i+1] < hts[i] for i in range(len(hts)-1))
    return v, ok, strictly_dec, hts

def seed(r):
    n = 4*r+1
    b = [0]*n
    b[0], b[1], b[2], b[4*r] = r-1, r, 2*r-1, r
    for j in range(3, 4*r):
        b[j] = 2*r - j//2
    return b

def gamma_delta(r):
    n = 4*r+1
    d = [0]*n; g = [0]*n
    for j, c in enumerate((2, 4, 6, 5, 4, 3, 2, 1)): d[j] = c
    d[4*r] = 3
    for j, c in enumerate((1, 2, 3, 2, 2, 1, 1, 0)): g[j] = c
    g[4*r] = 1
    return g, d

def independence_number(n):
    """Brute force maximum independent set in the E_n graph (n<=25 feasible via
    dynamic programming on the tree; here use simple tree DP)."""
    A = cartan(n)
    adj = {i: [j for j in range(n) if j != i and A[i][j] == -1] for i in range(n)}
    # tree DP rooted at 2
    import functools
    sys.setrecursionlimit(10000)
    def dp(v, parent):
        inc, exc = 1, 0
        for u in adj[v]:
            if u == parent: continue
            i_u, e_u = dp(u, v)
            inc += e_u; exc += max(i_u, e_u)
        return inc, exc
    return max(dp(2, -1))

def is_independent(A, I):
    return all(A[i][j] == 0 for i in I for j in I if i != j)

results = []
def check(label, cond):
    results.append((label, bool(cond)))
    print(("PASS " if cond else "FAIL ") + label)

for r in range(3, 7):
    n = 4*r + 1
    A = cartan(n)
    print(f"\n===== r={r}, n={n} =====")
    beta = seed(r)
    print("beta_r =", beta)
    g, d = gamma_delta(r)
    a = 2*r - 4
    I_r = sorted({0, 4*r} | set(range(3, 4*r, 2)))
    print("I_r =", I_r, "size", len(I_r))

    # determinant / nondegeneracy
    D = det_int(A)
    check(f"r={r}: det Cartan = {D} == 9-n = {9-n}", D == 9-n)
    check(f"r={r}: form nondegenerate", D != 0)

    # (beta_r, alpha_0) = r-2 and (beta,beta)=2
    check(f"r={r}: (beta_r,alpha_0) == r-2", pair(A, beta, simple(n, 0)) == r-2)
    check(f"r={r}: (beta_r,beta_r) == 2", pair(A, beta, beta) == 2)

    # ---- D_{4r} realization check (exact rationals). Coordinates: (z, e_1..e_{4r}) ----
    # inner product: z.z = 2-r, e_i.e_j = delta_ij, z.e_i = 0
    N = 4*r
    def vec(zc, es): return (Fr(zc), [Fr(x) for x in es])
    def ip(u, v): return u[0]*v[0]*(2-r) + sum(x*y for x, y in zip(u[1], v[1]))
    def e(i):   # 1-indexed
        es = [0]*N; es[i-1] = 1; return es
    def add(u, v, cu=1, cv=1): return (cu*u[0]+cv*v[0], [cu*x+cv*y for x, y in zip(u[1], v[1])])
    zero = vec(0, [0]*N)
    alph = {}
    alph[1] = vec(0, [x+y for x, y in zip(e(1), e(2))])
    alph[4*r] = vec(0, [y-x for x, y in zip(e(1), e(2))])
    for j in range(2, 4*r):
        alph[j] = vec(0, [y-x for x, y in zip(e(j), e(j+1))])
    alph[0] = vec(1, [Fr(-1, 2)]*N)
    gram_ok = all(ip(alph[i], alph[j]) == A[i][j] for i in range(n) for j in range(n))
    check(f"r={r}: D_4r realization has E_n Gram matrix", gram_ok)
    # expand beta_r
    B = zero
    for j in range(n):
        B = add(B, alph[j], 1, beta[j])
    exp_c = [Fr(-(r-1), 2) if (i % 2 == 1) else Fr(3-r, 2) for i in range(1, N+1)]
    check(f"r={r}: beta_r = (r-1)z + sum c_i e_i with stated c_i", B[0] == r-1 and B[1] == exp_c)
    # s_0 beta_r = beta_r - (r-2) alpha_0
    s0B = add(B, alph[0], 1, -(r-2))
    Falpha0 = vec(1, [Fr(-1, 2) if (i % 2 == 1) else Fr(1, 2) for i in range(1, N+1)])
    check(f"r={r}: s_0 beta_r == F alpha_0 (flip even coords)", s0B == Falpha0)
    # explicit F as product of r double sign changes r_{e_{4j-2}-e_{4j}} r_{e_{4j-2}+e_{4j}}
    def refl_vec(root, v):
        c = ip(v, root) / ip(root, root) * 2
        return add(v, root, 1, -c)
    v = alph[0]
    for j in range(1, r+1):
        rm = vec(0, [x-y for x, y in zip(e(4*j-2), e(4*j))])
        rp = vec(0, [x+y for x, y in zip(e(4*j-2), e(4*j))])
        v = refl_vec(rp, v); v = refl_vec(rm, v)
    check(f"r={r}: explicit product of r double sign changes sends alpha_0 to F alpha_0", v == Falpha0)
    # and s_0 F alpha_0 == beta_r
    check(f"r={r}: s_0 F alpha_0 == beta_r", add(Falpha0, alph[0], 1, -ip(Falpha0, alph[0])) == B)

    # ---- independent real-root check of beta_r by lowering ----
    ok, seq = lower_to_simple(A, beta)
    check(f"r={r}: beta_r lowers to a simple root (greedy), {len(seq)} steps", ok)

    # ---- gamma, gamma+delta real; delta null in affine parabolic ----
    check(f"r={r}: (gamma,gamma)==2", pair(A, g, g) == 2)
    check(f"r={r}: (delta,delta)==0", pair(A, d, d) == 0)
    check(f"r={r}: (gamma,delta)==0", pair(A, g, d) == 0)
    J = list(range(0, 8)) + [4*r]
    check(f"r={r}: (alpha_j,delta)==0 for j in affine subdiagram", all(pair(A, simple(n, j), d) == 0 for j in J))
    check(f"r={r}: (alpha_8,delta)==-1 (delta not in ambient radical)", pair(A, simple(n, 8), d) == -1)
    gd = [x+y for x, y in zip(g, d)]
    ok_g, _ = lower_to_simple(A, g); ok_gd, _ = lower_to_simple(A, gd)
    check(f"r={r}: gamma real root (greedy lowering)", ok_g)
    check(f"r={r}: gamma+delta real root (greedy lowering)", ok_gd)
    # appendix sequences
    seq_g = [2, 1, 0, 4, 3, 2, 1, 6, 5, 4, 3, 2]
    b = 4*r
    seq_gd = seq_g + [b, 2, 1, 0, 3, 2, 1, 4, 3, 2, 5, 4, 3, 6, 5, 4, 7, 6, 5, b, 2, 1, 0, 3, 2, 1, 4, 3, 2]
    vg, posg, decg, htg = apply_sequence(A, g, seq_g)
    check(f"r={r}: appendix gamma sequence ends at alpha_b, all intermediates >=0, heights strictly decrease",
          vg == simple(n, b) and posg and decg)
    vgd, posgd, decgd, htgd = apply_sequence(A, gd, seq_gd)
    check(f"r={r}: appendix gamma+delta sequence ends at alpha_b, intermediates >=0, heights strictly decrease",
          vgd == simple(n, b) and posgd and decgd)
    print("   gamma+delta heights along sequence:", htgd)

    # ---- pairings of seed with gamma and delta ----
    check(f"r={r}: (beta_r,gamma) == -r", pair(A, beta, g) == -r)
    check(f"r={r}: (beta_r,delta) == -a = {-a}", pair(A, beta, d) == -a)

    # ---- T^k beta_r vs closed form ----
    def T(v): return refl_root(A, g, refl_root(A, gd, v))
    v = list(beta)
    for k in range(0, 8):
        closed = [beta[j] - a*k*g[j] + (a*k*k + r*k)*d[j] for j in range(n)]
        check(f"r={r}: T^{k} beta_r == closed form", v == closed)
        check(f"r={r},k={k}: (beta_rk,beta_rk)==2", pair(A, v, v) == 2)
        check(f"r={r},k={k}: all coordinates strictly positive", all(c > 0 for c in v))
        # pairings m_j
        m = [pair(A, simple(n, j), v) for j in range(n)]
        exp_m = [None]*n
        exp_m[0] = r-2; exp_m[1] = 2-r; exp_m[2] = -1-a*k
        for j in range(3, 8): exp_m[j] = (-1)**(j+1)*(1+a*k)
        exp_m[8] = -1 - a*k*k - r*k
        for j in range(9, 4*r): exp_m[j] = (-1)**(j+1)
        exp_m[4*r] = 1 + a*k
        check(f"r={r},k={k}: pairings m_j match displayed formulas", m == exp_m)
        pos = sorted(j for j in range(n) if m[j] > 0)
        check(f"r={r},k={k}: m_j>0 exactly on I_r", pos == I_r)
        # descent sets via reflection matrix
        cols = [refl_root(A, v, simple(n, j)) for j in range(n)]
        R = sorted(j for j in range(n) if any(c < 0 for c in cols[j]))
        check(f"r={r},k={k}: R(b) == I_r (each column negative <=> all coords negative)",
              R == I_r and all((all(c < 0 for c in cols[j]) if j in I_r else all(c > 0 or (c == 0) for c in cols[j]) and sum(cols[j]) > 0) for j in range(n)))
        # left descents: b is an involution so L=R; check via b^{-1}=b anyway with matrix
        # two-sided terminal test
        term_ok = True
        for s in I_r:
            for t in range(n):
                if t != s and A[s][t] == -1:
                    img = refl_root(A, v, [x+y for x, y in zip(simple(n, s), simple(n, t))])
                    if not (all(c >= 0 for c in img) and any(c > 0 for c in img)): term_ok = False
                    if m[s] + m[t] > 0: term_ok = False
        check(f"r={r},k={k}: right terminal test passes (edge sums <=0, b(alpha_s+alpha_t)>0)", term_ok)
        # involution check: r_beta^2 = id on simple roots
        check(f"r={r},k={k}: r_beta is an involution", all(refl_root(A, v, cols[j]) == simple(n, j) for j in range(n)))
        # full support: every coordinate row changed by column 0
        check(f"r={r},k={k}: column 0 changes every coordinate row", all(cols[0][j] != simple(n, 0)[j] for j in range(n)))
        # distinctness formula
        check(f"r={r},k={k}: beta_rk[0] == r-1+4k+(4r-8)k^2", v[0] == r-1+4*k+(4*r-8)*k*k)
        # lowering certificate for beta_rk (independent real-root check)
        okk, seqk = lower_to_simple(A, v)
        check(f"r={r},k={k}: beta_rk real root by greedy lowering ({len(seqk)} steps, ht={sum(v)})", okk)
        v = T(v)

    # ---- independence number and maximality of I_r ----
    alpha_n = independence_number(n)
    check(f"r={r}: independence number of E_n = {alpha_n} == ceil(n/2) = {(n+1)//2}", alpha_n == (n+1)//2)
    check(f"r={r}: I_r independent and |I_r| == independence number == 2r+1", is_independent(A, I_r) and len(I_r) == alpha_n == 2*r+1)
    # matching stated in the proof
    M = [(0, 1), (2, 4*r)] + [(j, j+1) for j in range(3, 4*r-2, 2)]
    check(f"r={r}: stated matching has 2r edges, all are edges, vertex-disjoint, covers all but one vertex",
          len(M) == 2*r and all(A[i][j] == -1 for i, j in M) and len({x for p in M for x in p}) == 4*r)
    print("   matching:", M)

    # ---- (-1)-eigenspace dimensions (exact rational rank) ----
    def minus_one_eigdim(cols):
        # matrix of (g + I); nullity = dim of (-1)-eigenspace
        Mx = [[Fr(cols[j][i]) + (1 if i == j else 0) for j in range(n)] for i in range(n)]
        # rank via Gaussian elimination
        rank = 0; rows = Mx
        for c in range(n):
            piv = next((i for i in range(rank, n) if rows[i][c] != 0), None)
            if piv is None: continue
            rows[rank], rows[piv] = rows[piv], rows[rank]
            pv = rows[rank][c]
            rows[rank] = [x/pv for x in rows[rank]]
            for i in range(n):
                if i != rank and rows[i][c] != 0:
                    f = rows[i][c]; rows[i] = [x - f*y for x, y in zip(rows[i], rows[rank])]
            rank += 1
        return n - rank
    colsb = [refl_root(A, beta, simple(n, j)) for j in range(n)]
    # i(I_r) matrix
    def iI(vv):
        for s in I_r: vv = refl_simple(A, s, vv)
        return vv
    colsI = [iI(simple(n, j)) for j in range(n)]
    check(f"r={r}: (-1)-eigenspace dim of r_beta == 1", minus_one_eigdim(colsb) == 1)
    check(f"r={r}: (-1)-eigenspace dim of i(I_r) == 2r+1", minus_one_eigdim(colsI) == 2*r+1)

# ---- r=2 remark ----
print("\n===== r=2 remark (n=9) =====")
r = 2; n = 9; A = cartan(n); beta = seed(2); g, d = gamma_delta(2)
print("beta_2 =", beta, " delta-gamma =", [x-y for x, y in zip(d, g)])
check("r=2: (alpha_0,beta_2)==0", pair(A, beta, simple(n, 0)) == 0)
check("r=2: beta_2 == delta - gamma", beta == [x-y for x, y in zip(d, g)])
check("r=2: beta_2 is a real root (norm 2, lowering)", pair(A, beta, beta) == 2 and lower_to_simple(A, beta)[0])
m = [pair(A, simple(n, j), beta) for j in range(n)]
print("r=2 pairings:", m, " positive on", [j for j in range(n) if m[j] > 0], " (I_2 would be {0,3,5,7,8})")
check("r=2: det Cartan == 0 (affine)", det_int(A) == 0)

# ---- independence number formula for general n ----
for n in range(6, 30):
    an = independence_number(n)
    if an != (n+1)//2: print("INDEPENDENCE FORMULA FAILS at n=", n, an)
print("independence number formula ceil(n/2) verified for 6<=n<30")

nfail = sum(1 for _, ok in results if not ok)
print(f"\nTOTAL checks: {len(results)}, failures: {nfail}")
