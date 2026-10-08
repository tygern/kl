#!/usr/bin/env python3
"""Exact-integer check: are the printed words in Table 1 / Section 5 involutions,
reduced, and do their left/right descent sets equal the stated support?
Geometric representation in simple-root coordinates, W acting on the LEFT:
matrix of w has column j = w(alpha_j)."""
from fractions import Fraction

def edges(n):
    # chain 0-1-...-(n-2), node n-1 attached to node 2
    E = set()
    for i in range(n-2):
        E.add((i, i+1)); E.add((i+1, i))
    E.add((2, n-1)); E.add((n-1, 2))
    return E

def cartan(n):
    E = edges(n)
    return [[2 if i == j else (-1 if (i, j) in E else 0) for j in range(n)] for i in range(n)]

def ident(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]

def refl(n, A, i):
    # s_i(alpha_j) = alpha_j - A[i][j] alpha_i  (A symmetric here); column j
    M = ident(n)
    for j in range(n):
        M[i][j] -= A[i][j]
    return M

def mul(X, Y):
    n = len(X)
    return [[sum(X[i][k]*Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

def word_matrix(n, A, word, reverse=False):
    letters = [int(c) for c in word]
    if reverse:
        letters = letters[::-1]
    M = ident(n)
    for i in letters:
        M = mul(M, refl(n, A, i))   # right-multiply: product in the order written
    return M

def neg_col(M, j):
    col = [M[i][j] for i in range(len(M))]
    assert all(c <= 0 for c in col) or all(c >= 0 for c in col), ("mixed sign col", col)
    return all(c <= 0 for c in col)

def right_desc(M):
    return {j for j in range(len(M)) if neg_col(M, j)}

def length(n, A, M):
    # strip right descents
    M = [row[:] for row in M]
    L = 0
    while True:
        R = right_desc(M)
        if not R:
            return L
        i = min(R)
        M = mul(M, refl(n, A, i)); L += 1

def inverse(n, A, M):
    # M is an integer matrix of a Coxeter group element; invert via Fractions
    import copy
    aug = [[Fraction(M[i][j]) for j in range(n)] + [Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if aug[r][c] != 0)
        aug[c], aug[p] = aug[p], aug[c]
        pv = aug[c][c]
        aug[c] = [x / pv for x in aug[c]]
        for r in range(n):
            if r != c and aug[r][c] != 0:
                f = aug[r][c]
                aug[r] = [a - f*b for a, b in zip(aug[r], aug[c])]
    inv = [[aug[i][n+j] for j in range(n)] for i in range(n)]
    assert all(x.denominator == 1 for row in inv for x in row)
    return [[int(x) for x in row] for row in inv]

table = [
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
    (8, "7534231270123456210321432" + "5437210321432543721324357", "1357"),  # b_{8,6}
    (9, "312875645234123012856745231", "13578"),  # affine b_0
]

for n, w, x in table:
    A = cartan(n)
    M = word_matrix(n, A, w)
    Mr = word_matrix(n, A, w, reverse=True)
    I = ident(n)
    invol = (mul(M, M) == I)
    same = (M == Mr)
    ln = length(n, A, M)
    reduced = (ln == len(w))
    R = right_desc(M)
    Lset = right_desc(inverse(n, A, M))
    supp = set(int(c) for c in x)
    print(f"E{n} word={w[:12]+'...' if len(w)>14 else w:>14} len={len(w):2d} "
          f"reduced={reduced} involution={invol} fwd==rev={same} "
          f"R={sorted(R)} L={sorted(Lset)} R==L==supp(x)={R==Lset==supp}")
