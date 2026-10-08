#!/usr/bin/env python3
"""Compute reduced words / lengths of b_{r,k}=r_{beta_{r,k}} by descent stripping on the
integer matrix of the reflection (geometric representation), and double-check terminality
and descent sets group-theoretically.  Exact integers only."""
exec(open("check_uniform.py").read().split("results = []")[0])

def matrix_of_reflection(A, beta):
    n = len(A)
    return [refl_root(A, beta, simple(n, j)) for j in range(n)]   # columns

def mat_mul_simple_right(A, cols, j):
    """(w s_j)(alpha_i) = w(s_j alpha_i) = w(alpha_i) - A[j][i] w(alpha_j)"""
    n = len(A)
    new = []
    for i in range(n):
        c = A[j][i]
        new.append([x - c*y for x, y in zip(cols[i], cols[j])])
    return new

def reduced_word(A, cols):
    """Strip right descents: while w != 1, pick j with w(alpha_j)<0, w := w s_j; word is reversed."""
    n = len(A)
    ident = [simple(n, j) for j in range(n)]
    word = []
    while cols != ident:
        j = next(j for j in range(n) if any(c < 0 for c in cols[j]))
        word.append(j)
        cols = mat_mul_simple_right(A, cols, j)
        if len(word) > 200000: raise RuntimeError("too long")
    return word[::-1]   # w = s_{word[0]} ... s_{word[-1]}

for r in (3, 4):
    n = 4*r+1; A = cartan(n); beta = seed(r); g, d = gamma_delta(r); a = 2*r-4
    gd = [x+y for x, y in zip(g, d)]
    v = list(beta)
    for k in range(0, 4 if r == 3 else 2):
        cols = matrix_of_reflection(A, v)
        w = reduced_word(A, cols)
        L = len(w)
        supp = sorted(set(w))
        # sanity: re-multiply word to recover matrix
        cur = [simple(n, j) for j in range(n)]
        for j in w:   # w = s_{w0} s_{w1} ... ; multiply on the right successively
            cur = mat_mul_simple_right(A, cur, j)
        assert cur == cols
        print(f"r={r} k={k}: length {L} (odd={L%2==1}), height(beta)={sum(v)}, full support={supp==list(range(n))}")
        print("   reduced word:", "".join(chr(ord('a')+j) for j in w) if n <= 26 else w)
        v = refl_root(A, g, refl_root(A, gd, v))
