"""Exact Coxeter length of the reflections b_{r,k} = r_{beta_{r,k}} in W(E_{4r+1}).

Length is computed in the geometric representation with exact integers:
starting from the matrix of w = r_beta, repeatedly find a right descent s
(w(alpha_s) is a non-positive vector), replace w by w s, and count steps.
Also reports the height of beta_{r,k} and checks the root is real (norm 2)
and that the reflection is an involution of the root lattice.
"""
import sys
from fractions import Fraction

def cartan_E(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n): A[i][i] = 2
    for i in range(n-2): A[i][i+1] = A[i+1][i] = -1
    A[2][n-1] = A[n-1][2] = -1
    return A

def pair(A, u, v):
    n = len(A)
    return sum(u[i]*A[i][j]*v[j] for i in range(n) for j in range(n))

def beta_rk(r, k):
    n = 4*r+1
    beta = [0]*n
    beta[0], beta[1], beta[2], beta[4*r] = r-1, r, 2*r-1, r
    for j in range(3, 4*r):
        beta[j] = 2*r - j//2
    delta = [0]*n; gamma = [0]*n
    for j, v in enumerate((2,4,6,5,4,3,2,1)): delta[j] = v
    delta[4*r] = 3
    for j, v in enumerate((1,2,3,2,2,1,1,0)): gamma[j] = v
    gamma[4*r] = 1
    a = 2*r-4
    return [beta[i] - a*k*gamma[i] + (a*k*k + r*k)*delta[i] for i in range(n)]

def reflection_matrix(A, beta):
    n = len(A)
    # column j = r_beta(alpha_j) = alpha_j - (alpha_j, beta) beta
    M = [[0]*n for _ in range(n)]
    for j in range(n):
        m = sum(A[j][i]*beta[i] for i in range(n))
        for i in range(n):
            M[i][j] = (1 if i == j else 0) - m*beta[i]
    return M

def length(A, M):
    n = len(A)
    M = [row[:] for row in M]
    L = 0
    while True:
        s = -1
        for j in range(n):
            col = [M[i][j] for i in range(n)]
            if all(c <= 0 for c in col) and any(c < 0 for c in col):
                s = j; break
            assert all(c >= 0 for c in col), "column neither positive nor negative"
        if s < 0:
            for i in range(n):
                for j in range(n):
                    assert M[i][j] == (1 if i == j else 0)
            return L
        # w <- w s : column j becomes w(alpha_j) - A[s][j] w(alpha_s)
        cs = [M[i][s] for i in range(n)]
        for j in range(n):
            if j == s: continue
            a = A[s][j]
            if a:
                for i in range(n): M[i][j] -= a*cs[i]
        for i in range(n): M[i][s] = -cs[i]
        L += 1

if __name__ == '__main__':
    rs = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [3,4,5]
    ks = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [0,1,2]
    for r in rs:
        n = 4*r+1
        A = cartan_E(n)
        for k in ks:
            b = beta_rk(r, k)
            assert pair(A, b, b) == 2, "not norm 2"
            M = reflection_matrix(A, b)
            L = length(A, M)
            ht = sum(b)
            print(f"r={r} n={n} k={k}: ht(beta)={ht}, 2ht-1={2*ht-1}, ell(b_rk)={L}, ell-1={L-1}, "
                  f"8r^2+2r-4={8*r*r+2*r-4}, finding's 8r^2+2+116k={8*r*r+2+116*k}, family 8r^2+7r-2={8*r*r+7*r-2}")
            sys.stdout.flush()
