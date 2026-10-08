# Is the affine base element b_0 (Theorem affine) a reflection?  (-1)-eigenspace dimension = 9 - rank(B0 + I), exact rationals.
from fractions import Fraction
rows = [[0,0,1,-1,1,0,0,0,-1],[1,0,2,-2,2,-1,1,-1,-2],[2,-1,3,-3,3,-2,2,-1,-2],[2,-1,2,-2,3,-2,2,-1,-2],
        [2,-1,1,-1,2,-1,1,-1,-1],[2,-1,1,-1,1,0,1,-1,-1],[1,-1,1,0,0,0,1,-1,-1],[1,-1,1,0,0,0,0,0,-1],[1,-1,2,-2,2,-1,1,-1,-1]]
n = 9
def rank(M):
    M = [[Fraction(x) for x in r] for r in M]; rk = 0
    for c in range(n):
        p = next((i for i in range(rk, n) if M[i][c] != 0), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(n):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]; M[i] = [a - f*b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk
B = rows
BpI = [[B[i][j] + (1 if i == j else 0) for j in range(n)] for i in range(n)]
B2 = [[sum(B[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
print("B0^2 = I:", all(B2[i][j] == (1 if i == j else 0) for i in range(n) for j in range(n)))
print("dim (-1)-eigenspace of b_0 =", n - rank(BpI))
