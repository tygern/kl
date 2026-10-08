import sys, time
from fc_check import adjacency_E, is_fc_reduced, length_by_matrix
from my_family import Z2, Z2z, Z3

def family(n):
    """Uniform family for E_n (any n >= 9), modelled on the enumerated maxima:
    A | Z2(n-2) | Z3(n-4), Z3(n-7), ... (tops decreasing by 3 while top >= 6)
      | tail by residue | final run 2..n-2."""
    p = n-1
    w = list(range(n-2, -1, -1)) + [p]          # A: n letters
    w += Z2(n-2, p)                              # 2n-5 letters
    T = n-4
    while T >= 8:
        w += Z3(T, p); T -= 3
    # now T in {5,6,7}
    if T == 6:   w += Z3(6, p) + Z3(3, p)
    elif T == 5: w += Z2(5, p) + Z3(3, p)
    else:        w += Z2z(7, p) + Z2(5, p) + Z3(3, p)   # T == 7
    w += list(range(2, n-1))                     # final: n-3 letters
    return w

rows = []
for r in range(3, 26):
    n = 4*r+1
    w = family(n)
    adj = adjacency_E(n)
    ok = is_fc_reduced(w, adj)
    L = len(w)
    rows.append((r, n, L, ok))
    print(f"r={r:2d} n={n:3d}: len={L:5d} FC&reduced={ok}  8r^2+2={8*r*r+2:5d}  "
          f"finding 8r^2+7r-2={8*r*r+7*r-2:5d}  L-(8r^2+2)={L-8*r*r-2:4d}  L-(8r^2+2+116)={L-8*r*r-2-116:5d}")
# fit quadratic in r to lengths (check exactness)
import fractions
(r0,_,L0,_),(r1,_,L1,_),(r2,_,L2,_) = rows[0],rows[3],rows[6]   # r=3,6,9 all same residue mod 3
# general: lengths may depend on r mod 3; fit per residue class
for res in range(3):
    pts = [(r,L) for r,n,L,ok in rows if r%3==res][:3]
    (a0,b0),(a1,b1),(a2,b2) = pts
    # solve alpha r^2 + beta r + gamma
    import itertools
    from fractions import Fraction as F
    M = [[F(a*a),F(a),F(1),F(b)] for a,b in pts]
    # gaussian elimination
    for i in range(3):
        piv = M[i][i]
        M[i] = [x/piv for x in M[i]]
        for j in range(3):
            if j!=i:
                f = M[j][i]; M[j] = [x - f*y for x,y in zip(M[j],M[i])]
    al,be,ga = M[0][3],M[1][3],M[2][3]
    pred_ok = all(al*r*r+be*r+ga == L for r,n,L,ok in rows if r%3==res)
    print(f"residue r%3={res}: length = {al} r^2 + {be} r + {ga}; holds on all tested r in this class: {pred_ok}")
# independent exact-matrix reducedness check at n=13, 17, 73
for n in (13, 17, 73, 101):
    w = family(n); t=time.time()
    print(f"matrix length check n={n}: word length {len(w)}, matrix length {length_by_matrix(w, n)}  ({time.time()-t:.1f}s)")
