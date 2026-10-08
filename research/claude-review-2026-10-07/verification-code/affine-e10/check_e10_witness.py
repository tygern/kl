import sys
sys.path.insert(0, '.')
from cox import *
n = 10
edges = [(i, i + 1) for i in range(8)] + [(2, 9)]
A = cartan(n, edges)
gens = [simple_refl(A, s) for s in range(n)]
E = lambda j: tuple(1 if i == j else 0 for i in range(n))
delta = (2, 4, 6, 5, 4, 3, 2, 1, 0, 3)
gamma = (1, 2, 3, 2, 2, 1, 1, 0, 0, 1)
b = 9
w1 = [2,1,0,4,3,2,1,6,5,4,3,2]
w2 = w1 + [b,2,1,0,3,2,1,4,3,2,5,4,3,6,5,4,7,6,5, b,2,1,0,3,2,1,4,3,2]
def run(v, word):
    ok = True
    for s in word:
        w = apply(gens[s], v)
        if not all(x >= 0 for x in w) or sum(w) >= sum(v):
            ok = False
        v = w
    return v, ok
print("gamma witness:", run(gamma, w1), "ends at alpha_9:", run(gamma, w1)[0] == E(9))
gd = tuple(g + d for g, d in zip(gamma, delta))
print("gamma+delta witness:", run(gd, w2), "ends at alpha_9:", run(gd, w2)[0] == E(9))
# Also check the same witnesses in E9 numbering?  (the E9 section uses a different gamma; skip)
