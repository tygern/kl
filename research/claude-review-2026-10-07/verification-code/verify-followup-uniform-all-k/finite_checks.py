"""(1) Validate the heap FC test by counting FC elements of E6 (Stembridge: 662).
(2) Check: in E6, E7, E8 every reflection r_beta with beta non-simple is NOT FC,
    and l(r_beta) = 2 dp(beta) - 1 for every positive root beta."""
import sys
from elib import *


def positive_roots(A):
    n = len(A)
    simple = [tuple(1 if i == j else 0 for i in range(n)) for j in range(n)]
    seen = set(simple)
    frontier = list(simple)
    while frontier:
        new = []
        for b in frontier:
            for s in range(n):
                c = simple_refl(A, s, b)
                if is_positive(c) and c not in seen:
                    seen.add(c)
                    new.append(c)
        frontier = new
    return sorted(seen)


def bfs_depth(A, roots):
    """depth by BFS from simple roots along the root graph (independent of greedy)."""
    n = len(A)
    simple = [tuple(1 if i == j else 0 for i in range(n)) for j in range(n)]
    d = {b: 1 for b in simple}
    frontier = list(simple)
    while frontier:
        new = []
        for b in frontier:
            for s in range(n):
                c = simple_refl(A, s, b)
                if is_positive(c) and c not in d:
                    d[c] = d[b] + 1
                    new.append(c)
        frontier = new
    return d


# ---- (1) E6 FC count
A6 = cartan(6)
n = 6
I = identity(n)
seen = {I: []}
frontier = [I]
while frontier:
    new = []
    for M in frontier:
        w = seen[M]
        for s in range(n):
            if not is_negative(M[s]):  # s not a right descent: M s is longer
                N = mult_right_simple(A6, M, s)
                if N not in seen:
                    seen[N] = w + [s]
                    new.append(N)
    frontier = new
print("E6 order:", len(seen))
fc = sum(1 for M, w in seen.items() if is_fc_word(A6, w))
print("E6 FC count (heap test):", fc, "(expected 662)")
# cross-check reduced_word() lengths agree with BFS word lengths on a sample
bad = sum(1 for M, w in list(seen.items())[::97] if len(reduced_word(A6, M)) != len(w))
print("length mismatches on sample:", bad)

# ---- (2) reflections in E6, E7, E8
for n in (6, 7, 8):
    A = cartan(n)
    roots = positive_roots(A)
    bd = bfs_depth(A, roots)
    print(f"E{n}: {len(roots)} positive roots")
    fc_nonsimple = []
    lenviol = []
    depthviol = []
    for b in roots:
        M = reflection_matrix(A, b)
        w = reduced_word(A, M)
        L = len(w)
        dp = depth(A, b)
        if dp != bd[b]:
            depthviol.append(b)
        if L != 2 * dp - 1:
            lenviol.append((b, L, dp))
        simple = sum(b) == 1
        f = is_fc_word(A, w)
        if f and not simple:
            fc_nonsimple.append(b)
        if simple and (not f or L != 1):
            print("ERROR simple reflection", b, f, L)
    print(f"  FC non-simple reflections: {len(fc_nonsimple)}")
    print(f"  violations of l(r_beta)=2dp(beta)-1: {len(lenviol)}")
    print(f"  greedy-depth vs BFS-depth disagreements: {len(depthviol)}")
    print(f"  max reflection length: {max(len(reduced_word(A, reflection_matrix(A, b))) for b in roots)}")
