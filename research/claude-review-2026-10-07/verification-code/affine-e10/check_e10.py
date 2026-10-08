import sys
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/affine-e10')
from cox import *

# E10: chain 0-1-...-8, node 9 attached to node 2 (affine subdiagram on {0..7,9}).
n = 10
edges = [(i, i + 1) for i in range(8)] + [(2, 9)]
A = cartan(n, edges)
gens = [simple_refl(A, s) for s in range(n)]
E = lambda j: tuple(1 if i == j else 0 for i in range(n))

delta = (2, 4, 6, 5, 4, 3, 2, 1, 0, 3)
gamma = (1, 2, 3, 2, 2, 1, 1, 0, 0, 1)
beta0 = (3, 7, 10, 9, 7, 6, 4, 3, 1, 6)
I = {1, 3, 5, 7, 9}

print("(delta,delta) =", pair(A, delta, delta), " (delta,alpha_j) =", [pair(A, delta, E(j)) for j in range(n)])
print("(gamma,gamma) =", pair(A, gamma, gamma), " (gamma,delta) =", pair(A, gamma, delta))
print("(beta0,beta0) =", pair(A, beta0, beta0), " (beta0,gamma) =", pair(A, beta0, gamma), " (beta0,delta) =", pair(A, beta0, delta))


def greedy_lower(v, name):
    """Lower a positive vector by simple reflections with positive pairing; return final vector and word."""
    v = tuple(v)
    word = []
    while True:
        if sum(1 for x in v if x != 0) == 1 and sum(v) == 1:
            return v, word
        for s in range(n):
            if pair(A, v, E(s)) > 0 and apply(gens[s], v) != v:
                w = apply(gens[s], v)
                assert all(x >= 0 for x in w), (name, v, s, w)
                assert sum(w) < sum(v)
                v = w
                word.append(s)
                break
        else:
            raise SystemExit(f"{name}: stuck at {v}")


for name, v in [("gamma", gamma), ("gamma+delta", tuple(g + d for g, d in zip(gamma, delta))), ("beta0", beta0)]:
    end, word = greedy_lower(v, name)
    print(f"{name}: real root, greedy lowering ends at alpha_{end.index(1)} after {len(word)} steps")

# appendix lowering word for beta0
app_word = [1,3,5,4,7,6,5,9,2,1,0,3,2,1,4,3,2,
            5,6,7,8,9,2,1,0,3,2,1,4,3,2,5,4,3,6,
            5,4,7,6,5,9,2,1,0,3,2,1,4,3,2]
v = beta0
ok = True
for s in app_word:
    w = apply(gens[s], v)
    if not all(x >= 0 for x in w) or sum(w) >= sum(v):
        ok = False
        print("  appendix word FAIL at letter", s, v, "->", w)
    v = w
print("appendix lowering word: steps", len(app_word), " end vector", v, " = alpha_9:", v == E(9), " all intermediate nonneg & height decreasing:", ok)

# T and iterates
gd = tuple(g + d for g, d in zip(gamma, delta))
T = compose(reflection(A, gamma), reflection(A, gd))
cur = beta0
print("independence number check & pairings:")
for k in range(0, 6):
    formula = tuple(beta0[i] - k * gamma[i] + (k * k + 3 * k) * delta[i] for i in range(n))
    pairs = tuple(pair(A, E(j), cur) for j in range(n))
    paper_pairs = (-1, 1, -2 - k, 1 + k, -1 - k, 1 + k, -1 - k, 1 + k, -1 - k * k - 3 * k, 2 + k)
    print(f" k={k}: T^k beta0 == formula: {cur == formula}; beta_k={cur}; pairings {pairs} match paper: {pairs == paper_pairs}; (beta_k,beta_k)={pair(A, cur, cur)}; all coords >0: {all(x > 0 for x in cur)}; beta_k0=2k^2+5k+3: {cur[0] == 2*k*k+5*k+3}")
    cur = apply(T, cur)

# descents / terminal / support of reflections b_k for k = 0..4


def terminal_right(cols):
    ok = True
    for s in right_descents(cols):
        for t in range(n):
            if A[s][t] == -1:
                v = tuple(cols[s][i] + cols[t][i] for i in range(n))
                if not is_positive(v):
                    ok = False
    return ok


def support(cols):
    # generator j in support iff some column has coordinate j changed from identity
    return sorted(j for j in range(n) if any(cols[c][j] != (1 if c == j else 0) for c in range(n)))


prev = set()
for k in range(0, 5):
    bk_root = tuple(beta0[i] - k * gamma[i] + (k * k + 3 * k) * delta[i] for i in range(n))
    b = reflection(A, bk_root)
    inv_ok = compose(b, b) == identity(n)
    R = sorted(right_descents(b)); L = sorted(left_descents(A, b))
    term = terminal_right(b) and terminal_right(inverse(b))
    sup = support(b)
    assert b not in prev
    prev.add(b)
    msg = f" b_{k}: involution {inv_ok}; R={R}; L={L}; terminal {term}; support {sup}"
    if k <= 1:
        Lk, rw = length_by_stripping(A, b)
        msg += f"; length {Lk}; reduced-word support {sorted(set(rw))}"
    print(msg)

# independence number of E10 graph by brute force
import itertools
best = 0
for r in range(1, n + 1):
    found = False
    for S in itertools.combinations(range(n), r):
        if all(A[s][t] == 0 for s in S for t in S if s < t):
            found = True
            break
    if found:
        best = r
print("independence number of E10 graph:", best, " I independent:", all(A[s][t] == 0 for s in I for t in I if s < t), " |I| =", len(I))

# (-1)-eigenspace dimensions: reflection has dimension 1 (beta nondegenerate); i(I)
iI = identity(n)
for s in I:
    iI = compose(iI, gens[s])
# dimension of (-1)-eigenspace of iI = rank-nullity of (iI + 1)
from fractions import Fraction
def rank(M):
    M = [[Fraction(x) for x in row] for row in M]
    r = 0
    rows, cols = len(M), len(M[0])
    for c in range(cols):
        p = next((i for i in range(r, rows) if M[i][c] != 0), None)
        if p is None:
            continue
        M[r], M[p] = M[p], M[r]
        for i in range(rows):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        r += 1
    return r
Mplus = [[iI[j][i] + (1 if i == j else 0) for j in range(n)] for i in range(n)]
print("dim (-1)-eigenspace of i(I):", n - rank(Mplus))
b0m = reflection(A, beta0)
Mplus0 = [[b0m[j][i] + (1 if i == j else 0) for j in range(n)] for i in range(n)]
print("dim (-1)-eigenspace of r_beta0:", n - rank(Mplus0))
# Is the E10 form on span(alpha_s: s in I) nondegenerate? (Gram = 2*Id) yes trivially.
