import sys, itertools
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/affine-e10')
from cox import *

# E9 = affine E8 numbering from paper: chain 0-1-...-7, node 8 attached to 2.
n = 9
edges = [(i, i + 1) for i in range(7)] + [(2, 8)]
A = cartan(n, edges)
gens = [simple_refl(A, s) for s in range(n)]

delta = (2, 4, 6, 5, 4, 3, 2, 1, 3)
gamma = (1, 1, 1, 1, 1, 0, 0, 0, 0)
I = {1, 3, 5, 7, 8}
word = [int(c) for c in "312875645234123012856745231"]
print("word length", len(word))

# Paper's B_0 (rows as printed; columns = b0(alpha_j))
B0_rows = [
    [0, 0, 1, -1, 1, 0, 0, 0, -1],
    [1, 0, 2, -2, 2, -1, 1, -1, -2],
    [2, -1, 3, -3, 3, -2, 2, -1, -2],
    [2, -1, 2, -2, 3, -2, 2, -1, -2],
    [2, -1, 1, -1, 2, -1, 1, -1, -1],
    [2, -1, 1, -1, 1, 0, 1, -1, -1],
    [1, -1, 1, 0, 0, 0, 1, -1, -1],
    [1, -1, 1, 0, 0, 0, 0, 0, -1],
    [1, -1, 2, -2, 2, -1, 1, -1, -1],
]
B0_paper = tuple(tuple(B0_rows[i][j] for i in range(n)) for j in range(n))

# delta null root?
print("A*delta =", [sum(A[i][j] * delta[j] for j in range(n)) for i in range(n)])
print("(gamma,gamma) =", pair(A, gamma, gamma), " (gamma,delta) =", pair(A, gamma, delta))
# gamma a root: lower greedily
v = list(gamma)
while sum(v) > 1:
    for s in range(n):
        if pair(A, v, tuple(1 if i == s else 0 for i in range(n))) > 0:
            v = list(apply(gens[s], tuple(v)))
            break
    else:
        raise SystemExit("gamma not lowerable")
print("gamma lowers to", v)

ltr = word_product(A, word, True)
rtl = word_product(A, word, False)
print("left-to-right product equals B0_paper:", ltr == B0_paper)
print("right-to-left product equals B0_paper:", rtl == B0_paper)
b0 = ltr if ltr == B0_paper else rtl
if ltr != B0_paper and rtl != B0_paper:
    print("NEITHER matches!")
    print("ltr cols:", ltr)
    print("rtl cols:", rtl)
print("b0 is involution:", compose(b0, b0) == identity(n))
print("ltr involution:", compose(ltr, ltr) == identity(n), " rtl involution:", compose(rtl, rtl) == identity(n))
print("ltr == rtl:", ltr == rtl, " ltr == inverse(rtl):", ltr == inverse(rtl))
print("word reduced (ltr):", is_reduced_word(A, word), " (rtl):", is_reduced_word(A, list(reversed(word))))
L0, rw = length_by_stripping(A, b0)
print("length by stripping:", L0)
print("R(b0) =", sorted(right_descents(b0)), " L(b0) =", sorted(left_descents(A, b0)))
print("full support of word:", sorted(set(word)) == list(range(n)))

# Terminal test
def terminal_right(cols):
    ok = True
    for s in right_descents(cols):
        for t in range(n):
            if A[s][t] == -1:
                v = tuple(cols[s][i] + cols[t][i] for i in range(n))
                if not is_positive(v):
                    ok = False
                    print("  right-terminal FAIL at", s, t, v)
    return ok
print("two-sided terminal:", terminal_right(b0) and terminal_right(inverse(b0)))

# d vector
bg = apply(b0, gamma)
dvec = tuple(pair(A, tuple(1 if i == j else 0 for i in range(n)), tuple(gamma[i] - bg[i] for i in range(n))) for j in range(n))
print("d =", dvec, " paper d = (2,-1,1,-1,1,-1,1,-1,-1)", dvec == (2, -1, 1, -1, 1, -1, 1, -1, -1))

# T = r_gamma r_{gamma+delta}
gd = tuple(gamma[i] + delta[i] for i in range(n))
T = compose(reflection(A, gamma), reflection(A, gd))
Tinv = inverse(T)
# check T(alpha) = alpha - (alpha,gamma) delta
for j in range(n):
    e = tuple(1 if i == j else 0 for i in range(n))
    expect = tuple(e[i] - pair(A, e, gamma) * delta[i] for i in range(n))
    assert T[j] == expect
print("T formula verified on simple roots")

def conj(k):
    w = b0
    for _ in range(k):
        w = compose(T, compose(w, Tinv))
    return w

lengths = {}
for k in range(0, 4):
    bk = conj(k)
    pred = tuple(tuple(b0[j][i] + k * dvec[j] * delta[i] for i in range(n)) for j in range(n))
    inv_ok = compose(bk, bk) == identity(n)
    Lk, _ = length_by_stripping(A, bk)
    lengths[k] = Lk
    print(f"k={k}: eq:affinecolumns holds: {pred == bk}; involution: {inv_ok}; R={sorted(right_descents(bk))}; L={sorted(left_descents(A, bk))}; terminal={terminal_right(bk) and terminal_right(inverse(bk))}; length={Lk}; 27+92k={27+92*k}")
    # full support via stripping word
    _, rwk = length_by_stripping(A, bk)
    print("   support of a reduced word:", sorted(set(rwk)))

# Finite E8 roots (nodes != 7): enumerate via orbit of simple roots under finite parabolic
fin_nodes = [i for i in range(n) if i != 7]
def finite_roots():
    roots = set()
    frontier = [tuple(1 if i == s else 0 for i in range(n)) for s in fin_nodes]
    roots.update(frontier)
    while frontier:
        new = []
        for r in frontier:
            for s in fin_nodes:
                im = apply(gens[s], r)
                if im not in roots:
                    roots.add(im)
                    new.append(im)
        frontier = new
    return roots
R = finite_roots()
print("number of finite roots:", len(R))

# inversion count certificate
from collections import Counter
def eps(v):
    return 1 if is_negative(v) else 0
table = Counter()
total_zero = 0
def decompose(v):
    # v = eta' + h delta with eta' finite (coordinate 7 = 0)
    h = v[7]  # since delta_7 = 1
    etap = tuple(v[i] - h * delta[i] for i in range(n))
    assert etap[7] == 0
    return etap, h
for eta in R:
    v = apply(b0, eta)
    etap, h = decompose(v)
    assert etap in R
    Aeta = eps(etap) - h - eps(eta)
    Beta = -sum(eta[j] * dvec[j] for j in range(n))
    if Aeta <= 0 and Beta <= 0:
        total_zero += 1
    else:
        table[(Aeta, Beta)] += 1
print("nonzero-contribution (A,B) multiplicities:", dict(sorted(table.items())))
print("identically zero contributions:", total_zero)
for k in range(0, 4):
    s = sum(m * max(0, a + k * b) for (a, b), m in table.items())
    print(f"  inversion-count length k={k}: {s}  (stripping gave {lengths[k]})")
# sanity: also check the zero ones are really nonpositive for all k>=0
for eta in R:
    v = apply(b0, eta)
    etap, h = decompose(v)
    Aeta = eps(etap) - h - eps(eta)
    Beta = -sum(eta[j] * dvec[j] for j in range(n))
    if (Aeta, Beta) not in table:
        assert Aeta <= 0 and Beta <= 0

# Bruhat covers of b0 by letter deletion
covers = set()
for i in range(len(word)):
    w2 = word[:i] + word[i + 1:]
    el = word_product(A, w2, True)
    Lw, _ = length_by_stripping(A, el)
    if Lw == 26:
        covers.add(el)
print("distinct covers of length 26 from single deletions:", len(covers))
import pickle
with open('/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/affine-e10/covers.pkl', 'wb') as f:
    pickle.dump(sorted(covers), f)
