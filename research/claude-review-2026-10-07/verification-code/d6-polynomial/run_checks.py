#!/usr/bin/env python3
import sys, time, json
sys.path.insert(0, "/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/d6-polynomial")
from kl_indep import *

report = {}

# ---------------------------------------------------------------- self tests
print("== Self-test: type A3 (S4) via geometric model")
A3 = Geometric([(0, 1), (1, 2)], 3)
w0 = [0, 1, 0, 2, 1, 0]
I = Ideal(A3, w0, "A3")
assert I.N == 24
e = I.idx[A3.identity()]
# known: the only w in S4 with P_{e,w} != 1 are 3412 = s2 s1 s3 s2 and 4231 = s1 s2 s3 s2 s1 (1-based)
w3412 = I.idx[word_product(A3, [1, 0, 2, 1])]
w4231 = I.idx[word_product(A3, [0, 1, 2, 1, 0])]
nontriv = {}
for w in range(24):
    for x in range(24):
        p = I.kl(x, w)
        if p not in ((), (1,)):
            nontriv.setdefault(w, []).append((x, p))
assert set(nontriv) == {w3412, w4231}, nontriv
s2 = I.idx[word_product(A3, [1])]
assert sorted(nontriv[w3412]) == sorted([(e, (1, 1)), (s2, (1, 1))])
s1 = I.idx[word_product(A3, [0])]; s3 = I.idx[word_product(A3, [2])]; s13 = I.idx[word_product(A3, [0, 2])]
assert sorted(nontriv[w4231]) == sorted([(e, (1, 1)), (s1, (1, 1)), (s3, (1, 1)), (s13, (1, 1))])
# Poincare polynomial of S4
assert Counter(I.length).__eq__(Counter({0: 1, 1: 3, 2: 5, 3: 6, 4: 5, 5: 3, 6: 1}))
print("   S4 singular Schubert varieties 3412, 4231 recovered with P=1+q on the expected x; all others P=1. OK")

print("== Self-test: D4 Mongelli example P_{x,xcx} = 1+2q (x = product of leaves, c = centre)")
D4 = SignedD(4)   # code labels: 0,1 leaves (Gern s1,s2), 2 centre (s3), 3 leaf (s4)
Ib = Ideal(D4, [0, 1, 3, 2, 0, 1, 3], "D4")
x = Ib.idx[word_product(D4, [0, 1, 3])]
print("   P =", show(Ib.kl(x, Ib.top)))
assert Ib.kl(x, Ib.top) == (1, 2)
# same in the geometric D4 model
D4g = Geometric([(0, 2), (1, 2), (2, 3)], 4)
Ig = Ideal(D4g, [0, 1, 3, 2, 0, 1, 3], "D4g")
assert Ig.kl(Ig.idx[word_product(D4g, [0, 1, 3])], Ig.top) == (1, 2)
print("   OK (both models)")

# full D4 group: compare signed-permutation length formula with BFS length
w0D4 = [0, 1, 2, 0, 1, 2, 3, 2, 0, 1, 2, 3]
ID4 = Ideal(D4, w0D4, "D4 full")
assert ID4.N == 192
for i, w in enumerate(ID4.elems):
    assert D4.length(w) == ID4.length[i]
print("   D4: inversion-count length formula agrees with BFS length on all 192 elements. OK")

# ------------------------------------------------------ main object: D6 pair
print("== D6, signed permutations, Gern's conventions")
D6 = SignedD(6)
x6 = (-1, -2, 4, 3, 6, 5)
w6 = (-1, -6, 3, -4, 5, -2)
assert word_product(D6, [0, 1, 3, 5]) == x6      # Gern x_6 = s1 s2 s4 s6
print("   x6 = s1 s2 s4 s6 =", x6, " length", D6.length(x6))
print("   w6 =", w6, " length", D6.length(w6))
# Gern Thm 2.2.18 / Lemma 2.2.20 for n = 6: a1 = (-1)^{3} = -1, a_i = 2i-1, b_i = -(8-2i)
assert w6 == (-1, -6, 3, -4, 5, -2) == tuple([(-1)**3] + [v for i in range(1, 4) for v in ((-(8 - 2 * i),) if i == 1 else (2 * i - 1, -(8 - 2 * i)))])
assert sum(a < 0 for a in w6) % 2 == 0 and sum(a < 0 for a in x6) % 2 == 0
# get a reduced word of w6 by stripping descents
def reduced_word_of(G, w):
    word = []
    while G.length(w) > 0:
        s = next(s for s in G.gens if G.is_descent(w, s))
        word.append(s); w = G.right(w, s)
    return word[::-1]
rw6 = reduced_word_of(D6, w6)
print("   a reduced word of w6 (code labels):", rw6, " (Gern labels):", [s + 1 for s in rw6])
t0 = time.time()
I6 = Ideal(D6, rw6, "D6 w6")
e6 = I6.idx[D6.identity()]
xi, wi = I6.idx[x6], I6.top
print(f"   principal lower ideal |[e,w6]| = {I6.N}  (built in {I6.build_time:.1f}s)")
for i, w in enumerate(I6.elems):
    assert D6.length(w) == I6.length[i]
print("   BFS length agrees with inversion formula on the whole ideal")
assert I6.leq(xi, wi)
iv = I6.interval(xi, wi)
rv = I6.rank_vector(xi, wi)
print(f"   |[x6,w6]| = {len(iv)}  rank vector = {rv}")
fc, nwords = is_fully_commutative(I6, xi)
print(f"   x6 fully commutative: {fc}  (#reduced words {nwords})")
fcw, nw = is_fully_commutative(I6, wi)
print(f"   w6 fully commutative: {fcw}  (#reduced words {nw})")
P = I6.kl(xi, wi)
print(f"   P_{{x6,w6}} = {show(P)}   coefficients {P}   ({time.time()-t0:.1f}s, {len(I6.P)} memoised pairs)")
gap = I6.length[wi] - I6.length[xi]
mu = P[(gap - 1) // 2] if len(P) > (gap - 1) // 2 else 0
print(f"   gap = {gap}, mu(x6,w6) = {mu}")
descs = [s for s in D6.gens if I6.desc[wi][s]]
print("   right descents of w6 (code labels):", descs, "(Gern labels):", [s + 1 for s in descs])
for s in descs:
    assert I6.kl_with_descent(xi, wi, s) == P
print("   recursion with every right descent of w6 gives the same P")
Pe = I6.kl(e6, wi)
print(f"   P_{{e,w6}} = {show(Pe)}   (Gern Lemma 2.3.9 predicts P_{{x6,w6}} = P_{{e,w6}}: {Pe == P})")
print("   L(x6) = R(x6) in Gern labels:", [s + 1 for s in D6.gens if I6.desc[xi][s]])
# left descents of w6: descents of inverse; w6 is an involution
inv = tuple(sorted(range(1, 7), key=lambda k: abs(w6[k - 1])))  # not needed; w6 involution
assert all(w6[abs(w6[k]) - 1] * (1 if w6[k] > 0 else -1) == k + 1 for k in range(6)), "w6 involution"
report["D6"] = dict(x6=x6, w6=w6, len_x6=D6.length(x6), len_w6=D6.length(w6), ideal=I6.N,
                    interval=len(iv), rank_vector=rv, P=list(P), mu=mu, x6_FC=fc, P_e_w6=list(Pe),
                    reduced_word_w6_gern=[s + 1 for s in rw6])

# -------------------------------------------- D6 geometric model cross-check
print("== D6, geometric model (Gern diagram: s1,s2 attached to s3, chain s3-s4-s5-s6)")
D6g = Geometric([(0, 2), (1, 2), (2, 3), (3, 4), (4, 5)], 6)
I6g = Ideal(D6g, rw6, "D6g")
xg = I6g.idx[word_product(D6g, [0, 1, 3, 5])]
Pg = I6g.kl(xg, I6g.top)
print(f"   |ideal| = {I6g.N}, |[x,w]| = {len(I6g.interval(xg, I6g.top))}, rank vector {I6g.rank_vector(xg, I6g.top)}")
print(f"   P = {show(Pg)}  equal to signed model: {Pg == P}")
# compare the whole set of KL polynomials computed in both models via element bijection by reduced word
def canon_word(ideal, w):
    word = []
    while ideal.length[w] > 0:
        s = ideal.first_descent(w); word.append(s); w = ideal.r[w][s]
    return tuple(word[::-1])
map_s2g = {}
for i in range(I6.N):
    map_s2g[i] = I6g.idx[word_product(D6g, canon_word(I6, i))]
assert len(set(map_s2g.values())) == I6.N == I6g.N
# Bruhat order agreement
for i in range(I6.N):
    assert I6.low[i] == sum(1 << map_s2g[k] for k in I6.members(I6.low[i])) if False else True
mismatch = 0
for i in range(I6.N):
    mi = map_s2g[i]
    for k in I6.members(I6.low[i]):
        if not I6g.leq(map_s2g[k], mi):
            mismatch += 1
    if bin(I6.low[i]).count("1") != bin(I6g.low[mi]).count("1"):
        mismatch += 1
print(f"   Bruhat order agreement between the two D6 models: {'OK' if mismatch == 0 else mismatch}")
for (a, b), p in I6.P.items():
    assert I6g.kl(map_s2g[a], map_s2g[b]) == p
print(f"   all {len(I6.P)} memoised KL polynomials agree between the two D6 models")

# --------------------------------------------- E7 / E8 Table 1 length-15 rows
print("== E7 and E8: Table 1 length-15 rows")
def E_edges(n):  # chain 0-1-...-(n-2), node n-1 attached to node 2
    return [(i, i + 1) for i in range(n - 2)] + [(2, n - 1)]
for n, bword, xword in [(7, "132543621324356", "1356"), (8, "132543721324357", "1357")]:
    G = Geometric(E_edges(n), n)
    b = [int(ch) for ch in bword]
    xw = [int(ch) for ch in xword]
    I = Ideal(G, b, f"E{n}")
    xb = I.idx[word_product(G, xw)]
    Pb = I.kl(xb, I.top)
    print(f"   E{n}: b={bword} len {I.length[I.top]}, x={xword} len {I.length[xb]}, |ideal|={I.N}, |[x,b]|={len(I.interval(xb, I.top))}")
    print(f"        rank vector {I.rank_vector(xb, I.top)}")
    print(f"        P_(x,b) = {show(Pb)}   equals D6 value: {Pb == P}")
    fcx, _ = is_fully_commutative(I, xb)
    L_b = [s for s in G.gens if I.desc[I.top][s]]
    print(f"        x FC: {fcx}; R(b) = {L_b}; R(x) = {[s for s in G.gens if I.desc[xb][s]]}; support(b) = {sorted(set(b))}")
    # label mapping: Gern labels 1..6 -> nodes 1, n-1, 2, 3, 4, 5  (node -> Gern)
    node_to_gern = {1: 1, n - 1: 2, 2: 3, 3: 4, 4: 5, 5: 6}
    assert set(b) <= set(node_to_gern) and set(xw) <= set(node_to_gern)
    # check the induced subdiagram is Gern's D6 diagram under this relabelling
    sub_edges = {frozenset((node_to_gern[a], node_to_gern[c])) for a, c in E_edges(n) if a in node_to_gern and c in node_to_gern}
    assert sub_edges == {frozenset(p) for p in [(1, 3), (2, 3), (3, 4), (4, 5), (5, 6)]}, sub_edges
    bg = [node_to_gern[s] - 1 for s in b]     # code labels
    xg_ = [node_to_gern[s] - 1 for s in xw]
    wb = word_product(D6, bg); xx = word_product(D6, xg_)
    print(f"        mapped word (Gern labels) {[c+1 for c in bg]} -> signed perm {wb}  == w6: {wb == w6}")
    print(f"        mapped x    (Gern labels) {[c+1 for c in xg_]} -> signed perm {xx}  == x6: {xx == x6}")
    wb_rev = word_product(D6, bg[::-1])
    print(f"        reversed word gives {wb_rev} (w6 is an involution: {wb_rev == w6})")
    # also check the D6-parabolic element bijection: ideal elements of b are in W_J, J = supp(b)
    J = set(node_to_gern)
    assert all(all(all(c == (1 if k == j else 0) for k, c in enumerate(col) if k not in J) for j, col in enumerate(w)) for w in I.elems), 'ideal not inside W_J'
    report[f"E{n}"] = dict(b=bword, x=xword, len_b=I.length[I.top], len_x=I.length[xb], ideal=I.N,
                           interval=len(I.interval(xb, I.top)), rank_vector=I.rank_vector(xb, I.top),
                           P=list(Pb), mapped_b=wb, mapped_x=xx, x_FC=fcx)

json.dump(report, open("/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/d6-polynomial/report.json", "w"), indent=1, default=list)
print("done")
