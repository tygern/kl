import sys, json
sys.path.insert(0, ".")
from kl_indep import *
d = json.load(open("/Users/ic/workspace/kl/results/d6-recurrence-certificate.json"))
D6 = SignedD(6)
w6 = (-1, -6, 3, -4, 5, -2); x6 = (-1, -2, 4, 3, 6, 5)
def rw(G, w):
    word=[]
    while G.length(w)>0:
        s=next(s for s in G.gens if G.is_descent(w,s)); word.append(s); w=G.right(w,s)
    return word[::-1]
I = Ideal(D6, rw(D6, w6))
elems = [tuple(e) for e in d["elements"]]
assert set(elems) == set(I.elems), "certificate element set != my lower ideal"
print("certificate elements == my principal lower ideal of w6:", len(elems))
bad = 0; checked = 0; nonzero_strict = 0
for r in d["records"]:
    x, w = elems[r["x"]], elems[r["w"]]
    mine = I.kl(I.idx[x], I.idx[w])
    if list(mine) != list(r["polynomial"]):
        bad += 1
        if bad < 5: print("MISMATCH", x, w, mine, r["polynomial"])
    checked += 1
    if mine not in ((), (1,)): nonzero_strict += 1
print(f"records checked against my implementation: {checked}, mismatches: {bad}, nontrivial polys: {nonzero_strict}")
# also check the certificate's descent convention matches mine (first descent) and correction terms
cdiff = 0
for r in d["records"]:
    if "right_descent" in r:
        x, w = I.idx[elems[r["x"]]], I.idx[elems[r["w"]]]
        s = r["right_descent"]
        v = I.r[w][s]
        mine_terms = sorted((elems_i, e, m) for (z, e, m) in I.corrections(v, s) if I.leq(x, z) for elems_i in [I.elems[z]])
        theirs = sorted((elems[z], e, m) for z, e, m in r["correction_terms"])
        if mine_terms != theirs: cdiff += 1
print("records whose correction-term lists differ from mine:", cdiff)
