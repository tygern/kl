import sys
sys.path.insert(0, '.')
from en_tools import *
def fmt(r): return ''.join(str(c) for c in r)
cases = {
 'E6 w4 (Gern)': (6, '1325213', [[1,3],[2],[5]]),
 'E7 w6 (Gern)': (7, '132543621324356', [[1,3,5,6],[2,4],[1,3,6]]),
 'E7 new b_7':   (7, '1325436210321432543621324356', [[1,3,5,6],[2,4],[1,3,6],[0,2,4],[1,3,5,6]]),
 'E8 new b_8':   (8, '7534231270123456210321432'+'5437210321432543721324357',
                  [[1,3,5,7],[2,4],[1,3,7],[0,2,4],[1,3,5,7],[0,2,4,6],[1,3,5,7],[2,4]]),
}
for name,(n,wstr,layers) in cases.items():
    rs = RootSystem(n, en_edges(n))
    W = rs.word_matrix(parse(wstr))
    P_layers = layers[:-1]; K = layers[-1]
    Pword = [s for L in P_layers for s in L]
    P = rs.word_matrix(Pword)
    full = Pword + K + Pword[::-1]
    M = rs.word_matrix(full)
    ok = (M == W)
    lens = sum(len(L) for L in layers)*2 - len(K)
    # bipartition classes: colour by parity of distance from node 0 along chain; branch node n-1 adjacent to 2 -> colour of 3
    colour = {j: j % 2 for j in range(n-1)}; colour[n-1] = 1
    classes = [set(colour[s] for s in L) for L in layers]
    alternating = all(len(c)==1 for c in classes) and all(classes[i] != classes[i+1] for i in range(len(classes)-1))
    roots = [rs.apply(P, tuple(int(i==k) for i in range(n))) for k in K]
    neg = [r for r in rs.positive_roots() if rs.apply(W, r) == tuple(-c for c in r)]
    print(f"{name}: W == P i(K) P^-1 as words: {ok}; word length {len(full)} == l(W) {rs.length(W)}: {len(full)==rs.length(W)} (reduced)")
    print(f"   layers {layers}, colour classes {[sorted(c) for c in classes]}, alternating: {alternating}, l(P)={rs.length(P)}")
    print(f"   P(alpha_k) for k in K: {[fmt(r) for r in roots]}  == negated roots {[fmt(r) for r in neg]}: {sorted(roots)==sorted(neg)}")
    print(f"   so W is a product of {len(K)} orthogonal reflections (type {len(K)}A1): {all(rs.pairing(a,b)==0 for a in roots for b in roots if a!=b)}")
