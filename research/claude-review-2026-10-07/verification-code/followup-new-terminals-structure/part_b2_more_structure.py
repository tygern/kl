"""More structure: relation to Gern's w6, weak-order chain, i(I) d i(I) form, explicit w_0(J) factorisation."""
import sys
from collections import deque
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-new-terminals-structure')
from en_tools import *

E7W = '1325436210321432543621324356'
E8W = '7534231270123456210321432' + '5437210321432543721324357'
E7_IN_E8 = E7W.replace('6', '7')

def reflection_matrix(rs, beta):
    n = rs.n
    cols = []
    for j in range(n):
        e = tuple(int(i == j) for i in range(n))
        c = rs.pairing(e, beta)
        cols.append(tuple(a - c * b for a, b in zip(e, beta)))
    return tuple(cols)

def negated_roots(rs, M):
    return [r for r in rs.positive_roots() if rs.apply(M, r) == tuple(-c for c in r)]

def fmt(r):
    return ''.join(str(c) for c in r)

for n, wstr, gern_list, prev in [
        (7, E7W, {'w4': gern_wn_word(4), 'w6': gern_wn_word(6)}, None),
        (8, E8W, {'w4': gern_wn_word(4), 'w6': gern_wn_word(6)}, E7_IN_E8)]:
    rs = RootSystem(n, en_edges(n))
    W = rs.word_matrix(parse(wstr))
    I = sorted(rs.left_descents(W))
    print(f"\n==== E_{n} new element, I = {I} ====")
    negW = negated_roots(rs, W)
    print("negated roots of new element:", [fmt(r) for r in negW])
    for name, gw in gern_list.items():
        G = rs.word_matrix([gern_to_en(s, n) for s in gw])
        negG = negated_roots(rs, G)
        print(f"negated roots of Gern {name} (E_n labels): {[fmt(r) for r in negG]}  subset of new element's: {set(negG) <= set(negW)}")
        extra = [r for r in negW if r not in negG]
        P = G
        for r in extra:
            P = rs.mult(P, reflection_matrix(rs, r))
        print(f"   new element == Gern {name} * prod of reflections in {[fmt(r) for r in extra]} : {P == W}")
        # is Gern element a prefix (right weak order) of new element?
        Gi, Wi = rs.inverse(G), rs.inverse(W)
        pref = set(rs.inversions(Gi)) <= set(rs.inversions(Wi))
        suff = set(rs.inversions(G)) <= set(rs.inversions(W))
        print(f"   Gern {name} is a prefix of new element: {pref}; a suffix: {suff}")
        if pref:
            Y = rs.mult(Gi, W)
            print(f"   quotient y = {name}^-1 w: length {rs.length(Y)}, word {word_str(rs.reduced_word(Y))}, "
                  f"L(y)={sorted(rs.left_descents(Y))}, R(y)={sorted(rs.right_descents(Y))}, involution: {rs.mult(Y,Y)==rs.identity()}")
    if prev:
        Pm = rs.word_matrix(parse(prev))
        negP = negated_roots(rs, Pm)
        print(f"negated roots of E7 new element (relabelled in E8): {[fmt(r) for r in negP]}")
        Pi, Wi = rs.inverse(Pm), rs.inverse(W)
        pref = set(rs.inversions(Pi)) <= set(rs.inversions(Wi))
        print(f"   E7 new element is a prefix of E8 new element: {pref}")
        Z = rs.mult(Pi, W)
        print(f"   quotient z = w_E7^-1 w_E8: length {rs.length(Z)}, word {word_str(rs.reduced_word(Z))}, "
              f"L(z)={sorted(rs.left_descents(Z))}, R(z)={sorted(rs.right_descents(Z))}, involution: {rs.mult(Z,Z)==rs.identity()}")
        # eigen-structure relation
        print(f"   pairings of E8 element's negated roots with E7 element's negated roots:")
        for a in negW:
            print("     ", fmt(a), [rs.pairing(a, b) for b in negP])
    # i(I) d i(I)
    iI = rs.word_matrix(I)
    D = rs.mult(rs.mult(iI, W), iI)
    print(f"d = i(I) w i(I): length {rs.length(D)} (= l(w) - 2|I| = {rs.length(W) - 2*len(I)}: {rs.length(D) == rs.length(W) - 2*len(I)}), "
          f"word {word_str(rs.reduced_word(D))}, L(d)={sorted(rs.left_descents(D))}, R(d)={sorted(rs.right_descents(D))}, "
          f"terminal: {rs.terminal(D)}, negated roots: {[fmt(r) for r in negated_roots(rs, D)]}")
    # explicit factorisation w = u w_0(J) v with J = I + {0}
    J = sorted(set(I) | {0})
    w0J = rs.longest_element(J)
    print(f"J = {J}, l(w_0(J)) = {rs.length(w0J)}, word {word_str(rs.reduced_word(w0J))}")
    # BFS suffixes
    suffixes = {W: rs.identity()}  # V -> u with W = u V
    q = deque([W])
    hit = None
    while q:
        V = q.popleft()
        if set(J) <= rs.left_descents(V):
            hit = V
            break
        for s in rs.left_descents(V):
            Y = rs.left_mult(s, V)
            if Y not in suffixes:
                suffixes[Y] = rs.mult(suffixes[V], rs.word_matrix([s]))
                q.append(Y)
    if hit is not None:
        u = suffixes[hit]
        v = rs.mult(rs.inverse(w0J), hit)
        assert rs.mult(rs.mult(u, w0J), v) == W
        print(f"w = u * w_0(J) * v with l(u)={rs.length(u)}, l(v)={rs.length(v)}, additive: {rs.length(u)+rs.length(w0J)+rs.length(v)==rs.length(W)}")
        print(f"   u = {word_str(rs.reduced_word(u))}, w_0(J) = {word_str(rs.reduced_word(w0J))}, v = {word_str(rs.reduced_word(v))}")
    # Is there any factor w_0(J') with J' containing an A3 (length >= 6 connected)? Check all suffixes for L(V) containing a path of 3 nodes
    paths = []
    for V in suffixes:
        Lv = rs.left_descents(V)
        for a in Lv:
            for b in rs.adj[a] & Lv:
                for c in (rs.adj[b] & Lv) - {a}:
                    paths.append((a, b, c))
    print(f"   any suffix whose left descent set contains a 3-path (=> w_0(A3) factor): {bool(paths)}")
    # Also: max over all suffixes of l(w_0(L(V))) computed again, and which J achieve it (set)
    best = {}
    for V in suffixes:
        Lv = tuple(sorted(rs.left_descents(V)))
        if Lv:
            best[Lv] = rs.length(rs.longest_element(list(Lv)))
    mx = max(best.values())
    print(f"   max l(w_0(L(V))) over suffixes V: {mx}; achieved by J in {[J_ for J_, l in best.items() if l == mx]}")
    # Length of the inherited pair check: Gern's a(w6)=5 means bound (15-5)/2 = 5 = (15-4-1)/2 -> no vanishing (consistent with mu=1)
