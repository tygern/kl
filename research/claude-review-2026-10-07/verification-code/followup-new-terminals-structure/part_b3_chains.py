import sys
from collections import deque
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-new-terminals-structure')
from en_tools import *

E7W = '1325436210321432543621324356'
E8W = '7534231270123456210321432' + '5437210321432543721324357'

def fmt(r): return ''.join(str(c) for c in r)

def max_w0_factor(rs, W):
    suffixes = {W}
    q = deque([W]); best = 0; bestJ = None
    while q:
        V = q.popleft()
        Lv = tuple(sorted(rs.left_descents(V)))
        if Lv:
            l = rs.length(rs.longest_element(list(Lv)))
            if l > best: best, bestJ = l, Lv
        for s in rs.left_descents(V):
            Y = rs.left_mult(s, V)
            if Y not in suffixes:
                suffixes.add(Y); q.append(Y)
    return best, bestJ, len(suffixes)

def peel(rs, W, name):
    print(f"peeling {name}: l={rs.length(W)}")
    cur = W
    while True:
        L = sorted(rs.left_descents(cur)); R = sorted(rs.right_descents(cur))
        if not L: break
        indepL = not any(t in L for s in L for t in rs.adj[s])
        indepR = not any(t in R for s in R for t in rs.adj[s])
        print(f"   l={rs.length(cur):2d} L={L} R={R} commutingL={indepL} commutingR={indepR} involution={rs.mult(cur,cur)==rs.identity()}"
              f" word={word_str(rs.reduced_word(cur))}")
        if not (indepL and indepR) or L != R: break
        iI = rs.word_matrix(L)
        nxt = rs.mult(rs.mult(iI, cur), iI)
        if rs.length(nxt) != rs.length(cur) - 2 * len(L):
            print(f"   (double coset not free: l(i d i) = {rs.length(nxt)})"); break
        cur = nxt

for n, wstr in [(7, E7W), (8, E8W)]:
    rs = RootSystem(n, en_edges(n))
    W = rs.word_matrix(parse(wstr))
    print(f"\n==== E_{n} ====")
    fixed = [r for r in rs.positive_roots() if rs.apply(W, r) == r]
    print("fixed positive roots of new element:", [fmt(r) for r in fixed], "count", len(fixed))
    G6 = rs.word_matrix([gern_to_en(s, n) for s in gern_wn_word(6)])
    print("negated roots of Gern w6:", [fmt(r) for r in rs.positive_roots() if rs.apply(G6, r) == tuple(-c for c in r)])
    b, J, cnt = max_w0_factor(rs, G6)
    print(f"Gern w6 in E_{n}: max l(w_0(L(V))) over suffixes = {b} (J={J}); weak interval size {cnt}; a(w6)=5 per Gern, so must be <= 5: {b <= 5}")
    G4 = rs.word_matrix([gern_to_en(s, n) for s in gern_wn_word(4)])
    b4, J4, _ = max_w0_factor(rs, G4)
    print(f"Gern w4 in E_{n}: max l(w_0(L(V))) = {b4} (J={J4}); a(w4)=3 per Gern")
    peel(rs, W, f"E{n} new element")
    peel(rs, G6, "Gern w6")
    # sum of negated roots divisible by 2 (D4-frame test) for E7
    neg = [r for r in rs.positive_roots() if rs.apply(W, r) == tuple(-c for c in r)]
    S = [sum(c) for c in zip(*neg)]
    print("sum of negated roots:", S, "all even:", all(c % 2 == 0 for c in S))

# Gern chain in D10 (signed permutations): is w_{n} a prefix of w_{n+2}?  Use inversion sets in D10.
print("\n==== Gern's own chain in D_10 (right weak order = prefix order) ====")
m = 10
def image_root(w, root):
    out = {}
    for i, c in root.items():
        j = w[i - 1]; out[abs(j)] = out.get(abs(j), 0) + (c if j > 0 else -c)
    return out
pos_roots = []
for i in range(1, m + 1):
    for j in range(i + 1, m + 1):
        pos_roots.append({j: 1, i: -1}); pos_roots.append({i: 1, j: 1})
def root_is_positive(r):
    items = [(k, v) for k, v in r.items() if v != 0]
    (a, ca), (b, cb) = items
    if ca > 0 and cb > 0: return True
    if ca < 0 and cb < 0: return False
    hi = a if ca > 0 else b; lo = b if ca > 0 else a
    return hi > lo
def inv(w):
    out = [0] * m
    for i, j in enumerate(w, 1): out[abs(j) - 1] = i if j > 0 else -i
    return tuple(out)
def N(w):  # inversion set of w: positive roots sent negative
    return {tuple(sorted(r.items())) for r in pos_roots if not root_is_positive(image_root(w, r))}
ws = {k: signed_perm_of_gern_word(gern_wn_word(k), m) for k in (4, 6, 8, 10)}
for a, b in [(4, 6), (6, 8), (8, 10), (4, 8)]:
    # a prefix of b  iff N(a^{-1}) subset N(b^{-1}); a suffix iff N(a) subset N(b)
    print(f"w{a} prefix of w{b}: {N(inv(ws[a])) <= N(inv(ws[b]))}; suffix: {N(ws[a]) <= N(ws[b])}; lengths {len(N(ws[a]))},{len(N(ws[b]))}")
