"""Independent E8 two-sided terminal classification by one-step parabolic pruning.

Lemma: if w = a v (a in W^J minimal left-coset reps, v in W_J, lengths add) is
right-terminal in W, then v is right-terminal in W_J.  Take J = E7 = {0..5,7}.
The 576 right-terminals of E7 come from my brute-force C++ enumeration
(e7-right-terminals.txt, labels 0..6 with 6 the branch node) and are relabelled
6 -> 7 into E8.  W^J (240 elements) is enumerated by left-weak-order BFS.
"""
import sys, time
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-new-terminals-structure')
from en_tools import *

n = 8
rs = RootSystem(n, en_edges(n))
J = [0, 1, 2, 3, 4, 5, 7]
t0 = time.time()
# W^J: elements with no right descent in J.  BFS by left multiplication.
I = rs.identity()
WJ = {I}
frontier = [I]
while frontier:
    new = []
    for a in frontier:
        for s in range(n):
            b = rs.left_mult(s, a)
            if b in WJ:
                continue
            if rs.length(b) < rs.length(a):
                continue
            if all(not rs.is_negative(b[t]) for t in J):
                WJ.add(b)
                new.append(b)
    frontier = new
print(f"|W^J| = {len(WJ)} (expected 240)")
assert len(WJ) == 240

rt7 = [line.strip() for line in open('e7-right-terminals.txt')]
assert len(rt7) == 576
V = []
for w in rt7:
    word = [7 if c == '6' else int(c) for c in w]
    V.append(rs.word_matrix(word))
# sanity: each v right-terminal in E8 restricted to J (same test, J-generators only)
def right_terminal_in(M, K):
    for s in K:
        if rs.is_negative(M[s]):
            for t in rs.adj[s] & set(K):
                if not rs.is_positive(tuple(a + b for a, b in zip(M[s], M[t]))):
                    return False
    return True
assert all(right_terminal_in(v, J) for v in V)

cands = 0
right_terms = []
for a in WJ:
    for v in V:
        w = rs.mult(a, v)
        cands += 1
        if rs.right_terminal(w):
            right_terms.append(w)
print(f"candidates {cands}, right-terminals in E8: {len(right_terms)} (paper: 2160)  [{time.time()-t0:.1f}s]")
two_sided = [w for w in right_terms if rs.right_terminal(rs.inverse(w))]
print(f"two-sided terminals: {len(two_sided)} (paper: 64)")
noncomm = [w for w in two_sided if not rs.is_commuting_product(w)]
print(f"noncommuting: {len(noncomm)} (paper: 6); commuting: {len(two_sided)-len(noncomm)} (paper: 58)")
TABLE8 = ['1327213', '13257213', '61327213', '132543721324357', '1325437210321432543721324357',
          '7534231270123456210321432' + '5437210321432543721324357']
tab = {rs.word_matrix(parse(w)): w for w in TABLE8}
for w in sorted(noncomm, key=rs.length):
    print(f"  len {rs.length(w):2d} word {word_str(rs.reduced_word(w))} L={sorted(rs.left_descents(w))} R={sorted(rs.right_descents(w))} in table: {tab.get(w)}")
assert set(noncomm) == set(tab)
print("E8 noncommuting terminals == Table 1 rows: True")
print(f"total time {time.time()-t0:.1f}s")
