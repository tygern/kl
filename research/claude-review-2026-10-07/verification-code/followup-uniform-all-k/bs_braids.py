"""Locate the braid obstructions in b s (s in I_r): for a reduced word of x = b s, list the
generators t whose consecutive occurrences have exactly one neighbour occurrence between them
(convex chain t u t in the heap), with the neighbour u.  Also record multiplicities."""
import sys
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-uniform-all-k')
from uni import En, beta_rk, I_r
from collections import Counter

def braids(G, word):
    last = {}
    out = []
    for p, s in enumerate(word):
        if s in last:
            between = [t for t in word[last[s] + 1:p] if t in G.adjset[s]]
            if len(between) == 1:
                out.append((s, between[0]))
        last[s] = p
    return out

for r, k in [(3, 0), (3, 1), (4, 0), (4, 1), (5, 0)]:
    n = 4 * r + 1
    G = En(n)
    b = G.reflection(beta_rk(r, k))
    print(f"== r={r} k={k} E{n}")
    for s in I_r(r):
        x = G.rmul(b, s)
        w = G.word(x)
        br = braids(G, w)
        mult = Counter(w)
        print(f" s={s}: braids (t,u) in stripped word: {sorted(set(br))}; mult of 0,1,2,p,q: {mult[0]},{mult[1]},{mult[2]},{mult[n-1]},{mult[n-2]}")
    w = G.word(b)
    print(f" b itself: braids {sorted(set(braids(G, w)))}; word={w}")
