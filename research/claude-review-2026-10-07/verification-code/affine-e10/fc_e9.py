"""Enumerate the fully commutative elements of affine E8 (= E9) by weak-order
closure using the recurrence: w FC  <=>  R(w) pairwise commuting and ws FC
for all s in R(w).  Elements are stored as tuple-of-columns matrices.
Also check which Bruhat covers of b0 are FC and the FC maximum length."""
import sys, pickle, time
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/affine-e10')
from cox import *

n = 9
edges = [(i, i + 1) for i in range(7)] + [(2, 8)]
A = cartan(n, edges)
adj = [[t for t in range(n) if A[s][t] == -1] for s in range(n)]


def right_mult(cols, s):
    """(w s)(alpha_j) = w(alpha_j) - A[s][j] w(alpha_s)."""
    cs = cols[s]
    new = list(cols)
    new[s] = tuple(-x for x in cs)
    for t in adj[s]:
        new[t] = tuple(a + b for a, b in zip(cols[t], cs))
    return tuple(new)


def rdesc(cols):
    return [j for j in range(n) if cols[j][0] <= 0 and all(x <= 0 for x in cols[j]) and any(x < 0 for x in cols[j])]


def commuting(D):
    return all(A[s][t] == 0 for s in D for t in D if s < t)


t0 = time.time()
e = identity(n)
FC = {e: 0}
level = [e]
length = 0
maxlen = 0
while level:
    nxt = []
    seen = set()
    for w in level:
        D = set(rdesc(w))
        for s in range(n):
            if s in D:
                continue
            ws = right_mult(w, s)
            if ws in seen:
                continue
            seen.add(ws)
            D2 = rdesc(ws)
            if not commuting(D2):
                continue
            ok = True
            for t in D2:
                if right_mult(ws, t) not in FC:
                    ok = False
                    break
            if ok:
                FC[ws] = length + 1
                nxt.append(ws)
    level = nxt
    length += 1
    if nxt:
        maxlen = length
    print(f"length {length}: {len(nxt)} new FC elements (total {len(FC)})", flush=True)
print("TOTAL FC elements in affine E8:", len(FC), " max length:", maxlen, f" ({time.time()-t0:.1f}s)")

# length distribution check
from collections import Counter
dist = Counter(FC.values())
print("length distribution:", [dist[i] for i in range(maxlen + 1)])

with open('/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/affine-e10/covers.pkl', 'rb') as f:
    covers = pickle.load(f)
print("covers of b0:", len(covers), " FC among them:", sum(1 for c in covers if c in FC))
# Additional direct check of each cover: compute reduced word and check for braid via FC membership consistency
for c in covers:
    L, w = length_by_stripping(A, c)
    assert L == 26
print("all covers have length 26 (re-verified)")

# Sanity: finite E8 parabolic FC count (omit node 7) should be 10846
fin = sum(1 for w in FC if all(w[j][7] == 0 for j in range(n)) and all(w[7][i] == (1 if i == 7 else 0) for i in range(n)))
print("FC elements supported in finite E8 parabolic (omit node 7):", fin)
