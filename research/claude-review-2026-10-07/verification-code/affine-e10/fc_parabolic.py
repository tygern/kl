import sys, time
sys.path.insert(0, '.')
from cox import *
# finite E8 with the E9 numbering minus node 7: nodes 0..6, 8; node 8 attached to 2
n = 9
edges = [(i, i + 1) for i in range(6)] + [(2, 8)]   # no edge 6-7; node 7 isolated (acts as A1 factor)
A = cartan(n, edges)
adj = [[t for t in range(n) if A[s][t] == -1] for s in range(n)]
def right_mult(cols, s):
    cs = cols[s]; new = list(cols); new[s] = tuple(-x for x in cs)
    for t in adj[s]: new[t] = tuple(a + b for a, b in zip(cols[t], cs))
    return tuple(new)
def rdesc(cols):
    return [j for j in range(n) if all(x <= 0 for x in cols[j]) and any(x < 0 for x in cols[j])]
def commuting(D): return all(A[s][t] == 0 for s in D for t in D if s < t)
e = identity(n); FC = {e}; level = [e]; length = 0; mx = 0
active = [s for s in range(n) if s != 7]
while level:
    nxt = []; seen = set()
    for w in level:
        D = set(rdesc(w))
        for s in active:
            if s in D: continue
            ws = right_mult(w, s)
            if ws in seen: continue
            seen.add(ws)
            D2 = rdesc(ws)
            if not commuting(D2): continue
            if all(right_mult(ws, t) in FC for t in D2):
                FC.add(ws); nxt.append(ws)
    level = nxt; length += 1
    if nxt: mx = length
print("finite E8 FC count:", len(FC), "max length", mx)
