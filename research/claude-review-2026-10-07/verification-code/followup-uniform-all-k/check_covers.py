"""Enumerate all Bruhat covers of b_{r,k} (deleting one letter from a reduced word and
keeping length l-1; by the subword property this is exhaustive) and test each for full
commutativity with the heap criterion.  Also report the FC-ness via the recurrence-free
heap test on the stripped reduced word."""
import sys, time
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-uniform-all-k')
from uni import En, beta_rk

pairs = [(3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (5, 0), (6, 0), (7, 0)]
if len(sys.argv) > 1:
    pairs = [tuple(map(int, p.split(','))) for p in sys.argv[1:]]
for r, k in pairs:
    t0 = time.time()
    n = 4 * r + 1
    G = En(n)
    b = G.reflection(beta_rk(r, k))
    cov = G.covers(b)
    fc = []
    minsupp = n
    for v, sub in cov.items():
        w = G.word(v)
        assert len(w) == len(sub)
        minsupp = min(minsupp, len(set(w)))
        if G.is_fc_word(w):
            fc.append(w)
    print(f"r={r} k={k} E{n}: l(b)={8*r*r+3+116*k} distinct covers={len(cov)} FC covers={len(fc)} "
          f"min cover support={minsupp} [{time.time()-t0:.1f}s]")
    sys.stdout.flush()
