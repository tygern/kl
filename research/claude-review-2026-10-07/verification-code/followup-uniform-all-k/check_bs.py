"""For many (r,k): verify (i) L(b s) = I_r and R(b s) = I_r \\ {s} for every s in I_r,
(ii) b s is not FC (heap test on a reduced word), (iii) by inversion the same for s b.
These are exactly the elements that the lifting-property reduction leaves as the only
possible FC Bruhat covers of b = b_{r,k}."""
import sys, time
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-uniform-all-k')
from uni import En, beta_rk, I_r

pairs = [(r, k) for r in range(3, 9) for k in range(0, 6)] + [(9, 0), (9, 3), (10, 0), (10, 2), (12, 0), (12, 1)]
for r, k in pairs:
    t0 = time.time()
    n = 4 * r + 1
    G = En(n)
    b = G.reflection(beta_rk(r, k))
    I = I_r(r)
    allok = True
    nonfc = 0
    for s in I:
        x = G.rmul(b, s)
        word = G.word(x)
        assert len(word) == 8 * r * r + 2 + 116 * k
        L = G.L(x); R = G.R(x)
        ok = (L == I) and (R == sorted(set(I) - {s}))
        fc = G.is_fc_word(word)
        if not fc:
            nonfc += 1
        allok = allok and ok and not fc
    print(f"r={r} k={k} E{n}: all {len(I)} elements b s have L=I_r, R=I_r-{{s}}: {allok}; non-FC: {nonfc}/{len(I)} [{time.time()-t0:.1f}s]")
    sys.stdout.flush()
