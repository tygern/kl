"""Brute-force the specific Lemma maxI instances used in Theorem affine (E9) and
Proposition thm:indefinite (E10): over the complete FC catalogue, the only FC x with
I subset of L(x) and R(x) must be i(I).  Also sanity-check Stembridge's FC criterion
against the recurrence-generated catalogue, and the FC count / max length."""
import sys, time
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/completeness')
from en import En

def run(n, I, b0word=None):
    G = En(n)
    t0 = time.time()
    levels, fc = G.fc_catalogue()
    print(f"E{n}: FC count {len(fc)}, max length {len(levels)-1}, {time.time()-t0:.1f}s")
    iI = G.iI(I)
    hits = []
    for x in fc:
        Rx = set(G.R(x))
        if not set(I) <= Rx:
            continue
        Lx = set(G.L(x))
        if set(I) <= Lx:
            hits.append(x)
    print(f"  FC x with I={I} in L(x) and R(x): {len(hits)}; equals i(I): {hits == [iI]}")
    # Stembridge criterion consistency: every catalogue word passes; sample non-FC fails
    bad = 0
    for lvl in levels:
        for x in lvl:
            if not G.is_fc_word(G.word(x)):
                bad += 1
    print(f"  Stembridge-criterion failures on catalogue: {bad}")
    if b0word is not None:
        b0 = G.from_word(b0word)
        print(f"  b0 length {G.length(b0)}, R={G.R(b0)}, L={G.L(b0)}, terminal={G.terminal(b0)}, FC={G.is_fc_word(G.word(b0))}")
        # eligible FC x below b0 : count via descent inclusion only (Bruhat comparison skipped)
    return G, levels, fc

if __name__ == '__main__':
    run(9, [1, 3, 5, 7, 8], [int(c) for c in '312875645234123012856745231'])
    run(10, [1, 3, 5, 7, 9])
