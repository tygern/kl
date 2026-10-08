"""Verify l(b_{r,k}) = 8r^2 + 3 + 116k, descents, terminality, depth relation for many (r,k)."""
import sys, time
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-uniform-all-k')
from uni import En, beta_rk, I_r, uniform_data

pairs = [(r, k) for r in range(3, 8) for k in range(0, 5)] + [(8, 0), (8, 1), (9, 0), (10, 0)]
for r, k in pairs:
    t0 = time.time()
    n = 4 * r + 1
    G = En(n)
    beta = beta_rk(r, k)
    assert G.form(beta, beta) == 2, "norm"
    b = G.reflection(beta)
    word = G.word(b)
    l = len(word)
    R = G.R(b)
    Lb = G.L(b)
    term = G.terminal(b)
    fc = G.is_fc_word(word)
    full = len(set(word)) == n
    pred = 8 * r * r + 3 + 116 * k
    dp = G.depth(beta) if l < 400 else None
    print(f"r={r} k={k} E{n}: len={l} formula={pred} ok={l == pred} R=I_r:{R == I_r(r)} L=R:{Lb == R} "
          f"terminal={term} fullsupp={full} FC={fc} depth={dp} 2dp-1={None if dp is None else 2*dp-1} "
          f"[{time.time()-t0:.1f}s]")
    sys.stdout.flush()
