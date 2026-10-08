"""Lower bound for l(b_{r,0}) = |N(r_beta)| by counting positive roots alpha in the finite
parabolics D_{4r} (nodes 1..4r) and A_{4r} (nodes 0..4r-1) with (alpha, beta_r) > 0.
Any such alpha has a missing coordinate j with (r_beta alpha)_j = -(alpha,beta) beta_j < 0,
so r_beta alpha < 0 and alpha in N(r_beta).  Also check the translation-length criterion:
(eta, beta_r) <= -1 for all positive E8 roots eta (on J = {0..7,4r}) with (eta,gamma) >= 1,
and (eta, beta_r) <= a - 1 for negative such eta (i.e. (eta',beta_r) >= 1-a for positive
eta' with (eta',gamma) <= -1)."""
import sys
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-uniform-all-k')
from uni import En, beta_rk, uniform_data


def positive_roots(G, nodes):
    """Positive roots of the finite parabolic on `nodes` by closure from simple roots."""
    n = G.n
    roots = set()
    frontier = []
    for i in nodes:
        v = tuple(1 if j == i else 0 for j in range(n))
        roots.add(v)
        frontier.append(v)
    while frontier:
        nxt = []
        for v in frontier:
            for i in nodes:
                c = G.pair(v, i)
                if c < 0:
                    u = list(v); u[i] -= c; u = tuple(u)
                    if u not in roots:
                        roots.add(u); nxt.append(u)
        frontier = nxt
    return roots


for r in range(3, 8):
    n, beta, gamma, delta, a = uniform_data(r)
    G = En(n)
    D = positive_roots(G, list(range(1, n)))
    Achain = positive_roots(G, list(range(0, n - 1)))
    assert len(D) == (n - 1) * (n - 2) and len(Achain) == (n - 1) * n // 2
    inD = {v for v in D if G.form(v, beta) > 0}
    inA = {v for v in Achain if G.form(v, beta) > 0}
    union = inD | inA
    # sanity: every such root is an inversion
    b = G.reflection(beta)
    for v in union:
        img = tuple(sum(b[i][k] * v[i] for i in range(n)) for k in range(n))
        assert G.neg(img)
    print(f"r={r} E{n}: |N cap D_{n-1}^+|={len(inD)} |N cap A_{n-1}^+|={len(inA)} |union|={len(union)} "
          f"target 8r^2+3={8*r*r+3}")
    # translation criterion on the E8 parabolic J
    J = list(range(0, 8)) + [n - 1]
    E8 = positive_roots(G, J)
    assert len(E8) == 120
    bad = []
    for eta in E8:
        eg = G.form(eta, gamma)
        eb = G.form(eta, beta)
        if eg >= 1 and not eb <= -1:
            bad.append(('pos', eta, eg, eb))
        if eg <= -1 and not eb >= 1 - a:   # negative root -eta: (-eta,beta) <= a-1
            bad.append(('neg', eta, eg, eb))
    npos = sum(1 for eta in E8 if G.form(eta, gamma) >= 1)
    tot = sum(G.form(eta, gamma) for eta in E8 if G.form(eta, gamma) >= 1)
    print(f"   translation criterion violations: {len(bad)}; #pos E8 roots with (eta,gamma)>=1: {npos}, "
          f"sum of pairings (=l(T)): {tot}, delta pairing (delta,beta)={G.form(delta, beta)}, (gamma,beta)={G.form(gamma, beta)}")
    for x in bad[:5]:
        print('   ', x)
