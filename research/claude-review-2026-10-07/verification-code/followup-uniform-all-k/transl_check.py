"""Translation criterion for l(b_{r,k}) = l(b_{r,0}) + 116k.
Finite E8 parabolic J0 = {0,...,6, 4r} (affine node 7 omitted); affine E8 = J0 + {7}.
N(T^k) = { eta + m delta : eta in Phi(E8), (eta,gamma) >= 1, eps(eta) <= m < k(eta,gamma)+eps(eta) }.
Length additivity of T^k * b * T^{-k} holds iff (alpha, beta_r) <= -1 for all alpha in N(T^k).
Since (eta + m delta, beta_r) = (eta,beta_r) - a m with a = 2r-4 >= 2 and m >= eps(eta):
  need (eta,beta_r) <= -1 for positive eta with (eta,gamma) >= 1, and
       (eta,beta_r) <= a-1 for negative eta with (eta,gamma) >= 1.
(eta,beta_r) = (r-2)(eta_0 - eta_1) + f(eta) is affine-linear in r, so the check for all r >= 3
reduces to checking the coefficient and the value at r = 3."""
import sys
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-uniform-all-k')
from uni import En, uniform_data

def positive_roots(G, nodes):
    n = G.n
    roots = set(); frontier = []
    for i in nodes:
        v = tuple(1 if j == i else 0 for j in range(n)); roots.add(v); frontier.append(v)
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
        if len(roots) > 100000: raise RuntimeError("infinite")
    return roots

for r in range(3, 9):
    n, beta, gamma, delta, a = uniform_data(r)
    G = En(n)
    J0 = list(range(0, 7)) + [n - 1]
    E8 = positive_roots(G, J0)
    assert len(E8) == 120
    assert G.form(delta, delta) == 0 and G.form(gamma, delta) == 0 and G.form(gamma, gamma) == 2
    assert all(G.pair(delta, j) == 0 for j in J0 + [7])          # delta is the affine null root
    assert G.form(delta, beta) == -a and G.form(gamma, beta) == -r
    lT = sum(G.form(eta, gamma) for eta in E8 if G.form(eta, gamma) >= 1)
    viol = []
    for eta in E8:
        eg = G.form(eta, gamma); eb = G.form(eta, beta)
        if eg >= 1 and eb > -1: viol.append(('pos', eta, eg, eb))
        if eg <= -1 and -eb > a - 1: viol.append(('neg', eta, eg, eb))
    # r-coefficient check: (eta,beta_r) = (r-2)(eta_0-eta_1) + f(eta)
    coef_bad = [eta for eta in E8 if G.form(eta, gamma) >= 1 and eta[0] - eta[1] > 0]
    coef_bad2 = [eta for eta in E8 if G.form(eta, gamma) <= -1 and eta[0] - eta[1] < -2]  # -(eta,beta) grows faster than a-1 = 2r-5
    print(f"r={r}: #pos roots with (eta,gamma)>=1: {sum(1 for e in E8 if G.form(e,gamma)>=1)}, l(T)={lT}, violations={len(viol)}, "
          f"r-coefficient problems: {len(coef_bad)}, {len(coef_bad2)}")
    for v in viol[:5]: print('   ', v)
