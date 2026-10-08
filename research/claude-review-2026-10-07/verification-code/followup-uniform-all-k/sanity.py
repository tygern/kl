import sys
sys.path.insert(0, '.')
from uni import En, uniform_data, beta_rk, I_r
# (1) l(T^k) = 58k in E_{4r+1} for r=3,4 and k=1..3
for r in (3, 4):
    n, beta, gamma, delta, a = uniform_data(r)
    G = En(n)
    T = G.mul(G.reflection(gamma), G.reflection(tuple(g + d for g, d in zip(gamma, delta))))
    Tk = G.e
    for k in range(1, 4):
        Tk = G.mul(T, Tk)
        print(f"r={r} k={k}: l(T^k)={G.length(Tk)} (58k={58*k})")
# (2) every reflection of length >= 3 in E6, E7, E8 (all positive roots) is non-FC; l(r_beta)=2d+1
def positive_roots(G):
    n = G.n; roots=set(); fr=[]
    for i in range(n):
        v=tuple(1 if j==i else 0 for j in range(n)); roots.add(v); fr.append(v)
    while fr:
        nx=[]
        for v in fr:
            for i in range(n):
                c=G.pair(v,i)
                if c<0:
                    u=list(v); u[i]-=c; u=tuple(u)
                    if u not in roots: roots.add(u); nx.append(u)
        fr=nx
    return roots
for n in (6, 7, 8):
    G = En(n); R = positive_roots(G); bad = 0; depthbad = 0
    for v in R:
        t = G.reflection(v); w = G.word(t)
        if len(w) >= 3 and G.is_fc_word(w): bad += 1
        if len(w) != 2 * G.depth(v) + 1: depthbad += 1
    print(f"E{n}: {len(R)} positive roots; FC reflections of length>=3: {bad}; l != 2d+1: {depthbad}")
# (3) s b s = r_{s beta} is a non-simple reflection; s in L(bs) for all s in I_r (r=3..6,k=0..2)
for r in (3, 4, 5, 6):
    for k in (0, 1, 2):
        n = 4*r+1; G = En(n); beta = beta_rk(r, k); b = G.reflection(beta); ok = True
        for s in I_r(r):
            x = G.rmul(b, s); sbs = G.lmul(s, x)
            sbeta = G.refl_vec(beta, s)
            ok &= (sbs == G.reflection(sbeta)) and (sum(sbeta) > 1) and (s in G.L(x)) and (G.length(sbs) == G.length(x) - 1) and all(c >= 1 for c in sbeta)
        print(f"r={r} k={k}: s b s = r_(s beta) non-simple positive, s in L(bs), lengths consistent: {ok}")
