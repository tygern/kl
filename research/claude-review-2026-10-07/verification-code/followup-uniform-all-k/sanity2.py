import sys
sys.path.insert(0, '.')
from uni import En, beta_rk, I_r
for r in (3, 4, 5, 6):
    for k in (0, 1, 2):
        n = 4*r+1; G = En(n); beta = beta_rk(r, k); b = G.reflection(beta)
        c1=c2=c3=c4=c5=True; zeros=[]
        for s in I_r(r):
            x = G.rmul(b, s); sbs = G.lmul(s, x); sbeta = G.refl_vec(beta, s)
            c1 &= (sbs == G.reflection(sbeta)); c2 &= (sum(sbeta) > 1); c3 &= (s in G.L(x))
            c4 &= (G.length(sbs) == G.length(x) - 1); c5 &= all(c >= 0 for c in sbeta)
            zeros += [(s, j) for j, c in enumerate(sbeta) if c == 0]
        print(f"r={r} k={k}: sbs=r_(s beta):{c1} nonsimple:{c2} s in L(bs):{c3} l(sbs)=l(bs)-1:{c4} s beta>=0:{c5} zero coords:{zeros}")
