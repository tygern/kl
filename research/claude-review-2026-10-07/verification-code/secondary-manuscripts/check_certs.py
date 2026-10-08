import json, sys
sys.path.insert(0, '.')
from klfinite import *
d = json.load(open('/Users/ic/workspace/kl/results/computation.json'))
groups = {}
def G(t, n):
    if (t, n) not in groups: groups[(t, n)] = Group(t, n)
    return groups[(t, n)]
print("group orders:", {k: len(G(*k).elems) for k in [('A',4),('B',3),('D',4),('B',4)]})
for i, c in enumerate(d['certificates']):
    st, tt = c['source']['type'], c['target']['type']
    sn = int(st[1]) + (1 if st[0]=='A' else 0); tn = int(tt[1]) + (1 if tt[0]=='A' else 0)
    Gs, Gt = G(st[0], sn), G(tt[0], tn)
    xs, ws = tuple(c['source']['bottom']['one_line']), tuple(c['source']['top']['one_line'])
    xt, wt = tuple(c['target']['bottom']['one_line']), tuple(c['target']['top']['one_line'])
    pairs = [(tuple(m['source']['one_line']), tuple(m['target']['one_line'])) for m in c['map']]
    ok = check_bijection(Gs, xs, ws, Gt, xt, wt, pairs)
    Ps, Pt = Gs.P(xs, ws), Gt.P(xt, wt)
    rv = [0]*(Gs.len[ws]-Gs.len[xs]+1)
    for z in Gs.interval(xs, ws): rv[Gs.len[z]-Gs.len[xs]] += 1
    fc_s, fc_t = Gs.is_fc(xs), Gt.is_fc(xt)
    # every subinterval polynomial equal
    f = dict(pairs); I = Gs.interval(xs, ws)
    sub_ok = all(Gs.P(a, b) == Gt.P(f[a], f[b]) for a in I for b in I if Gs.leq(a, b))
    rec_ok = Gs.check_reciprocity(xs, ws) and Gt.check_reciprocity(xt, wt)
    print(f"cert {i}: {st} [{xs},{ws}] len ({Gs.len[xs]},{Gs.len[ws]}) -> {tt} [{xt},{wt}] len ({Gt.len[xt]},{Gt.len[wt]})")
    print(f"   bijection order-iso: {ok}; P_src={Ps} P_tgt={Pt} json={tuple(c['polynomial_coefficients'])}; mu_src={Gs.mu(xs,ws)} json mu={c['mu']}")
    print(f"   rank vector {rv} json {c['rank_vector']}; FC bottom src={fc_s} tgt={fc_t}; all subinterval P equal: {sub_ok}; reciprocity: {rec_ok}")
    print(f"   src top reduced word {Gs.reduced_word(ws)}, tgt bottom word {Gt.reduced_word(xt)}, tgt top word {Gt.reduced_word(wt)}")
# B3 cover reflection order-4 claim
Gb = G('B', 3)
xt, wt = (3,2,-1), (2,-1,-3)
covers = [z for z in Gb.interval(xt, wt) if Gb.len[z] == Gb.len[xt]+1]
def compose(u, v):  # (u*v)(i) = u(v(i)) as signed perms acting on positions? We use right-multiplication convention: w*s acts on positions. Define product of one-line permutations as functions.
    n = len(u)
    return tuple((1 if v[i] > 0 else -1) * u[abs(v[i])-1] for i in range(n))
def order(t):
    e = tuple(range(1, len(t)+1)); k = 1; w = t
    while w != e: w = compose(w, t); k += 1
    return k
print("B3 right cover labels x^{-1} z:", [(z, compose(Gb.inv(xt), z)) for z in covers])
labs = [compose(Gb.inv(xt), z) for z in covers]
for a, b in combinations(labs, 2):
    print("  product order of", a, b, "=", order(compose(a, b)))
# D4->B4 rank-5 example details
Gd, Gb4 = G('D', 4), G('B', 4)
print("D4 len (-4,-3,-2,-1):", Gd.len[(-4,-3,-2,-1)], "word", Gd.reduced_word((-4,-3,-2,-1)), "support all 4:", set(Gd.reduced_word((-4,-3,-2,-1)))=={0,1,2,3})
print("B4 words:", Gb4.reduced_word((2,4,3,-1)), Gb4.reduced_word((2,3,-1,-4)), "lens", Gb4.len[(2,4,3,-1)], Gb4.len[(2,3,-1,-4)])
# B4 right cover labels order check
xt, wt = (2,4,3,-1), (2,3,-1,-4)
covers = [z for z in Gb4.interval(xt, wt) if Gb4.len[z] == Gb4.len[xt]+1]
labs = [compose(Gb4.inv(xt), z) for z in covers]
print("B4 rank5 cover-label product orders:", sorted(order(compose(a,b)) for a,b in combinations(labs,2)))
# B2 butterfly: search all intervals in A5, D5, B3 for rank vector (1,2,2,2,1)
for t, n in [('A',5),('D',5),('B',3),('D',4)]:
    g = G(t, n); found = 0
    for w in g.elems:
        for x in g.ideal(w):
            if g.len[w]-g.len[x] == 4:
                I = g.interval(x, w)
                if len(I) == 8:
                    found += 1
    print(f"{t}{n}: rank-4 intervals with 8 elements: {found}")
# Coset lemma sanity in D4 with J={1,2,3} (type A3) and B3 with J={1,2}
for (t,n,J) in [('D',4,[1,2,3]),('B',3,[1,2]),('D',4,[0,2,3])]:
    g = G(t, n)
    WJ = {w for w in g.elems if set(g.reduced_word(w)) <= set(J)}
    minreps = [a for a in g.elems if all(g.len[act(t,a,s)] > g.len[a] for s in J)]
    bad = 0; tested = 0
    for a in minreps:
        for u in WJ:
            for v in WJ:
                au = compose(a, u); av = compose(a, v)  # right multiplication by u: a*u ; our compose(u,v)(i)=u(v(i)) corresponds to a*u in word sense? check below
                if g.leq(u, v) != g.leq(au, av): bad += 1
                if g.leq(u, v) and g.P(u, v) != g.P(au, av): bad += 1
                tested += 1
    print(f"coset check {t}{n} J={J}: tested {tested}, violations {bad}")
