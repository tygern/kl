"""Independent verification (exact integer arithmetic, no shared code) of finding
followup-uniform-all-k-D2: every Bruhat cover of b_{r,k} fails to be fully commutative,
and the supporting lemmas (sbs is a non-simple reflection; non-simple reflections in a
simply laced group are not FC; L(bs) = I_r).

Conventions follow results/exceptional-leading.tex, Section 'The rank-uniform
construction': E_{4r+1} numbered as a chain 0-1-...-(4r-1) with node 4r attached to 2.
Coordinates are in the simple-root basis; (alpha_i,alpha_i)=2, adjacent pairing -1.
"""
import sys, itertools

def En_edges(n, attach=2):
    # chain 0..n-2, node n-1 attached to `attach`
    E = set()
    for j in range(n - 2):
        E.add((j, j + 1)); E.add((j + 1, j))
    E.add((attach, n - 1)); E.add((n - 1, attach))
    return E

class Coxeter:
    def __init__(self, n, edges):
        self.n = n
        self.edges = edges
        self.A = [[2 if i == j else (-1 if (i, j) in edges else 0) for j in range(n)] for i in range(n)]
    def pair(self, u, v):
        n = self.n
        return sum(u[i] * self.A[i][j] * v[j] for i in range(n) for j in range(n))
    def simple(self, j):
        return tuple(1 if i == j else 0 for i in range(self.n))
    def reflect(self, beta, v):
        c = self.pair(v, beta)
        return tuple(v[i] - c * beta[i] for i in range(self.n))
    def refl_matrix(self, beta):
        # columns are images of simple roots
        return tuple(self.reflect(beta, self.simple(j)) for j in range(self.n))
    def apply(self, M, v):
        n = self.n
        return tuple(sum(M[j][i] * v[j] for j in range(n)) for i in range(n))
    def mul(self, M, N):
        # (MN)(alpha_j) = M(N alpha_j)
        return tuple(self.apply(M, N[j]) for j in range(self.n))
    def gen(self, s):
        return self.refl_matrix(self.simple(s))
    def identity(self):
        return tuple(self.simple(j) for j in range(self.n))
    def is_neg(self, v):
        return all(c <= 0 for c in v) and any(c < 0 for c in v)
    def is_pos(self, v):
        return all(c >= 0 for c in v) and any(c > 0 for c in v)
    def R(self, M):
        return [s for s in range(self.n) if self.is_neg(M[s])]
    def inverse(self, M):
        # Coxeter matrices in the root basis: inverse = transpose w.r.t. form; compute via word
        w = self.word(M)
        N = self.identity()
        for s in reversed(w):
            N = self.mul(N, self.gen(s))
        return N
    def L(self, M):
        return self.R(self.inverse(M))
    def word(self, M):
        """Reduced word by stripping right descents: returns list w with M = s_{w[0]}...s_{w[-1]}."""
        w = []
        cur = M
        while True:
            R = self.R(cur)
            if not R:
                break
            s = R[0]
            cur = self.mul(cur, self.gen(s))
            w.append(s)
        assert cur == self.identity()
        return list(reversed(w))
    def length(self, M):
        return len(self.word(M))
    def from_word(self, w):
        M = self.identity()
        for s in w:
            M = self.mul(M, self.gen(s))
        return M
    def word_is_fc(self, w):
        """Stembridge: a reduced word is in a braid-free commutation class iff between any two
        consecutive occurrences of s there are at least two letters not commuting with s."""
        n = len(w)
        last = {}
        for k, s in enumerate(w):
            if s in last:
                i = last[s]
                cnt = sum(1 for j in range(i + 1, k) if (w[j], s) in self.edges)
                assert cnt >= 1, "word not reduced"
                if cnt == 1:
                    return False
            last[s] = k
        return True
    def is_fc(self, M):
        return self.word_is_fc(self.word(M))

def beta_data(r):
    n = 4 * r + 1
    beta = [0] * n
    beta[0], beta[1], beta[2], beta[4 * r] = r - 1, r, 2 * r - 1, r
    for j in range(3, 4 * r):
        beta[j] = 2 * r - j // 2
    delta = [0] * n; gamma = [0] * n
    for j, c in enumerate((2, 4, 6, 5, 4, 3, 2, 1)):
        delta[j] = c
    delta[4 * r] = 3
    for j, c in enumerate((1, 2, 3, 2, 2, 1, 1, 0)):
        gamma[j] = c
    gamma[4 * r] = 1
    a = 2 * r - 4
    I = [0] + list(range(3, 4 * r, 2)) + [4 * r]
    return tuple(beta), tuple(gamma), tuple(delta), a, I

def beta_rk(r, k):
    beta, gamma, delta, a, I = beta_data(r)
    n = 4 * r + 1
    return tuple(beta[i] - a * k * gamma[i] + (a * k * k + r * k) * delta[i] for i in range(n))

def check_family(r, k, enumerate_covers=True):
    n = 4 * r + 1
    G = Coxeter(n, En_edges(n))
    beta, gamma, delta, a, I = beta_data(r)
    b_root = beta_rk(r, k)
    assert G.pair(b_root, b_root) == 2, "beta_{r,k} not norm 2"
    assert all(c >= 1 for c in b_root)
    B = G.refl_matrix(b_root)
    wB = G.word(B)
    lb = len(wB)
    assert lb % 2 == 1
    assert G.R(B) == I and G.L(B) == I, (G.R(B), G.L(B), I)
    assert len(set(wB)) == n, "not full support"
    print(f"r={r} k={k} E{n}: l(b)={lb} (8r^2+3+116k={8*r*r+3+116*k}), L=R=I_r={I}, |I|={len(I)}, l(i(I))={len(I)} odd, l(b)-1={lb-1} even")
    # Step (2): for s in I, bs: s in L(bs); sbs = r_{s beta}; s beta positive non-simple; bs not FC
    for s in I:
        ms = G.pair(G.simple(s), b_root)
        assert ms >= 1
        sbeta = G.reflect(G.simple(s), b_root)
        assert G.is_pos(sbeta) and sum(sbeta) >= 2, "s beta not positive non-simple"
        BS = G.mul(B, G.gen(s))
        assert G.length(BS) == lb - 1
        Lbs = G.L(BS); Rbs = G.R(BS)
        assert s in Lbs
        SBS = G.mul(G.gen(s), BS)
        assert SBS == G.refl_matrix(sbeta), "sbs != r_{s beta}"
        assert G.length(SBS) == lb - 2
        assert not G.is_fc(SBS), "sbs is FC?!"
        assert not G.is_fc(BS), "bs is FC?!"
        assert Lbs == I and Rbs == sorted(set(I) - {s}), (Lbs, Rbs)
    print(f"   all {len(I)} elements bs: s in L(bs), sbs = r_(s beta) non-simple reflection, bs and sbs non-FC, L(bs)=I_r, R(bs)=I_r-{{s}}: OK")
    if enumerate_covers:
        covers = {}
        for i in range(lb):
            w = wB[:i] + wB[i + 1:]
            M = G.from_word(w)
            if G.length(M) == lb - 1:
                covers[M] = w
        fc = [M for M in covers if G.is_fc(M)]
        # classify: which covers are of form sb or bs
        sb_bs = set()
        for s in I:
            sb_bs.add(G.mul(G.gen(s), B)); sb_bs.add(G.mul(B, G.gen(s)))
        others = [M for M in covers if M not in sb_bs]
        bad = [M for M in others if not (set(I) <= set(G.L(M)) and set(I) <= set(G.R(M)))]
        print(f"   distinct Bruhat covers={len(covers)}  FC covers={len(fc)}  covers of form sb/bs={len(sb_bs & set(covers))}/{len(sb_bs)}  other covers={len(others)}, of which with I not contained in L and R: {len(bad)}")
        assert not fc and not bad

def check_reflection_lemma_finite():
    """All reflections in E6, E7, E8 (nodes: chain 0..n-2, node n-1 attached to 2):
    the non-simple ones are not FC, and l(r_beta) = 2N+1 for every lowering chain."""
    for n in (6, 7, 8):
        G = Coxeter(n, En_edges(n))
        # generate positive roots by closure under simple reflections
        roots = set(G.simple(j) for j in range(n))
        frontier = list(roots)
        while frontier:
            new = []
            for v in frontier:
                for s in range(n):
                    u = G.reflect(G.simple(s), v)
                    if G.is_pos(u) and u not in roots:
                        roots.add(u); new.append(u)
            frontier = new
        nonfc = 0; nonsimple = 0
        memo = {}
        def Nset(v):
            if v in memo:
                return memo[v]
            if sum(v) == 1:
                memo[v] = frozenset({0}); return memo[v]
            out = set()
            for s in range(n):
                if G.pair(G.simple(s), v) > 0:
                    out |= {d + 1 for d in Nset(G.reflect(G.simple(s), v))}
            memo[v] = frozenset(out)
            return memo[v]
        for beta in roots:
            M = G.refl_matrix(beta)
            l = G.length(M)
            if sum(beta) == 1:
                assert l == 1 and G.is_fc(M)
                continue
            nonsimple += 1
            assert not G.is_fc(M), ("FC non-simple reflection", beta)
            nonfc += 1
            # all lowering chains have the same length N with l = 2N+1
            Ns = Nset(beta)
            assert len(Ns) == 1 and l == 2 * next(iter(Ns)) + 1, (beta, l, Ns)
        print(f"E{n}: {len(roots)} positive roots, {nonsimple} non-simple reflections, all non-FC: {nonfc}; l(r_beta)=2N+1 for every lowering chain: OK")

def check_reflection_lemma_affine(maxheight=14):
    """E9 (affine): all positive real roots of height <= maxheight: non-simple reflections not FC."""
    n = 9
    G = Coxeter(n, En_edges(n))
    roots = set(G.simple(j) for j in range(n))
    frontier = list(roots)
    while frontier:
        new = []
        for v in frontier:
            for s in range(n):
                u = G.reflect(G.simple(s), v)
                if G.is_pos(u) and sum(u) <= maxheight and u not in roots:
                    roots.add(u); new.append(u)
        frontier = new
    cnt = 0
    for beta in roots:
        if sum(beta) == 1:
            continue
        M = G.refl_matrix(beta)
        assert not G.is_fc(M), beta
        cnt += 1
    print(f"E9 affine: {cnt} non-simple real reflections of height <= {maxheight}, all non-FC: OK")

def check_e10():
    """Appendix E10 family: same cover argument."""
    n = 10
    G = Coxeter(n, En_edges(n))
    delta = (2, 4, 6, 5, 4, 3, 2, 1, 0, 3)
    gamma = (1, 2, 3, 2, 2, 1, 1, 0, 0, 1)
    beta0 = (3, 7, 10, 9, 7, 6, 4, 3, 1, 6)
    I = [1, 3, 5, 7, 9]
    for k in range(0, 3):
        beta = tuple(beta0[i] - k * gamma[i] + (k * k + 3 * k) * delta[i] for i in range(n))
        assert G.pair(beta, beta) == 2
        B = G.refl_matrix(beta)
        lb = G.length(B)
        assert G.L(B) == I and G.R(B) == I
        wB = G.word(B)
        covers = {}
        for i in range(lb):
            w = wB[:i] + wB[i + 1:]
            M = G.from_word(w)
            if G.length(M) == lb - 1:
                covers[M] = w
        fc = [M for M in covers if G.is_fc(M)]
        print(f"E10 k={k}: l(b_k)={lb}, distinct covers={len(covers)}, FC covers={len(fc)}")
        assert not fc

if __name__ == "__main__":
    sys.setrecursionlimit(10000)
    check_reflection_lemma_finite()
    check_reflection_lemma_affine()
    for (r, k) in [(3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (5, 0)]:
        check_family(r, k, enumerate_covers=True)
    for (r, k) in [(6, 0), (7, 0), (8, 0), (5, 3)]:
        check_family(r, k, enumerate_covers=False)
    check_e10()
    print("ALL CHECKS PASSED")
