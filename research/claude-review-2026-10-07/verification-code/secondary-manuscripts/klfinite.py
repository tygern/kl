"""Independent exact-integer KL toolkit for finite signed-permutation groups (A, B, D).
Conventions (chosen to match computations/README.md so certificates can be compared):
 words multiply on the right; labels start at 0.
 D: s_0 : (a,b,...) -> (-b,-a,...);  s_i (i>=1): swap one-based positions i,i+1.
 B: s_0 negates first entry; s_i swaps one-based positions i,i+1.
 A: s_i swaps zero-based positions i,i+1 (so A_{n-1} on n letters).
"""
from functools import lru_cache
from collections import deque
from itertools import combinations

def act(typ, w, s):
    p = list(w)
    if typ == 'A':
        p[s], p[s+1] = p[s+1], p[s]
    elif typ == 'B':
        if s == 0: p[0] = -p[0]
        else: p[s-1], p[s] = p[s], p[s-1]
    elif typ == 'D':
        if s == 0: p[0], p[1] = -p[1], -p[0]
        else: p[s-1], p[s] = p[s], p[s-1]
    return tuple(p)

def gens(typ, n):
    return list(range(n)) if typ in 'BD' else list(range(n-1))

def length(typ, w):
    n = len(w)
    inv = sum(1 for i in range(n) for j in range(i+1, n) if w[i] > w[j])
    if typ == 'A': return inv
    nsp = sum(1 for i in range(n) for j in range(i+1, n) if w[i] + w[j] < 0)
    if typ == 'D': return inv + nsp
    neg = sum(1 for a in w if a < 0)
    return inv + nsp + neg

class Group:
    def __init__(self, typ, n):
        self.typ, self.n = typ, n
        self.S = gens(typ, n)
        e = tuple(range(1, n+1))
        self.elems = [e]; seen = {e}
        dq = deque([e])
        while dq:
            w = dq.popleft()
            for s in self.S:
                v = act(typ, w, s)
                if v not in seen:
                    seen.add(v); self.elems.append(v); dq.append(v)
        self.len = {w: length(typ, w) for w in self.elems}
        # sanity: length via BFS distance equals formula
        dist = {e: 0}; dq = deque([e])
        while dq:
            w = dq.popleft()
            for s in self.S:
                v = act(typ, w, s)
                if v not in dist:
                    dist[v] = dist[w]+1; dq.append(v)
        assert all(dist[w] == self.len[w] for w in self.elems), "length formula mismatch"
        self._ideal = {}
        self._kl = {}
        self._R = {}
    def rdesc(self, w):
        return [s for s in self.S if self.len[act(self.typ, w, s)] < self.len[w]]
    def ldesc(self, w):
        wi = self.inv(w)
        return self.rdesc(wi)
    def inv(self, w):
        n = len(w); p = [0]*n
        for i, a in enumerate(w):
            p[abs(a)-1] = (i+1) * (1 if a > 0 else -1)
        return tuple(p)
    def ideal(self, w):
        """Bruhat lower ideal via lifting: [e,w] = [e,ws] U [e,ws]s."""
        if w in self._ideal: return self._ideal[w]
        if self.len[w] == 0:
            res = frozenset([w])
        else:
            s = self.rdesc(w)[0]
            I = self.ideal(act(self.typ, w, s))
            res = frozenset(I | {act(self.typ, y, s) for y in I})
        self._ideal[w] = res
        return res
    def leq(self, x, w):
        return x in self.ideal(w)
    def reduced_word(self, w):
        word = []
        while self.len[w] > 0:
            s = self.rdesc(w)[0]
            word.append(s); w = act(self.typ, w, s)
        return word[::-1]
    def is_fc(self, w, m=None):
        """FC iff no reduced word contains s t s (m=3) or s t s t (m=4). Check via all reduced words (DFS)."""
        # enumerate all reduced words of w
        if m is None: m = self.coxmat()
        @lru_cache(None)
        def words(v):
            if self.len[v] == 0: return [()]
            out = []
            for s in self.rdesc(v):
                for wd in words(act(self.typ, v, s)):
                    out.append(wd + (s,))
            return out
        for wd in words(w):
            for i in range(len(wd)):
                for s, t in [(wd[i], wd[i+1] if i+1 < len(wd) else None)]:
                    if t is None or s == t: continue
                    k = m[s][t]
                    if k >= 3 and i + k <= len(wd):
                        seg = wd[i:i+k]
                        if all(seg[j] == (s if j % 2 == 0 else t) for j in range(k)):
                            return False
        return True
    def coxmat(self):
        e = tuple(range(1, self.n+1))
        m = {}
        for s in self.S:
            m[s] = {}
            for t in self.S:
                if s == t: m[s][t] = 1; continue
                w = e; k = 0
                while True:
                    w = act(self.typ, w, s); w = act(self.typ, w, t); k += 1
                    if w == e: break
                m[s][t] = k
        return m
    # polynomials as tuples of ints ascending
    @staticmethod
    def padd(a, b, shift=0, scale=1):
        r = list(a) + [0]*max(0, len(b)+shift-len(a))
        for k, c in enumerate(b): r[k+shift] += scale*c
        while r and r[-1] == 0: r.pop()
        return tuple(r)
    def P(self, x, w):
        if (x, w) in self._kl: return self._kl[(x, w)]
        if not self.leq(x, w): res = ()
        elif x == w: res = (1,)
        else:
            s = self.rdesc(w)[0]
            ws = act(self.typ, w, s); xs = act(self.typ, x, s)
            c = 1 if self.len[xs] < self.len[x] else 0
            res = self.padd(self.padd((), self.P(xs, ws), 1-c), self.P(x, ws), c)
            d = self.len[w]
            for z in self.ideal(ws):
                if self.leq(x, z) and self.len[act(self.typ, z, s)] < self.len[z]:
                    mu = self.mu(z, ws)
                    if mu:
                        res = self.padd(res, self.P(x, z), (d - self.len[z])//2, -mu)
        self._kl[(x, w)] = res
        return res
    def mu(self, x, w):
        d = self.len[w] - self.len[x]
        if d <= 0 or d % 2 == 0: return 0
        p = self.P(x, w)
        k = (d-1)//2
        return p[k] if k < len(p) else 0
    def R(self, x, w):
        if (x, w) in self._R: return self._R[(x, w)]
        if not self.leq(x, w): res = ()
        elif x == w: res = (1,)
        else:
            s = self.rdesc(w)[0]
            ws = act(self.typ, w, s); xs = act(self.typ, x, s)
            if self.len[xs] < self.len[x]: res = self.R(xs, ws)
            else:
                res = self.padd(self.padd((), self.R(x, ws), 0, -1), self.padd(self.R(x, ws), self.R(xs, ws)), 1)
        self._R[(x, w)] = res
        return res
    def interval(self, x, w):
        return [z for z in self.ideal(w) if self.leq(x, z)]
    def check_reciprocity(self, x, w):
        d = self.len[w] - self.len[x]
        lhs = ()
        for z in self.interval(x, w):
            r = self.R(x, z); p = self.P(z, w)
            prod = ()
            for k, c in enumerate(p): prod = self.padd(prod, r, k, c)
            lhs = self.padd(lhs, prod)
        P = self.P(x, w)
        rhs = [0]*(d+1)
        for k, c in enumerate(P): rhs[d-k] += c
        while rhs and rhs[-1] == 0: rhs.pop()
        return lhs == tuple(rhs)

def poset_canon(G, x, w):
    """Canonical form data for interval: rank vector + sorted multiset of (rank, #down covers, #up covers) - used as invariant; full isomorphism check done separately."""
    I = G.interval(x, w)
    return sorted(G.len[z]-G.len[x] for z in I)

def interval_relations(G, x, w):
    I = sorted(G.interval(x, w), key=lambda z: (G.len[z], z))
    idx = {z: i for i, z in enumerate(I)}
    rel = {(idx[a], idx[b]) for a in I for b in I if G.leq(a, b)}
    return I, rel

def check_bijection(G1, x1, w1, G2, x2, w2, pairs):
    """pairs: list of (src_tuple, tgt_tuple). Verify it is a bijection of intervals preserving and reflecting order."""
    I1 = set(G1.interval(x1, w1)); I2 = set(G2.interval(x2, w2))
    src = [a for a, b in pairs]; tgt = [b for a, b in pairs]
    assert set(src) == I1 and len(set(src)) == len(I1), "source vertices mismatch"
    assert set(tgt) == I2 and len(set(tgt)) == len(I2), "target vertices mismatch"
    f = dict(pairs)
    for a in I1:
        for b in I1:
            if G1.leq(a, b) != G2.leq(f[a], f[b]):
                return False
    return True
