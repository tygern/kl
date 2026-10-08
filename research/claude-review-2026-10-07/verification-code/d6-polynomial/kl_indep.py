#!/usr/bin/env python3
"""Independent referee implementation of equal-parameter KL polynomials.

Written from scratch for the d6-polynomial review dimension.  Nothing is
imported from the repository under review.  Exact Python integers only.

Two group models:
  * SignedD(n): type D_n signed permutations with Gern's conventions
      s_1 : (a,b,...) -> (-b,-a,...)   (code label 0)
      s_i : swap one-based positions i-1,i (i>=2)  (code label i-1)
    so code label k corresponds to Gern's s_{k+1}.
  * Geometric(edges, n): simply-laced Coxeter group in the geometric
    representation; an element is the tuple of images of the simple roots
    (columns, integer coordinates in the simple-root basis).

Bruhat order is computed on the principal lower ideal of a fixed top
element via the subword property / Deodhar's property Z:
  L(w) = L(ws) u L(ws)s   for s a right descent of w.
KL polynomials use the standard recursion (Humphreys 7.11 / KL 1979):
  P_{x,w} = q^{1-c} P_{xs,v} + q^c P_{x,v}
            - sum_{x<=z<v, zs<z} mu(z,v) q^{(l(w)-l(z))/2} P_{x,z},
  v = ws < w,  c = 1 if xs < x else 0.
"""
import sys
import time
from collections import Counter

sys.setrecursionlimit(1000000)


def padd(a, b, shift=0, scale=1):
    """a + scale * q^shift * b, ascending coefficient tuples."""
    r = list(a) + [0] * max(0, len(b) + shift - len(a))
    for i, c in enumerate(b):
        r[i + shift] += scale * c
    while r and r[-1] == 0:
        r.pop()
    return tuple(r)


class SignedD:
    def __init__(self, n):
        self.n = n
        self.gens = list(range(n))

    def identity(self):
        return tuple(range(1, self.n + 1))

    def right(self, w, s):
        a = list(w)
        if s == 0:
            a[0], a[1] = -a[1], -a[0]
        else:
            a[s - 1], a[s] = a[s], a[s - 1]
        return tuple(a)

    def length(self, w):
        # Bjorner-Brenti Prop. 8.2.1 (= Gern Prop. 2.2.2)
        n = len(w)
        total = 0
        for i in range(n):
            for j in range(i + 1, n):
                if w[i] > w[j]:
                    total += 1
                if -w[i] > w[j]:
                    total += 1
        return total

    def is_descent(self, w, s):
        return self.length(self.right(w, s)) < self.length(w)


class Geometric:
    def __init__(self, edges, n):
        self.n = n
        self.gens = list(range(n))
        self.adj = [set() for _ in range(n)]
        for a, b in edges:
            self.adj[a].add(b)
            self.adj[b].add(a)

    def identity(self):
        return tuple(tuple(1 if i == j else 0 for i in range(self.n)) for j in range(self.n))

    def right(self, w, s):
        # (w s)(alpha_j) = w(alpha_j - A_{sj} alpha_s); simply laced: A_{sj} = -1 if adjacent
        ci = w[s]
        out = []
        for j in range(self.n):
            if j == s:
                out.append(tuple(-c for c in ci))
            elif j in self.adj[s]:
                out.append(tuple(w[j][k] + ci[k] for k in range(self.n)))
            else:
                out.append(w[j])
        return tuple(out)

    def is_descent(self, w, s):
        col = w[s]
        assert all(c <= 0 for c in col) or all(c >= 0 for c in col), "root neither positive nor negative"
        return all(c <= 0 for c in col)


def word_product(G, word):
    w = G.identity()
    for s in word:
        w = G.right(w, s)
    return w


class Ideal:
    """Principal lower Bruhat ideal of top = product of `word`, with
    lengths, Bruhat bitsets and memoised KL polynomials."""

    def __init__(self, G, word, name=""):
        self.G = G
        self.name = name
        t0 = time.time()
        # subword property: lower ideal = all subword products
        elems = {G.identity()}
        for s in word:
            elems |= {G.right(x, s) for x in elems}
        self.elems = sorted(elems)
        self.idx = {w: i for i, w in enumerate(self.elems)}
        self.top = self.idx[word_product(G, word)]
        N = len(self.elems)
        self.N = N
        # right multiplication table (index or -1 if outside the ideal)
        self.r = [[self.idx.get(G.right(w, s), -1) for s in G.gens] for w in self.elems]
        self.desc = [[G.is_descent(w, s) for s in G.gens] for w in self.elems]
        # lengths by BFS along ascents inside the ideal
        self.length = [None] * N
        e = self.idx[G.identity()]
        self.length[e] = 0
        frontier = [e]
        while frontier:
            nxt = []
            for i in frontier:
                for s in G.gens:
                    j = self.r[i][s]
                    if j >= 0 and not self.desc[i][s] and self.length[j] is None:
                        self.length[j] = self.length[i] + 1
                        nxt.append(j)
            frontier = nxt
        assert all(l is not None for l in self.length), "ideal not connected by ascents"
        # consistency: every descent lowers length by one and stays in the ideal
        for i in range(N):
            for s in G.gens:
                j = self.r[i][s]
                if self.desc[i][s]:
                    assert j >= 0 and self.length[j] == self.length[i] - 1
                elif j >= 0:
                    assert self.length[j] == self.length[i] + 1
        assert self.length[self.top] == len(word), ("word not reduced", self.length[self.top], len(word))
        # Bruhat lower sets as bitsets, by increasing length (property Z)
        order = sorted(range(N), key=lambda i: self.length[i])
        self.low = [0] * N
        for i in order:
            if self.length[i] == 0:
                self.low[i] = 1 << i
                continue
            s = next(s for s in G.gens if self.desc[i][s])
            base = self.low[self.r[i][s]]
            bits = base
            b = base
            while b:
                lsb = b & -b
                k = lsb.bit_length() - 1
                j = self.r[k][s]
                assert j >= 0
                bits |= 1 << j
                b ^= lsb
            self.low[i] = bits
        self.P = {}
        self.corr = {}
        self.build_time = time.time() - t0

    def leq(self, x, w):
        return (self.low[w] >> x) & 1 == 1

    def members(self, bits):
        out = []
        while bits:
            lsb = bits & -bits
            out.append(lsb.bit_length() - 1)
            bits ^= lsb
        return out

    def first_descent(self, w):
        return next(s for s in self.G.gens if self.desc[w][s])

    def mu(self, z, v):
        d = self.length[v] - self.length[z]
        if d <= 0 or d % 2 == 0:
            return 0
        p = self.kl(z, v)
        k = (d - 1) // 2
        return p[k] if len(p) > k else 0

    def corrections(self, v, s):
        """All z < v with zs < z, odd gap and mu(z,v) != 0."""
        key = (v, s)
        if key not in self.corr:
            out = []
            for z in self.members(self.low[v]):
                if z == v or not self.desc[z][s]:
                    continue
                d = self.length[v] - self.length[z]
                if d % 2 == 0:
                    continue
                m = self.mu(z, v)
                if m:
                    out.append((z, (d + 1) // 2, m))
            self.corr[key] = tuple(out)
        return self.corr[key]

    def kl_with_descent(self, x, w, s):
        assert self.desc[w][s]
        v = self.r[w][s]
        xs = self.r[x][s]
        assert xs >= 0
        c = 1 if self.length[xs] < self.length[x] else 0
        p = padd((), self.kl(xs, v), 1 - c)
        p = padd(p, self.kl(x, v), c)
        for z, e, m in self.corrections(v, s):
            if self.leq(x, z):
                p = padd(p, self.kl(x, z), e, -m)
        return p

    def kl(self, x, w):
        key = (x, w)
        if key in self.P:
            return self.P[key]
        if x == w:
            p = (1,)
        elif not self.leq(x, w):
            p = ()
        else:
            p = self.kl_with_descent(x, w, self.first_descent(w))
            assert p and p[0] == 1 and min(p) >= 0, (x, w, p)
            assert 2 * (len(p) - 1) <= self.length[w] - self.length[x] - 1
        self.P[key] = p
        return p

    def interval(self, x, w):
        return [z for z in self.members(self.low[w]) if self.leq(x, z)]

    def rank_vector(self, x, w):
        c = Counter(self.length[z] - self.length[x] for z in self.interval(x, w))
        return [c[i] for i in range(self.length[w] - self.length[x] + 1)]


def reduced_words(ideal, w):
    """All reduced words of element index w (as code-label tuples)."""
    if ideal.length[w] == 0:
        return {()}
    out = set()
    for s in ideal.G.gens:
        if ideal.desc[w][s]:
            for u in reduced_words(ideal, ideal.r[w][s]):
                out.add(u + (s,))
    return out


def commutes(G, s, t):
    if isinstance(G, Geometric):
        return t not in G.adj[s]
    # SignedD: check on identity
    e = G.identity()
    return G.right(G.right(e, s), t) == G.right(G.right(e, t), s)


def commutation_class(G, word):
    seen = {tuple(word)}
    stack = [tuple(word)]
    while stack:
        u = stack.pop()
        for i in range(len(u) - 1):
            if u[i] != u[i + 1] and commutes(G, u[i], u[i + 1]):
                v = u[:i] + (u[i + 1], u[i]) + u[i + 2:]
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
    return seen


def is_fully_commutative(ideal, w):
    """Stembridge: w is FC iff its reduced words form one commutation class."""
    words = reduced_words(ideal, w)
    cls = commutation_class(ideal.G, next(iter(words)))
    return words == cls, len(words)


def show(p):
    terms = []
    for i, c in enumerate(p):
        if c:
            terms.append(f"{c}" if i == 0 else (f"{c}q" if i == 1 else f"{c}q^{i}"))
    return " + ".join(terms) if terms else "0"
