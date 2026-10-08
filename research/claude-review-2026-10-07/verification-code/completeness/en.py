"""Exact-integer geometric representation of generalized E_n (paper numbering):
chain 0-1-...-(n-2), node n-1 attached to node 2.  Elements are stored as a tuple of
n columns, column j = w(alpha_j) in the simple-root basis.  No repo code is imported."""
import sys

def adjacency(n):
    adj = [set() for _ in range(n)]
    for i in range(n - 2):
        adj[i].add(i + 1); adj[i + 1].add(i)
    adj[n - 1].add(2); adj[2].add(n - 1)
    return adj

class En:
    def __init__(self, n):
        self.n = n
        self.adj = adjacency(n)
        self.e = tuple(tuple(1 if i == j else 0 for i in range(n)) for j in range(n))

    def pair(self, v, i):
        """(v, alpha_i) with (alpha_i,alpha_i)=2, adjacent -1."""
        return 2 * v[i] - sum(v[j] for j in self.adj[i])

    def refl_simple(self, v, i):
        c = self.pair(v, i)
        if c == 0:
            return v
        v = list(v); v[i] -= c
        return tuple(v)

    def rmul(self, w, i):
        """w * s_i: a_i -> -a_i, a_j -> a_j + a_i for j ~ i."""
        cols = list(w)
        ai = cols[i]
        for j in self.adj[i]:
            cols[j] = tuple(x + y for x, y in zip(cols[j], ai))
        cols[i] = tuple(-x for x in ai)
        return tuple(cols)

    def lmul(self, i, w):
        return tuple(self.refl_simple(col, i) for col in w)

    @staticmethod
    def neg(col):
        return any(x < 0 for x in col) and all(x <= 0 for x in col)

    @staticmethod
    def pos(col):
        return any(x > 0 for x in col) and all(x >= 0 for x in col)

    def R(self, w):
        return [i for i in range(self.n) if self.neg(w[i])]

    def inverse(self, w):
        # w^{-1} via descent stripping: w = s_{i1}...s_{ik}; w^{-1} = s_{ik}...s_{i1}
        word = self.word(w)
        v = self.e
        for i in reversed(word):
            v = self.rmul(v, i)
        return v

    def L(self, w):
        return self.R(self.inverse(w))

    def word(self, w):
        """A reduced word (left-to-right product) by right descent stripping."""
        out = []
        while True:
            R = self.R(w)
            if not R:
                break
            i = R[0]
            out.append(i)
            w = self.rmul(w, i)
        return list(reversed(out))

    def length(self, w):
        return len(self.word(w))

    def from_word(self, word):
        w = self.e
        for i in word:
            w = self.rmul(w, i)
        return w

    def commuting(self, I):
        return all(j not in self.adj[i] for i in I for j in I)

    def iI(self, I):
        return self.from_word(sorted(I))

    def is_fc_word(self, word):
        """Stembridge: a reduced word of a simply-laced FC element has, between any two
        consecutive occurrences of s, at least two letters not commuting with s."""
        last = {}
        n = self.n
        for pos, s in enumerate(word):
            if s in last:
                between = word[last[s] + 1:pos]
                k = sum(1 for t in between if t in self.adj[s])
                if k < 2:
                    return False
            last[s] = pos
        return True

    def right_terminal(self, w):
        for s in self.R(w):
            for t in self.adj[s]:
                col = tuple(x + y for x, y in zip(w[s], w[t]))
                if not self.pos(col):
                    return False
        return True

    def terminal(self, w):
        return self.right_terminal(w) and self.right_terminal(self.inverse(w))

    def fc_catalogue(self, maxlen=None):
        """Level-by-level closure under right ascents with the recurrence test
        w in FC  <=>  R(w) commuting and w s in FC for all s in R(w)."""
        levels = [{self.e}]
        allfc = {self.e}
        while True:
            prev = levels[-1]
            prevprev = levels[-2] if len(levels) >= 2 else set()
            nxt = set()
            for u in prev:
                Ru = set(self.R(u))
                for s in range(self.n):
                    if s in Ru:
                        continue
                    w = self.rmul(u, s)
                    if w in nxt:
                        continue
                    Rw = self.R(w)
                    if not self.commuting(Rw):
                        continue
                    ok = True
                    for t in Rw:
                        if self.rmul(w, t) not in prev:
                            ok = False; break
                    if ok:
                        nxt.add(w)
            if not nxt:
                break
            levels.append(nxt)
            allfc |= nxt
            if maxlen is not None and len(levels) - 1 >= maxlen:
                break
        return levels, allfc
