"""Exact crystallographic word surgery experiments; no third-party packages.

Faithful geometric representation in the simple root basis.  Generalized
Cartan entries use products 0,1,2,3,4 for braid labels 2,3,4,6,infinity.
All interval comparisons are by the standard lifting recursion.
"""
from functools import lru_cache
from collections import Counter
import json


class Coxeter:
    def __init__(self, n, edges):
        self.n = n
        self.cartan = [[2 * (i == j) for j in range(n)] for i in range(n)]
        for i, j, p in edges:
            self.cartan[i][j] = -p
            self.cartan[j][i] = -1
        self.e = tuple(tuple(int(i == j) for i in range(n)) for j in range(n))

    @lru_cache(None)
    def right(self, a, s):
        root = a[s]
        return tuple(tuple(col[i] - self.cartan[s][j] * root[i]
                           for i in range(self.n)) for j, col in enumerate(a))

    def positive(self, a, s):
        return all(c >= 0 for c in a[s])

    @lru_cache(None)
    def word(self, a):
        if a == self.e:
            return ()
        for s in range(self.n):
            if not self.positive(a, s):
                return self.word(self.right(a, s)) + (s,)
        raise AssertionError(a)

    def length(self, a):
        return len(self.word(a))

    def element(self, word):
        a = self.e
        for s in word:
            a = self.right(a, s)
        return a

    @lru_cache(None)
    def le(self, a, b):
        if a == self.e:
            return True
        if self.length(a) > self.length(b):
            return False
        if self.length(a) == self.length(b):
            return a == b
        s = self.word(b)[-1]
        if not self.positive(a, s):
            a = self.right(a, s)
        return self.le(a, self.right(b, s))

    def lower(self, a):
        vertices = {self.e}
        for s in self.word(a):
            vertices |= {self.right(x, s) for x in vertices}
        return vertices

    def interval(self, a, b):
        return {x for x in self.lower(b) if self.le(a, x)}

    def ranks(self, vertices):
        count = Counter(map(self.length, vertices))
        if not count:
            return []
        return [count[i] for i in range(min(count), max(count) + 1)]


def d_edges(n):
    return [(0, 2, 1), (1, 2, 1)] + [(i, i + 1, 1) for i in range(2, n - 1)]


def bad_word(n):
    # Independent signed-permutation reduction supplies the fixed D_n word.
    identity = tuple(range(1, n + 1))
    top = tuple(((-1) ** (n // 2)) if i == 1 else
                i if i % 2 else -(n + 2 - i) for i in range(1, n + 1))
    def length(a):
        return sum(a[i] > a[j] for i in range(n) for j in range(i+1,n)) + sum(
            -a[i] > a[j] for i in range(n) for j in range(i+1,n))
    reverse = []
    while top != identity:
        for s in range(n):
            b = list(top)
            if s == 0:
                b[0], b[1] = -b[1], -b[0]
            else:
                b[s-1], b[s] = b[s], b[s-1]
            b = tuple(b)
            if length(b) < length(top):
                top = b
                reverse.append(s)
                break
    return tuple(reversed(reverse))


def test_inflation(n):
    word = bad_word(n)
    bottom = (0, 1) + tuple(range(3, n, 2))
    for edge in [None] + list(range(n - 1)):
        edges = d_edges(n)
        if edge is not None:
            i,j,_ = edges[edge]
            edges[edge] = (i,j,2)
        group = Coxeter(n, edges)
        top = group.element(word)
        lo = group.element(bottom)
        vertices = group.interval(lo, top)
        print(json.dumps({"n":n,"inflated":edge,"word":word,
                          "rank_sizes":group.ranks(vertices),
                          "size":len(vertices)}), flush=True)


if __name__ == "__main__":
    test_inflation(4)
    test_inflation(6)
