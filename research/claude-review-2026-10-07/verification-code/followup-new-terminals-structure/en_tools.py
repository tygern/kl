"""Independent exact tools for generalized E_n (paper labelling) and Gern's D_m.

Paper labelling: chain 0-1-...-(n-2), node n-1 attached to node 2.
Elements are integer matrices M with columns M[:,j] = w(alpha_j) in the
simple-root basis.  Words act by right multiplication: word "abc" = s_a s_b s_c.
"""
from fractions import Fraction
from itertools import combinations
import sys

sys.setrecursionlimit(10000)


def en_edges(n):
    return [(i, i + 1) for i in range(n - 2)] + [(2, n - 1)]


class RootSystem:
    def __init__(self, n, edges):
        self.n = n
        self.edges = edges
        self.adj = {i: set() for i in range(n)}
        for a, b in edges:
            self.adj[a].add(b)
            self.adj[b].add(a)
        # Cartan matrix
        self.C = [[2 if i == j else (-1 if j in self.adj[i] else 0)
                   for j in range(n)] for i in range(n)]
        self.pos_roots = None

    def pairing(self, u, v):
        return sum(u[i] * self.C[i][j] * v[j] for i in range(self.n) for j in range(self.n))

    def simple_reflect(self, s, v):
        """s_s(v) = v - (alpha_s, v) alpha_s."""
        c = sum(self.C[s][j] * v[j] for j in range(self.n))
        w = list(v)
        w[s] -= c
        return tuple(w)

    def positive_roots(self):
        if self.pos_roots is None:
            simple = [tuple(int(i == j) for i in range(self.n)) for j in range(self.n)]
            seen = set(simple)
            frontier = list(simple)
            while frontier:
                new = []
                for r in frontier:
                    for s in range(self.n):
                        q = self.simple_reflect(s, r)
                        if all(c >= 0 for c in q) and q not in seen:
                            seen.add(q)
                            new.append(q)
                frontier = new
            self.pos_roots = sorted(seen, key=lambda r: (sum(r), r))
        return self.pos_roots

    # ---- matrices ----
    def identity(self):
        return tuple(tuple(int(i == j) for i in range(self.n)) for j in range(self.n))

    def col(self, M, j):
        return M[j]

    def right_mult(self, M, s):
        """M * s_s : columns change: a_s -> -a_s, a_t -> a_t + a_s for t ~ s."""
        cols = [list(c) for c in M]
        old = cols[s][:]
        cols[s] = [-v for v in old]
        for t in self.adj[s]:
            cols[t] = [u + v for u, v in zip(cols[t], old)]
        return tuple(tuple(c) for c in cols)

    def left_mult(self, s, M):
        """s_s * M : apply s_s to every column."""
        return tuple(self.simple_reflect(s, c) for c in M)

    def word_matrix(self, word):
        M = self.identity()
        for s in word:
            M = self.right_mult(M, s)
        return M

    def apply(self, M, v):
        out = [0] * self.n
        for j in range(self.n):
            if v[j]:
                for i in range(self.n):
                    out[i] += M[j][i] * v[j]
        return tuple(out)

    def mult(self, A, B):
        """(A*B)(alpha_j) = A(B(alpha_j))."""
        return tuple(self.apply(A, B[j]) for j in range(self.n))

    def inverse(self, M):
        """Inverse via orthogonality: M^{-1} = C^{-1} M^T C (exact rationals)."""
        n = self.n
        # Solve: we want N with N(alpha_j) = w^{-1}(alpha_j). Use Fractions.
        # Build matrix of M in standard (rows=coords, cols=j): Mat[i][j] = M[j][i]
        Mat = [[Fraction(M[j][i]) for j in range(n)] for i in range(n)]
        # Gaussian elimination inverse
        aug = [row[:] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(Mat)]
        for c in range(n):
            p = next(r for r in range(c, n) if aug[r][c] != 0)
            aug[c], aug[p] = aug[p], aug[c]
            piv = aug[c][c]
            aug[c] = [x / piv for x in aug[c]]
            for r in range(n):
                if r != c and aug[r][c] != 0:
                    f = aug[r][c]
                    aug[r] = [x - f * y for x, y in zip(aug[r], aug[c])]
        inv = [[aug[i][n + j] for j in range(n)] for i in range(n)]
        for row in inv:
            for x in row:
                assert x.denominator == 1
        return tuple(tuple(int(inv[i][j]) for i in range(n)) for j in range(n))

    def is_negative(self, v):
        return all(c <= 0 for c in v) and any(c < 0 for c in v)

    def is_positive(self, v):
        return all(c >= 0 for c in v) and any(c > 0 for c in v)

    def right_descents(self, M):
        return {s for s in range(self.n) if self.is_negative(M[s])}

    def left_descents(self, M):
        return self.right_descents(self.inverse(M))

    def length(self, M):
        return sum(1 for r in self.positive_roots() if self.is_negative(self.apply(M, r)))

    def inversions(self, M):
        return [r for r in self.positive_roots() if self.is_negative(self.apply(M, r))]

    def right_terminal(self, M):
        R = self.right_descents(M)
        for s in R:
            for t in self.adj[s]:
                v = tuple(a + b for a, b in zip(M[s], M[t]))
                if not self.is_positive(v):
                    return False
        return True

    def terminal(self, M):
        return self.right_terminal(M) and self.right_terminal(self.inverse(M))

    def is_commuting_product(self, M):
        """w is a product of commuting generators iff w = i(R(w)) with R(w) independent."""
        R = self.right_descents(M)
        if any(t in R for s in R for t in self.adj[s]):
            return False
        return self.word_matrix(sorted(R)) == M

    def reduced_word(self, M):
        """Greedy reduced word via right descents."""
        word = []
        Mi = M
        while True:
            R = self.right_descents(Mi)
            if not R:
                break
            s = min(R)
            word.append(s)
            Mi = self.right_mult(Mi, s)
        return word[::-1]

    def support(self, M):
        return set(self.reduced_word(M))

    def longest_element(self, J):
        """w_0(W_J) by greedy ascent inside W_J."""
        M = self.identity()
        while True:
            asc = [s for s in J if not self.is_negative(M[s])]
            if not asc:
                return M
            M = self.right_mult(M, asc[0])


def word_str(word):
    return ''.join(str(s) for s in word)


def parse(wstr):
    return [int(c) for c in wstr]


# ---------- Gern's D_m in his labels 1..m ----------
# Gern: nodes 1,2 attached to 3; chain 3-4-...-m.  Signed permutations:
# s1 = (1,-2)(-1,2), s_i = (i-1,i) for i>=2.

def gern_interval(i, j):
    """Gern Definition 2.3.1, returns a list of generator labels."""
    if 0 <= j < i and i >= 2:
        return gern_interval(j, i)[::-1]
    if i == 1 and j >= 3:
        return [1] + list(range(3, j + 1))
    if i == 0 and j >= 2:
        return list(range(1, j + 1))
    if 2 <= i <= j:
        return list(range(i, j + 1))
    raise ValueError((i, j))


def gern_wn_word(n):
    """Lemma 2.3.4 bracket form (even n; for odd n use n-1)."""
    if n % 2 == 1:
        return gern_wn_word(n - 1)
    k = n // 2 - 2
    word = []
    for i in range(2, n + 1, 2):
        word += gern_interval(i, 0)
    # tail [n-k, n-2k] [n-k+1, n-2k+2] ... [n-1, n-2] [n, n]
    for i in range(0, k + 1):
        word += gern_interval(n - k + i, n - 2 * k + 2 * i)
    return word


def gern_to_en(label, n_en):
    """Gern D_{n_en-1} label -> paper E_{n_en} label (Section 3 of the paper)."""
    if label == 1:
        return 1
    if label == 2:
        return n_en - 1
    return label - 1


def signed_perm_of_gern_word(word, m):
    """One-line signed permutation of a Gern word, composing left to right
    (word = s_a s_b ... acting on the left: w(i) = s_a(s_b(...(i)))."""
    def gen(label):
        p = list(range(1, m + 1))
        if label == 1:
            p[0], p[1] = -2, -1
        else:
            p[label - 2], p[label - 1] = label, label - 1
        return p

    def compose(u, v):  # u after v
        return [(1 if a > 0 else -1) * u[abs(a) - 1] for a in v]

    w = list(range(1, m + 1))
    for s in word:
        w = compose(w, gen(s))
    return tuple(w)
