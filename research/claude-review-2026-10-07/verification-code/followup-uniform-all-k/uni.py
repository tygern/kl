"""Independent toolkit for generalized E_n (manuscript numbering: chain 0-1-...-(n-2),
node n-1 attached to node 2).  Elements are integer matrices stored as a tuple of n
columns, column j = w(alpha_j) in the simple-root basis.  Exact integer arithmetic only.
Nothing from the repository or from other scratch directories is imported."""


def adjacency(n):
    adj = [[] for _ in range(n)]
    for i in range(n - 2):
        adj[i].append(i + 1)
        adj[i + 1].append(i)
    adj[n - 1].append(2)
    adj[2].append(n - 1)
    return adj


class En:
    def __init__(self, n):
        self.n = n
        self.adj = adjacency(n)
        self.adjset = [set(a) for a in self.adj]
        self.e = tuple(tuple(1 if i == j else 0 for i in range(n)) for j in range(n))

    # bilinear form: (alpha_i, alpha_i) = 2, adjacent -1
    def pair(self, v, i):
        return 2 * v[i] - sum(v[j] for j in self.adj[i])

    def form(self, v, w):
        return sum(v[i] * self.pair(w, i) for i in range(self.n))

    def refl_vec(self, v, i):
        c = self.pair(v, i)
        if c == 0:
            return v
        v = list(v)
        v[i] -= c
        return tuple(v)

    def rmul(self, w, i):
        cols = list(w)
        ai = cols[i]
        for j in self.adj[i]:
            cols[j] = tuple(x + y for x, y in zip(cols[j], ai))
        cols[i] = tuple(-x for x in ai)
        return tuple(cols)

    def lmul(self, i, w):
        return tuple(self.refl_vec(col, i) for col in w)

    def mul(self, u, v):
        # (uv)(alpha_j) = u(v(alpha_j)) = sum_i v_ij u(alpha_i)
        n = self.n
        out = []
        for j in range(n):
            col = [0] * n
            for i in range(n):
                c = v[j][i]
                if c:
                    ui = u[i]
                    for k in range(n):
                        col[k] += c * ui[k]
            out.append(tuple(col))
        return tuple(out)

    @staticmethod
    def neg(col):
        return any(x < 0 for x in col) and all(x <= 0 for x in col)

    @staticmethod
    def pos(col):
        return any(x > 0 for x in col) and all(x >= 0 for x in col)

    def R(self, w):
        return [i for i in range(self.n) if self.neg(w[i])]

    def word(self, w):
        """A reduced word (left-to-right product) found by right descent stripping."""
        out = []
        while True:
            d = -1
            for i in range(self.n):
                if self.neg(w[i]):
                    d = i
                    break
            if d < 0:
                break
            out.append(d)
            w = self.rmul(w, d)
        out.reverse()
        return out

    def length(self, w):
        return len(self.word(w))

    def inverse(self, w):
        v = self.e
        for i in reversed(self.word(w)):
            v = self.rmul(v, i)
        return v

    def L(self, w):
        return self.R(self.inverse(w))

    def from_word(self, word):
        w = self.e
        for i in word:
            w = self.rmul(w, i)
        return w

    def reflection(self, beta):
        """r_beta(alpha_j) = alpha_j - (alpha_j,beta) beta."""
        n = self.n
        cols = []
        for j in range(n):
            m = self.pair(beta, j)
            cols.append(tuple((1 if i == j else 0) - m * beta[i] for i in range(n)))
        return tuple(cols)

    def right_terminal(self, w):
        """No reduced word ends in a noncommuting pair: for s in R(w) and t ~ s, w(alpha_s+alpha_t) > 0."""
        for s in self.R(w):
            for t in self.adj[s]:
                if not self.pos(tuple(a + b for a, b in zip(w[s], w[t]))):
                    return False
        return True

    def terminal(self, w):
        return self.right_terminal(w) and self.right_terminal(self.inverse(w))

    def is_fc_word(self, word):
        """Stembridge heap test (simply laced): a reduced word is FC iff between any two
        consecutive occurrences of a generator s there are at least two occurrences of
        neighbours of s (all neighbour occurrences between them in the word are between them
        in the heap, since adjacent labels are comparable)."""
        last = {}
        for p, s in enumerate(word):
            if s in last:
                cnt = sum(1 for t in word[last[s] + 1:p] if t in self.adjset[s])
                if cnt < 2:
                    return False
            last[s] = p
        return True

    def covers(self, w):
        """All Bruhat covers x < w with l(x) = l(w) - 1 (single deletions from one reduced word)."""
        word = self.word(w)
        l = len(word)
        seen = {}
        for p in range(l):
            sub = word[:p] + word[p + 1:]
            v = self.from_word(sub)
            if v in seen:
                continue
            if self.length(v) == l - 1:
                seen[v] = sub
        return seen

    def depth(self, beta):
        """Depth: least number of simple reflections lowering beta (positive root) to a simple
        root, by greedy BFS on heights (exact BFS over the lowering tree)."""
        from collections import deque
        start = tuple(beta)
        dist = {start: 0}
        dq = deque([start])
        while dq:
            v = dq.popleft()
            if sum(v) == 1:
                return dist[v]
            for i in range(self.n):
                c = self.pair(v, i)
                if c > 0:
                    u = list(v)
                    u[i] -= c
                    u = tuple(u)
                    if u not in dist:
                        dist[u] = dist[v] + 1
                        dq.append(u)
        return None


def uniform_data(r):
    n = 4 * r + 1
    beta = [0] * n
    beta[0], beta[1], beta[2], beta[4 * r] = r - 1, r, 2 * r - 1, r
    for j in range(3, 4 * r):
        beta[j] = 2 * r - j // 2
    delta = [0] * n
    gamma = [0] * n
    for j, v in enumerate((2, 4, 6, 5, 4, 3, 2, 1)):
        delta[j] = v
    delta[4 * r] = 3
    for j, v in enumerate((1, 2, 3, 2, 2, 1, 1, 0)):
        gamma[j] = v
    gamma[4 * r] = 1
    return n, tuple(beta), tuple(gamma), tuple(delta), 2 * r - 4


def beta_rk(r, k):
    n, beta, gamma, delta, a = uniform_data(r)
    return tuple(beta[j] - a * k * gamma[j] + (a * k * k + r * k) * delta[j] for j in range(n))


def I_r(r):
    return sorted({0, 4 * r} | set(range(3, 4 * r, 2)))
