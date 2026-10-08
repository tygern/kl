"""Independent exact-integer Coxeter geometric representation utilities.

Vectors are tuples of ints in the simple-root basis.  A group element is
stored as a tuple of n column-vectors: col[j] = w(alpha_j).
Simply laced only; (alpha_i,alpha_i)=2, (alpha_i,alpha_j)=-1 for edges.
"""
from fractions import Fraction


def cartan(n, edges):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    for (i, j) in edges:
        A[i][j] = -1
        A[j][i] = -1
    return A


def pair(A, u, v):
    n = len(A)
    return sum(u[i] * A[i][j] * v[j] for i in range(n) for j in range(n))


def apply(cols, v):
    """Apply element (columns) to vector v."""
    n = len(cols)
    out = [0] * n
    for j in range(n):
        if v[j]:
            c = cols[j]
            for i in range(n):
                out[i] += v[j] * c[i]
    return tuple(out)


def compose(a, b):
    """Return the element a*b (first apply b, then a)."""
    return tuple(apply(a, b[j]) for j in range(len(b)))


def identity(n):
    return tuple(tuple(1 if i == j else 0 for i in range(n)) for j in range(n))


def simple_refl(A, s):
    n = len(A)
    cols = []
    for j in range(n):
        v = [1 if i == j else 0 for i in range(n)]
        # s_s(alpha_j) = alpha_j - A[s][j] alpha_s
        v[s] -= A[s][j]
        cols.append(tuple(v))
    return tuple(cols)


def reflection(A, beta):
    """r_beta(v) = v - (v,beta) beta ; requires (beta,beta)=2."""
    assert pair(A, beta, beta) == 2, pair(A, beta, beta)
    n = len(A)
    cols = []
    for j in range(n):
        e = tuple(1 if i == j else 0 for i in range(n))
        p = pair(A, e, beta)
        cols.append(tuple(e[i] - p * beta[i] for i in range(n)))
    return tuple(cols)


def word_product(A, word, left_to_right=True):
    """Product of simple reflections.  left_to_right: s_{w0} s_{w1} ... as
    an operator product (so the rightmost letter acts first)."""
    n = len(A)
    gens = [simple_refl(A, s) for s in range(n)]
    w = identity(n)
    seq = word if left_to_right else list(reversed(word))
    for s in seq:
        w = compose(w, gens[s])
    return w


def is_negative(v):
    return all(x <= 0 for x in v) and any(x < 0 for x in v)


def is_positive(v):
    return all(x >= 0 for x in v) and any(x > 0 for x in v)


def right_descents(cols):
    return frozenset(j for j in range(len(cols)) if is_negative(cols[j]))


def inverse(cols):
    """Inverse via exact rational Gaussian elimination (works for any
    invertible integer matrix; result must be integral for Weyl elements)."""
    n = len(cols)
    M = [[Fraction(cols[j][i]) for j in range(n)] + [Fraction(1 if i == k else 0) for k in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    inv = [[M[i][n + j] for j in range(n)] for i in range(n)]
    for row in inv:
        for x in row:
            assert x.denominator == 1
    return tuple(tuple(int(inv[i][j]) for i in range(n)) for j in range(n))


def left_descents(A, cols):
    return right_descents(inverse(cols))


def length_by_stripping(A, cols, max_steps=10**6):
    """Length via repeated right-descent removal; also returns a reduced word
    (read as letters w = s_{a_1} ... s_{a_m})."""
    n = len(A)
    gens = [simple_refl(A, s) for s in range(n)]
    w = cols
    word_rev = []
    steps = 0
    while True:
        D = right_descents(w)
        if not D:
            break
        s = min(D)
        w = compose(w, gens[s])
        word_rev.append(s)
        steps += 1
        assert steps < max_steps
    assert w == identity(n), "did not reach identity: not a Weyl group element?"
    return len(word_rev), list(reversed(word_rev))


def is_reduced_word(A, word):
    """Word is reduced iff each successive prefix gains length."""
    n = len(A)
    gens = [simple_refl(A, s) for s in range(n)]
    w = identity(n)
    for s in word:
        # w*s longer than w iff w(alpha_s) > 0
        if not is_positive(w[s]):
            return False
        w = compose(w, gens[s])
    return True


def height(v):
    return sum(v)
