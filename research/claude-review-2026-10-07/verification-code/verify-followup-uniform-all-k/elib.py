"""Exact-integer tools for generalized E_n Coxeter groups (simple-root basis).

Numbering (manuscript): chain 0-1-...-(n-2), node n-1 attached to node 2.
Group elements are integer matrices (tuple of tuples), column j = image of alpha_j.
"""
from functools import lru_cache
import sys

sys.setrecursionlimit(10000)


def cartan(n):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        A[i][i] = 2
    for i in range(n - 2):
        A[i][i + 1] = A[i + 1][i] = -1
    A[2][n - 1] = A[n - 1][2] = -1
    return A


def adjacent(A, s, t):
    return s != t and A[s][t] != 0


def pairing(A, u, v):
    n = len(A)
    return sum(u[i] * A[i][j] * v[j] for i in range(n) for j in range(n) if u[i] and v[j])


def simple_refl(A, s, v):
    """s_s(v) = v - (alpha_s, v) alpha_s"""
    c = sum(A[s][j] * v[j] for j in range(len(A)))
    w = list(v)
    w[s] -= c
    return tuple(w)


def reflection_matrix(A, beta):
    """r_beta(alpha_j) = alpha_j - (alpha_j, beta) beta. Returns matrix as tuple of columns."""
    n = len(A)
    cols = []
    for j in range(n):
        c = sum(A[j][i] * beta[i] for i in range(n))
        col = [-c * beta[i] for i in range(n)]
        col[j] += 1
        cols.append(tuple(col))
    return tuple(cols)  # cols[j] = image of alpha_j


def identity(n):
    return tuple(tuple(1 if i == j else 0 for i in range(n)) for j in range(n))


def mult_right_simple(A, M, s):
    """(M s_s)(alpha_j) = M(s_s alpha_j) = M(alpha_j - A[s][j] alpha_s) = M alpha_j - A[s][j] M alpha_s."""
    n = len(A)
    ms = M[s]
    cols = []
    for j in range(n):
        if j == s:
            cols.append(tuple(-x for x in ms))
        elif A[s][j] != 0:
            a = A[s][j]
            cols.append(tuple(M[j][i] - a * ms[i] for i in range(n)))
        else:
            cols.append(M[j])
    return tuple(cols)


def mult_left_simple(A, s, M):
    """(s_s M)(alpha_j) = s_s(M alpha_j)."""
    return tuple(simple_refl(A, s, M[j]) for j in range(len(A)))


def is_negative(v):
    return all(x <= 0 for x in v) and any(x < 0 for x in v)


def is_positive(v):
    return all(x >= 0 for x in v) and any(x > 0 for x in v)


def right_descents(M):
    return [s for s in range(len(M)) if is_negative(M[s])]


def reduced_word(A, M):
    """Return a reduced word (list of generators) for M by stripping right descents.
    If M = s_{i1}...s_{im} then stripping finds i_m first; we return [i1,...,im]."""
    n = len(A)
    word = []
    I = identity(n)
    cur = M
    while cur != I:
        D = right_descents(cur)
        assert D, "no right descent but not identity: not a Coxeter group element?"
        s = D[0]
        word.append(s)
        cur = mult_right_simple(A, cur, s)
    word.reverse()
    return word


def word_to_matrix(A, word):
    M = identity(len(A))
    for s in word:
        M = mult_right_simple(A, M, s)
    return M


def inverse(A, M):
    w = reduced_word(A, M)
    return word_to_matrix(A, list(reversed(w)))


def left_descents(A, M):
    return right_descents(inverse(A, M))


def length(A, M):
    return len(reduced_word(A, M))


def is_fc_word(A, word):
    """Heap test (Stembridge): a reduced word is a reduced word of an FC element iff
    the heap has no convex chain i<j<k with labels s,t,s. For a reduced word in a
    simply laced group this is: no two occurrences i<k of the same letter whose open
    heap interval (i,k) consists of exactly one element."""
    m = len(word)
    # above[i] = bitmask of j>i with i<j in heap order (transitive closure of noncommuting precedence)
    above = [0] * m
    for i in range(m - 1, -1, -1):
        mask = 0
        for j in range(i + 1, m):
            if word[j] == word[i] or adjacent(A, word[i], word[j]):
                mask |= (1 << j) | above[j]
        above[i] = mask
    for i in range(m):
        for k in range(i + 1, m):
            if word[k] == word[i] and (above[i] >> k) & 1:
                # open interval = {x : i<x<k} = above[i] & below[k]
                # below[k]: compute via above: x<k iff (above[x]>>k)&1
                cnt = 0
                mask = above[i]
                x = 0
                while mask:
                    if mask & 1 and x < k and (above[x] >> k) & 1:
                        cnt += 1
                        if cnt > 1:
                            break
                    mask >>= 1
                    x += 1
                if cnt == 1:
                    return False
    return True


def is_fc(A, M):
    return is_fc_word(A, reduced_word(A, M))


def depth(A, beta):
    """dp(beta) = min{ l(w) : w beta < 0 }, computed by greedy lowering
    (dp(s beta) = dp(beta)-1 when (beta, alpha_s) > 0)."""
    n = len(A)
    d = 1
    b = tuple(beta)
    while True:
        if sum(b) == 1 and all(x in (0, 1) for x in b):
            return d
        for s in range(n):
            if sum(A[s][j] * b[j] for j in range(n)) > 0:
                b = simple_refl(A, s, b)
                d += 1
                break
        else:
            raise ValueError("no lowering generator; not a real root?")


def uniform_root(r, k):
    """beta_{r,k} of the manuscript, in E_{4r+1}, as a tuple."""
    n = 4 * r + 1
    beta = [0] * n
    beta[0], beta[1], beta[2], beta[n - 1] = r - 1, r, 2 * r - 1, r
    for j in range(3, 4 * r):
        beta[j] = 2 * r - (j // 2)
    delta = [0] * n
    gamma = [0] * n
    for j, v in enumerate((2, 4, 6, 5, 4, 3, 2, 1)):
        delta[j] = v
    delta[n - 1] = 3
    for j, v in enumerate((1, 2, 3, 2, 2, 1, 1, 0)):
        gamma[j] = v
    gamma[n - 1] = 1
    a = 2 * r - 4
    return tuple(beta[i] - a * k * gamma[i] + (a * k * k + r * k) * delta[i] for i in range(n))


def I_r(r):
    return sorted({0} | set(range(3, 4 * r, 2)) | {4 * r})
