#!/usr/bin/env python3
"""Independent reproduction: E9 (affine E8) reflection in root (1,2,3,3,2,2,1,1,2).

Numbering per manuscript Section 5: chain 0-1-2-3-4-5-6-7, node 8 attached to node 2.
Simple roots alpha_j normalised (a_j,a_j)=2, adjacent (a_s,a_t)=-1.
All arithmetic exact integers. Elements represented by their matrix in the simple-root
basis (columns = images of simple roots).
"""
import sys
from collections import deque

n = 9
edges = [(i, i + 1) for i in range(7)] + [(2, 8)]
adj = [set() for _ in range(n)]
for a, b in edges:
    adj[a].add(b); adj[b].add(a)
# Cartan / bilinear form
A = [[2 if i == j else (-1 if j in adj[i] else 0) for j in range(n)] for i in range(n)]

def form(u, v):
    return sum(u[i] * A[i][j] * v[j] for i in range(n) for j in range(n))

def matmul(X, Y):
    return tuple(tuple(sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)) for i in range(n))

def apply(X, v):
    return tuple(sum(X[i][k] * v[k] for k in range(n)) for i in range(n))

ID = tuple(tuple(1 if i == j else 0 for j in range(n)) for i in range(n))
# simple reflection matrices: s_j(alpha_k) = alpha_k - A[j][k] alpha_j  (column k)
def simple(j):
    M = [[1 if i == k else 0 for k in range(n)] for i in range(n)]
    for k in range(n):
        M[j][k] -= A[j][k]
    return tuple(tuple(r) for r in M)
S = [simple(j) for j in range(n)]

def col(X, j):
    return tuple(X[i][j] for i in range(n))

def is_neg(v):
    return all(c <= 0 for c in v) and any(c < 0 for c in v)

def is_pos(v):
    return all(c >= 0 for c in v) and any(c > 0 for c in v)

def right_descents(X):
    return frozenset(j for j in range(n) if is_neg(col(X, j)))

def transpose_inv(X):
    # inverse via reduced word; cheaper: compute reduced word then reverse
    w = reduced_word(X)
    Y = ID
    for s in reversed(w):
        Y = matmul(Y, S[s])
    return Y

def reduced_word(X):
    """Right-strip descents: returns word w1...wk with X = S[w1]...S[wk]."""
    word = []
    Y = X
    while True:
        D = right_descents(Y)
        if not D:
            break
        s = min(D)
        Y = matmul(Y, S[s])
        word.append(s)
    assert Y == ID, "did not reach identity"
    return list(reversed(word))

def left_descents(X):
    return right_descents(transpose_inv(X))

def terminal_right(X):
    """Manuscript eq. (terminaltest): s in R(w) => w(alpha_s+alpha_t) > 0 for all t ~ s."""
    for s in right_descents(X):
        for t in adj[s]:
            v = tuple(X[i][s] + X[i][t] for i in range(n))
            if not is_pos(v):
                return False
    return True

def reflection_matrix(beta):
    assert form(beta, beta) == 2, "not a real root of norm 2"
    cols = []
    for j in range(n):
        e = tuple(1 if i == j else 0 for i in range(n))
        c = form(e, beta)
        cols.append(tuple(e[i] - c * beta[i] for i in range(n)))
    return tuple(tuple(cols[j][i] for j in range(n)) for i in range(n))

def bruhat_leq(x, w):
    """Deodhar Property Z recursion: for s in D_R(w): x<=w iff min(x,xs) <= ws."""
    while True:
        Dw = right_descents(w)
        if not Dw:
            return x == ID
        s = min(Dw)
        w = matmul(w, S[s])
        if s in right_descents(x):
            x = matmul(x, S[s])

# ---------- FC test via heap convex chains ----------
def is_fc_word(word):
    """word is a reduced word. Returns True iff its commutation class has no factor sts (m=3)."""
    L = len(word)
    # heap order: p<q if p before q and labels don't commute (equal or adjacent); transitive closure
    less = [[False] * L for _ in range(L)]
    for p in range(L):
        for q in range(p + 1, L):
            a, b = word[p], word[q]
            if a == b or b in adj[a]:
                less[p][q] = True
    for k in range(L):
        for i in range(L):
            if less[i][k]:
                for j in range(L):
                    if less[k][j]:
                        less[i][j] = True
    # for consecutive occurrences i<k of same label s, open interval must not be exactly {j} with label adjacent to s
    last = {}
    for k in range(L):
        s = word[k]
        if s in last:
            i = last[s]
            between = [p for p in range(L) if less[i][p] and less[p][k]]
            if len(between) == 1 and word[between[0]] in adj[s]:
                return False
        last[s] = k
    return True

def enumerate_fc(cap=2_000_000):
    """BFS over FC elements via right multiplication. Returns dict matrix->word."""
    seen = {ID: []}
    q = deque([ID])
    while q:
        X = q.popleft()
        wd = seen[X]
        D = right_descents(X)
        for s in range(n):
            if s in D:
                continue
            Y = matmul(X, S[s])
            if Y in seen:
                continue
            if is_fc_word(wd + [s]):
                seen[Y] = wd + [s]
                q.append(Y)
                if len(seen) > cap:
                    raise RuntimeError("cap")
    return seen

if __name__ == "__main__":
    beta = (1, 2, 3, 3, 2, 2, 1, 1, 2)
    print("form(beta,beta) =", form(beta, beta))
    R = reflection_matrix(beta)
    w = reduced_word(R)
    print("length of r_beta =", len(w))
    print("reduced word:", "".join(map(str, w)))
    print("support:", sorted(set(w)), "full support:", set(w) == set(range(n)))
    Rd = right_descents(R); Ld = left_descents(R)
    print("R(b) =", sorted(Rd), " L(b) =", sorted(Ld))
    print("involution (R^2=I):", matmul(R, R) == ID)
    print("right-terminal:", terminal_right(R), " left-terminal:", terminal_right(transpose_inv(R)))
    print("is b FC?", is_fc_word(w))

    # extensions by s0, s1
    exts = {}
    for s in (0, 1):
        E = matmul(R, S[s])
        we = reduced_word(E)
        exts[s] = E
        print(f"b*s{s}: length {len(we)}, R={sorted(right_descents(E))}, L={sorted(left_descents(E))},"
              f" right-terminal={terminal_right(E)}, left-terminal={terminal_right(transpose_inv(E))}")

    # Independence number of E9 graph (brute force)
    import itertools
    best = 0
    for r in range(n, 0, -1):
        found = False
        for sub in itertools.combinations(range(n), r):
            if all(b not in adj[a] for a, b in itertools.combinations(sub, 2)):
                best = r; found = True; break
        if found:
            break
    print("independence number of E9 diagram:", best, " |I4| =", len(Rd))

    # FC catalogue
    print("enumerating FC elements of E9 ...", flush=True)
    fc = enumerate_fc()
    maxlen = max(len(v) for v in fc.values())
    print("number of FC elements:", len(fc), " max length:", maxlen)

    def eligible(bmat):
        Lb = left_descents(bmat); Rb = right_descents(bmat)
        out = []
        for X, wd in fc.items():
            if X == bmat:
                continue
            if not Rb <= right_descents(X):
                continue
            if not Lb <= left_descents(X):
                continue
            if bruhat_leq(X, bmat):
                out.append(wd)
        return sorted(out, key=lambda v: (len(v), v))

    lb = len(w)
    el = eligible(R)
    print(f"eligible FC bottoms for b (len {lb}): {len(el)}")
    for v in el:
        gap = lb - len(v)
        print("  ", "".join(map(str, v)), "len", len(v), "gap", gap, "odd" if gap % 2 else "even")
    for s in (0, 1):
        E = exts[s]; le = len(reduced_word(E))
        el = eligible(E)
        print(f"eligible FC bottoms for b*s{s} (len {le}): {len(el)}")
        for v in el:
            gap = le - len(v)
            print("  ", "".join(map(str, v)), "len", len(v), "gap", gap, "odd" if gap % 2 else "even")
