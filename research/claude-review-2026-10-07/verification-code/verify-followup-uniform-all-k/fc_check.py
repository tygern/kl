"""FC / reducedness checks for words in E_n (chain 0..n-2, node n-1 attached to 2).

is_fc_reduced(word): Stembridge's heap criterion -- between any two consecutive
occurrences of a letter s there must be at least two letters adjacent to s.
This is equivalent to: the word is a reduced word of a fully commutative element.

length_by_matrix(word): independent check of reducedness by computing the element
in the geometric representation (exact integers) and reading off its length by
descent stripping.
"""
import itertools

def adjacency_E(n):
    adj = [set() for _ in range(n)]
    for i in range(n-2):
        adj[i].add(i+1); adj[i+1].add(i)
    adj[2].add(n-1); adj[n-1].add(2)
    return adj

def is_fc_reduced(word, adj):
    last = {}
    between = {}
    for pos, s in enumerate(word):
        if s in last:
            cntN = sum(1 for t in word[last[s]+1:pos] if t in adj[s])
            if cntN < 2:
                return False
        last[s] = pos
    return True

def cartan_E(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n): A[i][i] = 2
    for i in range(n-2): A[i][i+1] = A[i+1][i] = -1
    A[2][n-1] = A[n-1][2] = -1
    return A

def length_by_matrix(word, n):
    """Multiply simple reflections, then strip descents to compute the length."""
    A = cartan_E(n)
    M = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    for s in word:
        # M <- M * s : column j becomes M(alpha_j) - A[s][j] M(alpha_s)
        cs = [M[i][s] for i in range(n)]
        for j in range(n):
            if j == s: continue
            a = A[s][j]
            if a:
                for i in range(n): M[i][j] -= a*cs[i]
        for i in range(n): M[i][s] = -cs[i]
    L = 0
    while True:
        s = -1
        for j in range(n):
            col = [M[i][j] for i in range(n)]
            if all(c <= 0 for c in col) and any(c < 0 for c in col):
                s = j; break
            assert all(c >= 0 for c in col)
        if s < 0:
            assert all(M[i][j] == (1 if i == j else 0) for i in range(n) for j in range(n))
            return L
        cs = [M[i][s] for i in range(n)]
        for j in range(n):
            if j == s: continue
            a = A[s][j]
            if a:
                for i in range(n): M[i][j] -= a*cs[i]
        for i in range(n): M[i][s] = -cs[i]
        L += 1
