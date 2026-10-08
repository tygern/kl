"""Test of the conjectured lemma on small trees (FC-finite): for a tree Coxeter graph G and a
maximum independent set I, every FC x with L(x) ⊇ I and |R(x) ∩ I| ≥ |I|-1 has
l(x) ≤ |I|+1 (and is i(I) or i(I)t for a single generator t).  Height-vector representation
(h_i = height of w(alpha_i)); FC closure by levels with the descent recurrence."""
import sys, itertools

def fc_elements(adj):
    n = len(adj)
    e = tuple([1] * n)
    def rmul(h, j):
        h = list(h)
        for i in adj[j]:
            h[i] += h[j]
        h[j] = -h[j]
        return tuple(h)
    levels = [{e}]
    while True:
        prev = levels[-1]
        nxt = set()
        for u in prev:
            for s in range(n):
                if u[s] < 0:
                    continue
                v = rmul(u, s)
                if v in nxt:
                    continue
                R = [j for j in range(n) if v[j] < 0]
                ok = all(b not in adj[a] for a in R for b in R if a < b)
                if not ok:
                    continue
                for t in R:
                    if t != s and rmul(v, t) not in prev:
                        ok = False
                        break
                if ok:
                    nxt.add(v)
        if not nxt:
            break
        levels.append(nxt)
    return levels, rmul

def word_of(h, adj, rmul):
    n = len(adj)
    w = []
    while True:
        d = next((j for j in range(n) if h[j] < 0), None)
        if d is None:
            break
        w.append(d); h = rmul(h, d)
    return w[::-1]

def max_independent_sets(adj):
    n = len(adj)
    best = 0; sets = []
    for mask in range(1 << n):
        S = [i for i in range(n) if mask >> i & 1]
        if all(b not in adj[a] for a in S for b in S):
            if len(S) > best:
                best = len(S); sets = [S]
            elif len(S) == best:
                sets.append(S)
    return sets

def test(name, adj):
    n = len(adj)
    levels, rmul = fc_elements(adj)
    total = sum(len(l) for l in levels)
    e = tuple([1] * n)
    results = []
    for I in max_independent_sets(adj):
        Iset = set(I)
        worst = []
        for lv in levels:
            for h in lv:
                R = {j for j in range(n) if h[j] < 0}
                if len(R & Iset) < len(I) - 1:
                    continue
                w = word_of(h, adj, rmul)
                inv = e
                for t in reversed(w):
                    inv = rmul(inv, t)
                L = {j for j in range(n) if inv[j] < 0}
                if len(L & Iset) < len(I) - 1:
                    continue
                if not (Iset <= L or Iset <= R):
                    continue
                worst.append((len(w), w, sorted(L), sorted(R)))
        mx = max(x[0] for x in worst)
        bad = [x for x in worst if x[0] > len(I) + 1]
        results.append((I, mx, len(worst), bad[:3]))
    print(f"{name}: n={n} FC={total} maxlen={len(levels)-1}")
    for I, mx, cnt, bad in results:
        flag = "" if not bad else "  VIOLATION"
        print(f"   I={I} |I|={len(I)}: {cnt} elements, max length {mx} (bound {len(I)+1}){flag}")
        for b in bad:
            print("      ", b)

def path(n):
    return [[j for j in (i - 1, i + 1) if 0 <= j < n] for i in range(n)]

def tree_T(a, b, c):
    """Tree with centre 0 and three arms of lengths a,b,c."""
    adj = [[]]
    for arm in (a, b, c):
        prev = 0
        for _ in range(arm):
            k = len(adj); adj.append([prev]); adj[prev].append(k); prev = k
    return adj

if __name__ == '__main__':
    for n in range(2, 9):
        test(f"A{n}", path(n))
    for n in range(4, 9):
        test(f"D{n}", tree_T(1, 1, n - 3))
    for n in range(6, 10):
        test(f"E{n}", tree_T(1, 2, n - 4))
    test("B3-like star K1,3 plus? (T(1,1,1)=D4)", tree_T(1, 1, 1))
    test("T(1,1,5)=D8", tree_T(1, 1, 5))
    test("T(2,2,2)=~E6 (FC-infinite: skip)", path(1)) if False else None
    test("T(1,2,2)=E6 again", tree_T(1, 2, 2))
    # a tree with two branch nodes (FC-finite? H-shaped tree: check quickly up to a cap is not possible; skip)
