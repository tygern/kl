#!/usr/bin/env python3
"""Independent check of Fan (JAMS 1997) / Green (arXiv:0704.0283, Prop 4.4) monomial
cell statements in type E_n, n = 6..10, directly from the defining relations of the
Temperley-Lieb quotient:  b_s^2 = delta b_s,  b_s b_t = b_t b_s (st=ts),  b_s b_t b_s = b_s (m=3).

Labelling (Fan 1997 / Green 2007): nodes 0..n-1; 1-2-...-(n-1) is a path; 0 is joined to 3.

Checks performed for each rank:
  C1  #FC elements matches Fan's table (p.156).
  C2  a(w) := heap width equals the exponent in b_{w^{-1}} b_w = delta^{a(w)} b_d (Fan Lemma 5.2.6),
      and d is an involution.
  C3  a is constant on left, right and two-sided (monomial) cells (Green 4.4(ii)).
  C4  Prop 4.4.3: right cell of i(S) == {w : S subset of L(w), a(w)=|S|} for every independent S.
  C5  dual of 4.4.3 (left cells).
  C6  Every left cell and right cell in the same two-sided cell meet in exactly one element
      (Green 4.4(v), Fan p.140 / proof of Thm 6.1.2).
  C7  Two-sided cells are indexed by P/~ (Fan Thm 4.5.1): number of two-sided cells and the
      left-cell sizes match Fan's dimension table (p.156); #T = (#left cells in T)^2.
  C8  Green 4.4(iii): w' <=_R w, w' not ~_R w  implies a(w') > a(w)  (checked on single-step edges).
  C9  The exact Lemma maxI instance: for every FC x, if L(x) cap R(x) contains an independent set of
      size alpha(supp x) (independence number of induced subgraph) then x = i(L(x) cap R(x)).
"""
import sys, time
from collections import defaultdict, deque

sys.setrecursionlimit(100000)

def en_graph(n):
    adj = {i: set() for i in range(n)}
    for i in range(1, n - 1):
        adj[i].add(i + 1); adj[i + 1].add(i)
    adj[0].add(3); adj[3].add(0)
    return adj

class TL:
    def __init__(self, n):
        self.n = n
        self.adj = en_graph(n)
        self.comm = [[(i != j and j not in self.adj[i]) for j in range(n)] for i in range(n)]

    # ---------- words / canonical forms ----------
    def canon(self, word):
        """Cartier-Foata normal form of a reduced FC word (unique in commutation class)."""
        comm = self.comm
        w = list(word); out = []
        while w:
            layer = []
            rest = []
            blocked = set()  # letters that have a non-commuting predecessor already seen
            seen = []
            for s in w:
                if s in blocked:
                    rest.append(s)
                else:
                    layer.append(s)
                # everything not commuting with s that comes later is blocked
                for t in range(self.n):
                    if not comm[s][t]:
                        blocked.add(t)
            layer.sort()
            out.extend(layer)
            w = rest
        return tuple(out)

    def is_fc_reduced_append(self, word, s):
        """word reduced FC; is word+(s,) reduced FC? (property R3 on the new letter)"""
        comm = self.comm
        cnt = 0
        for t in reversed(word):
            if t == s:
                return cnt >= 2
            if not comm[s][t]:
                cnt += 1
        return True

    def left_desc(self, word):
        comm = self.comm
        L = set(); blocked = set()
        for s in word:
            if s not in blocked:
                L.add(s)
            for t in range(self.n):
                if not comm[s][t]:
                    blocked.add(t)
        return frozenset(L)

    def right_desc(self, word):
        return self.left_desc(tuple(reversed(word)))

    # ---------- heap width = a-function (Fan Def 2.3.1 / remark p.142) ----------
    def heap_width(self, word):
        k = len(word)
        comm = self.comm
        # comparability: i<j comparable iff some chain of non-commuting letters; compute transitive closure
        below = [0] * k  # bitmask of elements below i
        for i in range(k):
            m = 0
            for j in range(i):
                if not comm[word[i]][word[j]]:
                    m |= (1 << j) | below[j]
            below[i] = m
        # width via Dilworth: k - max matching in bipartite graph (i -> j if i < j comparable)
        succ = [[j for j in range(i + 1, k) if (below[j] >> i) & 1] for i in range(k)]
        match = [-1] * k
        def try_aug(i, seen):
            for j in succ[i]:
                if j in seen: continue
                seen.add(j)
                if match[j] == -1 or try_aug(match[j], seen):
                    match[j] = i
                    return True
            return False
        mm = 0
        for i in range(k):
            if try_aug(i, set()):
                mm += 1
        return k - mm

    # ---------- monomial multiplication ----------
    def lmul(self, s, w):
        """b_s b_w = delta^m b_{w'}; w any reduced FC word; returns (m, w') with w' a reduced FC word."""
        comm = self.comm
        idx = -1
        for i, t in enumerate(w):
            if t == s:
                idx = i; break
        if idx == -1:
            return 0, (s,) + tuple(w)
        nc = [j for j in range(idx) if not comm[s][w[j]]]
        if not nc:
            return 1, tuple(w)
        if len(nc) >= 2:
            return 0, (s,) + tuple(w)
        j = nc[0]
        w1 = w[:j]; w2 = w[j + 1:idx] + w[idx + 1:]
        m, u = self.lmul(s, w2)
        for t in reversed(w1):
            m2, u = self.lmul(t, u)
            m += m2
        return m, u

    def rmul(self, w, s):
        m, u = self.lmul(s, tuple(reversed(w)))
        return m, tuple(reversed(u))

def scc(nodes, out_edges):
    """iterative Tarjan; returns comp id per node"""
    index = {}; low = {}; onstack = set(); stack = []; comp = {}; cid = 0; idx = 0
    for root in nodes:
        if root in index: continue
        work = [(root, iter(out_edges[root]))]
        index[root] = low[root] = idx; idx += 1; stack.append(root); onstack.add(root)
        while work:
            v, it = work[-1]
            advanced = False
            for w in it:
                if w not in index:
                    index[w] = low[w] = idx; idx += 1; stack.append(w); onstack.add(w)
                    work.append((w, iter(out_edges[w]))); advanced = True; break
                elif w in onstack:
                    low[v] = min(low[v], index[w])
            if advanced: continue
            work.pop()
            if work:
                u = work[-1][0]; low[u] = min(low[u], low[v])
            if low[v] == index[v]:
                while True:
                    x = stack.pop(); onstack.discard(x); comp[x] = cid
                    if x == v: break
                cid += 1
    return comp, cid

def independent_sets(adj, nodes):
    nodes = sorted(nodes)
    res = []
    def rec(i, cur):
        if i == len(nodes):
            res.append(frozenset(cur)); return
        rec(i + 1, cur)
        v = nodes[i]
        if all(u not in adj[v] for u in cur):
            cur.append(v); rec(i + 1, cur); cur.pop()
    rec(0, [])
    return res

def indep_number(adj, nodes):
    return max(len(S) for S in independent_sets(adj, nodes))

FAN_TABLE = {
    11: (737762, [1, 11, 65, 220, 506, 527, 341, 187]),  # rank: (#Wc, sorted left-cell sizes per two-sided cell = irreducible dims)
    6: (662, [1, 6, 20, 15]),
    7: (2670, [1, 7, 27, 35, 21, 15]),
    8: (10846, [1, 8, 35, 84, 50]),
    9: (44199, [1, 9, 44, 120, 135, 84, 50]),
    10: (180438, [1, 10, 54, 165, 340, 186]),
}

def run(n):
    t0 = time.time()
    A = TL(n)
    # enumerate FC elements by right-extension BFS on canonical words
    start = ()
    elems = {start}
    q = deque([start])
    while q:
        w = q.popleft()
        for s in range(n):
            if A.is_fc_reduced_append(w, s):
                u = A.canon(w + (s,))
                if u not in elems:
                    elems.add(u); q.append(u)
    elems = sorted(elems, key=lambda w: (len(w), w))
    N = len(elems)
    report = {}
    report['C1_count'] = (N, FAN_TABLE[n][0], N == FAN_TABLE[n][0])
    print(f"E{n}: {N} FC elements (Fan table {FAN_TABLE[n][0]}) [{time.time()-t0:.1f}s]", flush=True)

    aval = {w: A.heap_width(w) for w in elems}
    Ld = {w: A.left_desc(w) for w in elems}
    Rd = {w: A.right_desc(w) for w in elems}

    # multiplication tables
    lout = defaultdict(set); rout = defaultdict(set)
    lstep = {}; rstep = {}
    for w in elems:
        for s in range(n):
            m, u = A.lmul(s, w); u = A.canon(u)
            assert u in aval, (s, w, u)
            lstep[(s, w)] = (m, u); lout[w].add(u)
            m2, u2 = A.rmul(w, s); u2 = A.canon(u2)
            assert u2 in aval
            rstep[(w, s)] = (m2, u2); rout[w].add(u2)
    print(f"  multiplication tables done [{time.time()-t0:.1f}s]", flush=True)

    # C2: b_{w^-1} b_w = delta^{a(w)} b_d, d involution
    bad2 = 0
    for w in elems:
        m = 0; u = w
        for t in w:  # b_{w^{-1}} = b_{s_k}...b_{s_1}; apply s_1 first
            m2, u = lstep[(t, u)]; m += m2
        inv_u = A.canon(tuple(reversed(u)))
        if m != aval[w] or inv_u != u:
            bad2 += 1
    report['C2_lemma526'] = bad2
    print(f"  C2 Lemma 5.2.6 failures: {bad2}", flush=True)

    # cells
    lcomp, nl = scc(elems, lout)
    rcomp, nr = scc(elems, rout)
    both = {w: lout[w] | rout[w] for w in elems}
    tcomp, nt = scc(elems, both)
    print(f"  left cells {nl}, right cells {nr}, two-sided cells {nt} [{time.time()-t0:.1f}s]", flush=True)

    # C3 a constant on cells
    def const_on(comp):
        vals = defaultdict(set)
        for w in elems: vals[comp[w]].add(aval[w])
        return all(len(v) == 1 for v in vals.values())
    report['C3_a_constant'] = (const_on(lcomp), const_on(rcomp), const_on(tcomp))
    print(f"  C3 a constant on (L,R,LR) cells: {report['C3_a_constant']}", flush=True)

    # C4/C5 Prop 4.4.3 and dual
    P = independent_sets(A.adj, range(n))
    bad4 = []; bad5 = []
    rcell_members = defaultdict(set); lcell_members = defaultdict(set)
    for w in elems:
        rcell_members[rcomp[w]].add(w); lcell_members[lcomp[w]].add(w)
    for S in P:
        iS = A.canon(tuple(sorted(S)))
        pred_R = {w for w in elems if S <= Ld[w] and aval[w] == len(S)}
        if rcell_members[rcomp[iS]] != pred_R: bad4.append(S)
        pred_L = {w for w in elems if S <= Rd[w] and aval[w] == len(S)}
        if lcell_members[lcomp[iS]] != pred_L: bad5.append(S)
    report['C4_prop443'] = (len(P), bad4)
    report['C5_dual443'] = (len(P), bad5)
    print(f"  C4 Prop 4.4.3 checked for {len(P)} independent sets, failures {bad4}", flush=True)
    print(f"  C5 dual 4.4.3 failures {bad5}", flush=True)

    # C6 singleton intersections
    tcells = defaultdict(list)
    for w in elems: tcells[tcomp[w]].append(w)
    bad6 = 0; dims = []
    for T, ws in tcells.items():
        ls = {lcomp[w] for w in ws}; rs = {rcomp[w] for w in ws}
        cnt = defaultdict(int)
        for w in ws: cnt[(lcomp[w], rcomp[w])] += 1
        if len(cnt) != len(ls) * len(rs) or any(c != 1 for c in cnt.values()):
            bad6 += 1
        if len(ls) != len(rs) or len(ws) != len(ls) ** 2:
            bad6 += 1
        dims.append(len(ls))
    report['C6_singleton'] = bad6
    report['C7_dims'] = (sorted(dims), sorted(FAN_TABLE[n][1]), sorted(dims) == sorted(FAN_TABLE[n][1]), nt)
    # number of ~ classes of P (Fan 2.3): neighbours relation
    Pset = set(P)
    pid = {S: i for i, S in enumerate(P)}
    parent = list(range(len(P)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for S in P:
        for s in S:
            for t in A.adj[s]:
                T = (S - {s}) | {t}
                if T in Pset:
                    a_, b_ = find(pid[S]), find(pid[T])
                    if a_ != b_: parent[a_] = b_
    nPbar = len({find(i) for i in range(len(P))})
    report['C7_Pbar'] = (nPbar, nt, nPbar == nt)
    print(f"  C6 singleton-intersection failures: {bad6}; dims {sorted(dims)} vs Fan {sorted(FAN_TABLE[n][1])}; |P/~|={nPbar}", flush=True)

    # C8 Green 4.4(iii) on single steps
    bad8 = 0
    for (w, s), (m, u) in rstep.items():
        if rcomp[u] != rcomp[w] and not (aval[u] > aval[w]): bad8 += 1
    for (s, w), (m, u) in lstep.items():
        if lcomp[u] != lcomp[w] and not (aval[u] > aval[w]): bad8 += 1
    report['C8_443iii'] = bad8
    print(f"  C8 4.4(iii) single-step failures: {bad8}", flush=True)

    # C9 Lemma maxI exact instance
    bad9 = []
    alpha_cache = {}
    hits = 0
    for w in elems:
        supp = frozenset(w)
        if supp not in alpha_cache:
            alpha_cache[supp] = indep_number(A.adj, supp) if supp else 0
        al = alpha_cache[supp]
        common = Ld[w] & Rd[w]
        if len(common) >= al:  # common descents are automatically independent (R6)
            hits += 1
            if not (len(w) == len(common) == al):
                bad9.append(w)
    report['C9_maxI'] = (hits, bad9)
    print(f"  C9 Lemma maxI: {hits} elements with |L∩R| >= alpha(supp); violations {len(bad9)} {bad9[:5]}", flush=True)
    print(f"E{n} done in {time.time()-t0:.1f}s", flush=True)
    return report

if __name__ == '__main__':
    ranks = [int(a) for a in sys.argv[1:]] or [6, 7, 8, 9]
    import json
    allrep = {}
    for n in ranks:
        allrep[n] = run(n)
    with open(f"report_{'_'.join(map(str, ranks))}.json", 'w') as f:
        json.dump({str(k): {kk: (list(vv) if isinstance(vv, (set, frozenset)) else repr(vv)) for kk, vv in v.items()} for k, v in allrep.items()}, f, indent=1)
