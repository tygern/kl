"""Part (b)/(c): structural invariants of the two new exceptional terminals."""
from fractions import Fraction
from itertools import combinations
from collections import deque
import sys, json
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-new-terminals-structure')
from en_tools import *

E7W = '1325436210321432543621324356'
E8W = '7534231270123456210321432' + '5437210321432543721324357'
E7_IN_E8 = E7W.replace('6', '7')
GERN_W6_E7 = [gern_to_en(s, 7) for s in gern_wn_word(6)]
GERN_W6_E8 = [gern_to_en(s, 8) for s in gern_wn_word(6)]
GERN_W4_E7 = [gern_to_en(s, 7) for s in gern_wn_word(4)]
GERN_W4_E8 = [gern_to_en(s, 8) for s in gern_wn_word(4)]


def charpoly(M, n):
    """Faddeev-LeVerrier, exact. M as columns M[j][i]."""
    A = [[Fraction(M[j][i]) for j in range(n)] for i in range(n)]
    I = [[Fraction(int(i == j)) for j in range(n)] for i in range(n)]
    def mm(X, Y):
        return [[sum(X[i][k] * Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    coeffs = [Fraction(1)]
    Mk = [[Fraction(0)] * n for _ in range(n)]
    for k in range(1, n + 1):
        Mk = mm(A, Mk)
        for i in range(n):
            Mk[i][i] += coeffs[-1]
        AM = mm(A, Mk)
        c = -sum(AM[i][i] for i in range(n)) / k
        coeffs.append(c)
    return [int(c) for c in coeffs]  # x^n + c1 x^{n-1} + ...


def rank_of(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    ncol = len(rows[0]) if rows else 0
    for c in range(ncol):
        p = next((r for r in range(rk, len(rows)) if rows[r][c] != 0), None)
        if p is None:
            continue
        rows[rk], rows[p] = rows[p], rows[rk]
        piv = rows[rk][c]
        rows[rk] = [x / piv for x in rows[rk]]
        for r in range(len(rows)):
            if r != rk and rows[r][c] != 0:
                f = rows[r][c]
                rows[r] = [x - f * y for x, y in zip(rows[r], rows[rk])]
        rk += 1
    return rk


def matrix_order(rs, M):
    I = rs.identity()
    P = M
    for k in range(1, 100000):
        if P == I:
            return k
        P = rs.mult(P, M)
    return None


def components(rs, roots):
    """Connected components of the Dynkin graph of a set of roots (edges when pairing != 0)."""
    comps = []
    seen = set()
    for i in range(len(roots)):
        if i in seen:
            continue
        comp = [i]
        seen.add(i)
        q = [i]
        while q:
            a = q.pop()
            for b in range(len(roots)):
                if b not in seen and rs.pairing(roots[a], roots[b]) != 0:
                    seen.add(b)
                    comp.append(b)
                    q.append(b)
        comps.append(comp)
    return comps


def subsystem_type(rs, roots):
    """Name the root subsystem spanned by the given positive roots (closure in Phi)."""
    pos = rs.positive_roots()
    span_rank = rank_of([list(r) for r in roots])
    # closure: all positive roots in the Q-span
    def in_span(v):
        return rank_of([list(r) for r in roots] + [list(v)]) == span_rank
    sub = [r for r in pos if in_span(r)]
    # simple roots of sub: those not a sum of two others in sub
    subset = set(sub)
    simple = [r for r in sub if not any(tuple(a - b for a, b in zip(r, s)) in subset for s in sub if s != r)]
    comps = components(rs, simple)
    names = []
    for comp in comps:
        k = len(comp)
        npos = sum(1 for r in sub if rank_of([list(simple[i]) for i in comp] + [list(r)]) == k)
        # identify by rank and number of positive roots
        table = {(1, 1): 'A1', (2, 3): 'A2', (3, 6): 'A3', (4, 10): 'A4', (5, 15): 'A5', (6, 21): 'A6', (7, 28): 'A7',
                 (4, 12): 'D4', (5, 20): 'D5', (6, 30): 'D6', (7, 42): 'D7', (8, 56): 'D8',
                 (6, 36): 'E6', (7, 63): 'E7', (8, 120): 'E8'}
        names.append(table.get((k, npos), f'rank{k}/{npos}'))
    return sorted(names), len(sub), simple


def analyze(n, wstr, label):
    rs = RootSystem(n, en_edges(n))
    word = parse(wstr)
    M = rs.word_matrix(word)
    Minv = rs.inverse(M)
    ell = rs.length(M)
    print(f"\n######## {label}: E_{n}, word {wstr}")
    print(f"length (inversions) = {ell}, word length = {len(word)}, reduced: {ell == len(word)}")
    L, R = rs.left_descents(M), rs.right_descents(M)
    print(f"L = {sorted(L)}, R = {sorted(R)}, terminal = {rs.terminal(M)}, support = {sorted(rs.support(M))}")
    print(f"involution (w^2 = 1): {rs.mult(M, M) == rs.identity()};  w == w^{{-1}}: {M == Minv}")
    order = matrix_order(rs, M)
    print(f"order of w: {order}")
    cp = charpoly(M, n)
    print(f"characteristic polynomial coefficients (x^n + c1 x^(n-1) + ...): {cp}")
    MI = [[M[j][i] + int(i == j) for j in range(n)] for i in range(n)]
    MmI = [[M[j][i] - int(i == j) for j in range(n)] for i in range(n)]
    dim_minus = n - rank_of(MI)
    dim_plus = n - rank_of(MmI)
    print(f"dim (-1)-eigenspace = {dim_minus}; dim (+1)-eigenspace = {dim_plus}")
    # roots negated / fixed
    pos = rs.positive_roots()
    negated = [r for r in pos if rs.apply(M, r) == tuple(-c for c in r)]
    fixed = [r for r in pos if rs.apply(M, r) == r]
    print(f"positive roots negated by w: {len(negated)}; fixed: {len(fixed)}")
    if negated:
        names, npos, simple = subsystem_type(rs, negated)
        print(f"  root subsystem Phi_- spanned by negated roots: type {names}, {npos} positive roots, rank {rank_of([list(r) for r in negated])}")
        print(f"  simple roots of Phi_-: {[list(r) for r in simple]}")
        # orthogonal decomposition into reflections if involution
        if rs.mult(M, M) == rs.identity():
            # greedy maximal orthogonal set within negated roots, check product
            best = None
            for start in negated:
                chosen = [start]
                for r in negated:
                    if all(rs.pairing(r, c) == 0 for c in chosen):
                        chosen.append(r)
                if len(chosen) == dim_minus:
                    best = chosen
                    break
            if best:
                P = rs.identity()
                for r in best:
                    # reflection matrix r_beta
                    Rm = tuple(tuple(a - rs.pairing(tuple(int(i == j) for i in range(n)), r) * b
                                     for a, b in zip(tuple(int(i == j) for i in range(n)), r)) for j in range(n))
                    P = rs.mult(P, Rm)
                print(f"  w == product of {len(best)} mutually orthogonal reflections {[list(r) for r in best]}: {P == M}")
                print(f"  heights of those roots: {[sum(r) for r in best]}; reflection lengths 2ht-1: {[2*sum(r)-1 for r in best]}")
    if fixed:
        names, npos, simple = subsystem_type(rs, fixed)
        print(f"  root subsystem Phi_+ of fixed roots: type {names}, {npos} positive roots")
    inv = rs.inversions(M)
    print(f"inversion set N(w) (positive roots sent negative), by height:")
    byh = {}
    for r in inv:
        byh.setdefault(sum(r), []).append(''.join(str(c) for c in r))
    for h in sorted(byh):
        print(f"  ht {h:2d}: {' '.join(byh[h])}")
    # positive roots not inverted (complement) heights
    print(f"  max inverted height {max(sum(r) for r in inv)}; highest root height {max(sum(r) for r in pos)}")

    # ---- double coset maximality ----
    print("double coset check: is w the longest element of W_I w W_J for I<=L, J<=R?")
    found = []
    for i in range(len(L) + 1):
        for I in combinations(sorted(L), i):
            for j in range(len(R) + 1):
                for J in combinations(sorted(R), j):
                    if not I and not J:
                        continue
                    # enumerate W_I w W_J
                    seen = {M}
                    q = deque([M])
                    while q:
                        X = q.popleft()
                        for s in I:
                            Y = rs.left_mult(s, X)
                            if Y not in seen:
                                seen.add(Y); q.append(Y)
                        for s in J:
                            Y = rs.right_mult(X, s)
                            if Y not in seen:
                                seen.add(Y); q.append(Y)
                    maxlen = max(rs.length(X) for X in seen)
                    if maxlen == ell:
                        found.append((I, J, len(seen)))
    print(f"  w is longest in its double coset for (I,J) = {found if found else 'NONE'}")
    # Since I<=L(w), J<=R(w) always, also report the full one:
    I, J = tuple(sorted(L)), tuple(sorted(R))
    seen = {M}; q = deque([M])
    while q:
        X = q.popleft()
        for s in I:
            Y = rs.left_mult(s, X)
            if Y not in seen: seen.add(Y); q.append(Y)
        for s in J:
            Y = rs.right_mult(X, s)
            if Y not in seen: seen.add(Y); q.append(Y)
    lens = sorted(rs.length(X) for X in seen)
    print(f"  |W_L w W_R| = {len(seen)}, lengths range {lens[0]}..{lens[-1]}, w has {ell}")

    # ---- weak-order suffix BFS: factorizations w = u * v reduced ----
    # suffixes v with w = u v, l(w)=l(u)+l(v): v obtained by stripping left descents of w.
    suffixes = {M: 0}
    q = deque([M])
    while q:
        V = q.popleft()
        for s in rs.left_descents(V):
            Y = rs.left_mult(s, V)
            if Y not in suffixes:
                suffixes[Y] = suffixes[V] + 1
                q.append(Y)
    print(f"right weak order interval [e,w]: {len(suffixes)} elements (= number of distinct suffixes)")
    # maximal parabolic longest-element factors: for suffix v, J = L(v) gives factor w_0(J)
    best = {}
    for V in suffixes:
        J = tuple(sorted(rs.left_descents(V)))
        if not J:
            continue
        w0 = rs.longest_element(list(J))
        lw0 = rs.length(w0)
        key = J
        if key not in best or best[key] < lw0:
            best[key] = lw0
    top = sorted(best.items(), key=lambda kv: -kv[1])[:8]
    print(f"  longest-element factors w = u w_0(J) v reduced, top J by l(w_0(J)): {top}")
    maxJ, maxl = top[0]
    print(f"  => a(w) >= {maxl} via w_0({maxJ}); parity-free vanishing needs a(w) >= |I|+2 = {len(L)+2}: {maxl >= len(L)+2}")
    # components of top J
    rsJ_edges = [(a, b) for a, b in rs.edges if a in maxJ and b in maxJ]
    print(f"  edges inside that J: {rsJ_edges}")

    # ---- inner-factor containment tests ----
    def is_prefix(F, V):
        """F is a prefix of V (V = F V' reduced) iff N(F^{-1}) subset N(V^{-1})."""
        Fi, Vi = rs.inverse(F), rs.inverse(V)
        NF = set(rs.inversions(Fi)); NV = set(rs.inversions(Vi))
        return NF <= NV
    def contains_factor(Fword, name):
        F = rs.word_matrix(Fword)
        hits = [V for V in suffixes if is_prefix(F, V)]
        print(f"  factor {name} (len {rs.length(F)}) occurs as inner length-additive factor: {bool(hits)}"
              + (f" (e.g. after prefix of length {min(suffixes[V] for V in hits)})" if hits else ''))
        return hits
    return rs, M, suffixes, contains_factor


rs7, M7, suf7, cf7 = analyze(7, E7W, 'E7 new terminal')
print("inner factors in E7 element:")
cf7(GERN_W4_E7, 'Gern w4 (E7 labels)')
cf7(GERN_W6_E7, 'Gern w6 (E7 labels)')
cf7(parse('1326213'), 'table row 1326213')
cf7(parse('132543621324356'), 'table row 132543621324356')

rs8, M8, suf8, cf8 = analyze(8, E8W, 'E8 new terminal b_{8,6}')
print("inner factors in E8 element:")
cf8(GERN_W4_E8, 'Gern w4 (E8 labels)')
cf8(GERN_W6_E8, 'Gern w6 (E8 labels)')
cf8(parse(E7_IN_E8), 'E7 new element relabelled 6->7')
cf8(parse(E7_IN_E8[3:]), 'E7 new element minus prefix 132')
cf8(parse(E7W), 'E7 word with label 6 (A7-type subdiagram) literal')

# ---- parabolic decompositions w = a v, a in W^J, v in W_J ----
def parabolic_decomp(rs, M, J):
    """Strip right descents in J: v accumulates; a = M v^{-1}."""
    V = rs.identity(); A = M
    vword = []
    while True:
        R = rs.right_descents(A) & set(J)
        if not R:
            break
        s = min(R)
        A = rs.right_mult(A, s)
        vword.append(s)
    vword = vword[::-1]
    V = rs.word_matrix(vword)
    assert rs.mult(A, V) == M
    return A, V, vword

print("\n######## parabolic decompositions w = a v (a minimal coset rep, v in W_J)")
for (n, rs, M, label, Js) in [
        (7, rs7, M7, 'E7 new', {'D6={1..6}': [1,2,3,4,5,6], 'E6={0,1,2,3,4,6}': [0,1,2,3,4,6], 'A6={0..5}': [0,1,2,3,4,5]}),
        (8, rs8, M8, 'E8 new', {'E7={0..5,7}': [0,1,2,3,4,5,7], 'D7={1..7}': [1,2,3,4,5,6,7], 'A7={0..6}': [0,1,2,3,4,5,6], 'E6={0,1,2,3,4,7}':[0,1,2,3,4,7]})]:
    for name, J in Js.items():
        A, V, vword = parabolic_decomp(rs, M, J)
        print(f"{label} wrt {name}: l(a)={rs.length(A)}, l(v)={rs.length(V)}, v word = {word_str(vword)}, "
              f"v terminal in W_J: {rs.terminal(V)}, a word = {word_str(rs.reduced_word(A))}")
        # and the left version w = v' a'
        Ai, Vi, vw = parabolic_decomp(rs, rs.inverse(M), J)
        print(f"   left version w = v' a': l(v')={rs.length(Vi)}, v' word = {word_str(vw[::-1])}, v' terminal: {rs.terminal(Vi)}")

# ---- compare with Gern's w8 in D8 (non-parabolic subsystem of E8): eigen-data ----
print("\n######## Gern's w_n as signed permutations: (-1)-eigenspace dimensions")
for n in (4, 6, 8, 10):
    sp = signed_perm_of_gern_word(gern_wn_word(n), n)
    # matrix in e-basis
    Mat = [[0]*n for _ in range(n)]
    for i, j in enumerate(sp, 1):
        Mat[abs(j)-1][i-1] = 1 if j > 0 else -1
    MI = [[Mat[i][j] + int(i == j) for j in range(n)] for i in range(n)]
    print(f"w_{n} = {sp}: dim(-1)-eigenspace = {n - rank_of(MI)}, ell = {len(gern_wn_word(n))}")
