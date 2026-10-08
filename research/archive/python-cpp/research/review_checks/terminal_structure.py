"""Structure of the noncommuting terminal elements of E6, E7 and E8.

Certifies, with exact integer matrices and no imports from the repository's
enumeration engines, the statements of the manuscript's remarks on the
finite exceptional terminals (Section 3):

1. Gern's list plus two.  Under the label map 1 -> 1, 2 -> n-1, j -> j-1
   (j = 3, ..., n-1) from Gern's D_{n-1} labels into the paper's E_n labels,
   Gern's bad elements w_4, w_4 u (u a commuting generator neither in nor
   adjacent to supp w_4) and w_6 are exactly the rows of Table 1 of lengths
   7, 8 and 15, as matrices; the E8 row of length 28 is the E7 element w_7
   inside the parabolic on {0,...,5,7}; the only rows not supported on a
   proper type-D or type-E parabolic are w_7 (E7, length 28) and w_8 (E8,
   length 50).  Gern's Theorem 2.3.6 is re-derived for D5, D6 (and D7 with
   --full) by exhaustive enumeration of signed permutations.
2. Each of w_4, w_6, w_7, w_8 is an involution with L = R = I commuting, a
   reduced palindrome of independent-set layers alternating between the two
   colour classes of the bipartite Dynkin diagram, and a product of mutually
   orthogonal reflections (1, 3, 4, 2 of them) in the listed roots; the fixed
   root subsystem of w_8 has rank 6 with 30 positive roots (type D6).
3. The chain w_4 < w_6 < w_7 < w_8 in right weak order (prefix order), inside
   E8, and w_4 < w_6 < w_7 inside E7.
4. Length-additive factorizations b = u w_0(J) v with J = I cup {0} of type
   A2 x A1^3, l(w_0(J)) = 6, for b = w_7 and w_8, and the maximum length of a
   longest parabolic element occurring as a length-additive inner factor:
   3 for w_4, 5 for w_6, 6 for w_7 and w_8.

Conventions.  E_n: chain 0-1-...-(n-2), node n-1 attached to node 2.  A
string of labels denotes the product of the simple reflections in the order
written; W acts on the left and the matrix of w has columns w(alpha_j) in the
simple-root basis.  Gern's D_m: nodes 1 and 2 attached to 3, chain 3-...-m;
s_1 = (1,-2)(-1,2) and s_i = (i-1, i) as signed permutations.

Run: python3 research/review_checks/terminal_structure.py
Writes results/terminal-structure-certificate.json (repository-relative).
Python 3.10+, standard library only.
"""
from collections import deque
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]

TABLE = [  # (n, word, common descent set) as printed in Table 1
    (6, '1325213', '135'),
    (7, '1326213', '136'),
    (7, '13256213', '1356'),
    (7, '132543621324356', '1356'),
    (7, '1325436210321432543621324356', '1356'),
    (8, '1327213', '137'),
    (8, '13257213', '1357'),
    (8, '61327213', '1367'),
    (8, '132543721324357', '1357'),
    (8, '1325437210321432543721324357', '1357'),
    (8, '7534231270123456210321432' + '5437210321432543721324357', '1357'),
]
LAYERS = {  # layered palindromic words of the four maximal noncommuting terminals
    'w4': (6, '1325213', [[1, 3], [2], [5], [2], [1, 3]]),
    'w6': (7, '132543621324356', [[1, 3, 5, 6], [2, 4], [1, 3, 6], [2, 4], [1, 3, 5, 6]]),
    'w7': (7, '1325436210321432543621324356',
           [[1, 3, 5, 6], [2, 4], [1, 3, 6], [0, 2, 4], [1, 3, 5, 6], [0, 2, 4], [1, 3, 6], [2, 4], [1, 3, 5, 6]]),
    'w8': (8, '7534231270123456210321432' + '5437210321432543721324357',
           [[1, 3, 5, 7], [2, 4], [1, 3, 7], [0, 2, 4], [1, 3, 5, 7], [0, 2, 4, 6], [1, 3, 5, 7], [2, 4],
            [1, 3, 5, 7], [0, 2, 4, 6], [1, 3, 5, 7], [0, 2, 4], [1, 3, 7], [2, 4], [1, 3, 5, 7]]),
}
EXPECTED_NEGATED = {
    'w4': ['011101'], 'w6': ['0011001', '0111000', '0111111'],
    'w7': ['0011111', '0111101', '0121001', '1222111'], 'w8': ['12332212', '13432102'],
}
EXPECTED_FACTORIZATION = {
    'w7': ('31265234312', '653010', '42312645231'),
    'w8': ('31275234312', '753010', '276453423127563452341230127345231'),
}


def parse(word):
    return [int(c) for c in word]


def word_str(word):
    return ''.join(map(str, word))


class RootSystem:
    """Simply laced root system given by Dynkin edges; elements as column tuples."""

    def __init__(self, n, edges):
        self.n, self.edges = n, edges
        self.adj = {i: set() for i in range(n)}
        for a, b in edges:
            self.adj[a].add(b)
            self.adj[b].add(a)
        self.identity = tuple(tuple(int(i == j) for i in range(n)) for j in range(n))
        self._roots = None

    def pairing(self, u, v):
        return sum(2 * u[i] * v[i] - u[i] * sum(v[t] for t in self.adj[i]) for i in range(self.n))

    def simple_reflect(self, s, v):
        c = 2 * v[s] - sum(v[t] for t in self.adj[s])
        w = list(v)
        w[s] -= c
        return tuple(w)

    def positive_roots(self):
        if self._roots is None:
            simple = [tuple(int(i == j) for i in range(self.n)) for j in range(self.n)]
            seen, frontier = set(simple), list(simple)
            while frontier:
                new = []
                for r in frontier:
                    for s in range(self.n):
                        q = self.simple_reflect(s, r)
                        if all(c >= 0 for c in q) and q not in seen:
                            seen.add(q)
                            new.append(q)
                frontier = new
            self._roots = sorted(seen, key=lambda r: (sum(r), r))
        return self._roots

    def right_mult(self, M, s):
        cols = [list(c) for c in M]
        old = cols[s][:]
        cols[s] = [-v for v in old]
        for t in self.adj[s]:
            cols[t] = [u + v for u, v in zip(cols[t], old)]
        return tuple(tuple(c) for c in cols)

    def left_mult(self, s, M):
        return tuple(self.simple_reflect(s, c) for c in M)

    def word_matrix(self, word):
        M = self.identity
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
        return tuple(self.apply(A, B[j]) for j in range(self.n))

    def inverse(self, M):
        word = self.reduced_word(M)
        return self.word_matrix(word[::-1])

    @staticmethod
    def is_negative(v):
        return all(c <= 0 for c in v) and any(c < 0 for c in v)

    @staticmethod
    def is_positive(v):
        return all(c >= 0 for c in v) and any(c > 0 for c in v)

    def right_descents(self, M):
        return {s for s in range(self.n) if self.is_negative(M[s])}

    def left_descents(self, M):
        return self.right_descents(self.inverse(M))

    def length(self, M):
        return sum(1 for r in self.positive_roots() if self.is_negative(self.apply(M, r)))

    def inversions(self, M):
        return {r for r in self.positive_roots() if self.is_negative(self.apply(M, r))}

    def right_terminal(self, M):
        for s in self.right_descents(M):
            for t in self.adj[s]:
                if not self.is_positive(tuple(a + b for a, b in zip(M[s], M[t]))):
                    return False
        return True

    def terminal(self, M):
        return self.right_terminal(M) and self.right_terminal(self.inverse(M))

    def reduced_word(self, M):
        word = []
        while True:
            R = self.right_descents(M)
            if not R:
                break
            s = min(R)
            word.append(s)
            M = self.right_mult(M, s)
        return word[::-1]

    def reflection(self, beta):
        cols = []
        for j in range(self.n):
            e = tuple(int(i == j) for i in range(self.n))
            c = self.pairing(e, beta)
            cols.append(tuple(a - c * b for a, b in zip(e, beta)))
        return tuple(cols)

    def longest_element(self, J):
        M = self.identity
        while True:
            ascents = [s for s in J if not self.is_negative(M[s])]
            if not ascents:
                return M
            M = self.right_mult(M, ascents[0])

    def is_prefix(self, U, W):
        """U <= W in right weak order (W = U V length-additively) iff N(U^-1) subset N(W^-1)."""
        return self.inversions(self.inverse(U)) <= self.inversions(self.inverse(W))

    def is_suffix(self, V, W):
        return self.inversions(V) <= self.inversions(W)

    def suffixes(self, W):
        """All V with W = U V length-additively, by stripping left descents; returns {V: U}."""
        out = {W: self.identity}
        queue = deque([W])
        while queue:
            V = queue.popleft()
            for s in self.left_descents(V):
                Y = self.left_mult(s, V)
                if Y not in out:
                    out[Y] = self.mult(out[V], self.word_matrix([s]))
                    queue.append(Y)
        return out


def en(n):
    return RootSystem(n, [(i, i + 1) for i in range(n - 2)] + [(2, n - 1)])


def rank_of(rows):
    rows = [[Fraction(x) for x in r] for r in rows]
    rk = 0
    for c in range(len(rows[0]) if rows else 0):
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


# ---------- Gern's D_m ----------
def gern_interval(i, j):
    """Gern Definition 2.3.1."""
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
    """Bracket form of w_n, Gern Lemma 2.3.4 (n even)."""
    assert n % 2 == 0
    k = n // 2 - 2
    word = []
    for i in range(2, n + 1, 2):
        word += gern_interval(i, 0)
    for i in range(k + 1):
        word += gern_interval(n - k + i, n - 2 * k + 2 * i)
    return word


def gern_to_en(label, n):
    return 1 if label == 1 else (n - 1 if label == 2 else label - 1)


def signed_perm(word, m):
    """One-line signed permutation of a Gern word (left action, composed left to right)."""
    def gen(label):
        p = list(range(1, m + 1))
        if label == 1:
            p[0], p[1] = -2, -1
        else:
            p[label - 2], p[label - 1] = label, label - 1
        return p

    w = list(range(1, m + 1))
    for s in word:
        g = gen(s)
        w = [(1 if a > 0 else -1) * w[abs(a) - 1] for a in g]
    return tuple(w)


def gern_cor_2_2_19(n, m):
    """Signed permutation of w_n inside D_m, Gern Corollary 2.2.19 (n even)."""
    w = []
    for i in range(1, m + 1):
        if i > n:
            w.append(i)
        elif i == 1:
            w.append((-1) ** (n // 2))
        elif i % 2 == 1:
            w.append(i)
        else:
            w.append(-(n + 2 - i))
    return tuple(w)


def enumerate_bad_dm(m):
    """All noncommuting terminal elements of D_m as signed permutations (Gern labels)."""
    simple = {1: {1: 1, 2: 1}}
    for i in range(2, m + 1):
        simple[i] = {i: 1, i - 1: -1}
    adj = {i: set() for i in range(1, m + 1)}
    for a, b in [(1, 3), (2, 3)] + [(i, i + 1) for i in range(3, m)]:
        adj[a].add(b)
        adj[b].add(a)
    pos_roots = []
    for i in range(1, m + 1):
        for j in range(i + 1, m + 1):
            pos_roots.append({j: 1, i: -1})
            pos_roots.append({i: 1, j: 1})

    def image(w, root):
        out = {}
        for i, c in root.items():
            j = w[i - 1]
            out[abs(j)] = out.get(abs(j), 0) + (c if j > 0 else -c)
        return out

    def positive(r):
        items = [(k, v) for k, v in r.items() if v]
        (a, ca), (b, cb) = items
        if ca > 0 and cb > 0:
            return True
        if ca < 0 and cb < 0:
            return False
        hi, lo = (a, b) if ca > 0 else (b, a)
        return hi > lo

    def descents(w):
        return {s for s in range(1, m + 1) if not positive(image(w, simple[s]))}

    def inv(w):
        out = [0] * m
        for i, j in enumerate(w, 1):
            out[abs(j) - 1] = i if j > 0 else -i
        return tuple(out)

    def right_terminal(w):
        for s in descents(w):
            for t in adj[s]:
                r = {k: simple[s].get(k, 0) + simple[t].get(k, 0) for k in set(simple[s]) | set(simple[t])}
                if not positive(image(w, {k: v for k, v in r.items() if v})):
                    return False
        return True

    def length(w):
        return sum(1 for r in pos_roots if not positive(image(w, r)))

    def commuting_product(w):
        R = descents(w)
        if any(t in R for s in R for t in adj[s]):
            return False
        return length(w) == len(R) and signed_perm(sorted(R), m) == w

    bad, count = [], 0
    for perm in permutations(range(1, m + 1)):
        for signs in product((1, -1), repeat=m):
            if sum(1 for s in signs if s < 0) % 2:
                continue
            w = tuple(s * p for s, p in zip(signs, perm))
            count += 1
            if right_terminal(w) and right_terminal(inv(w)) and not commuting_product(w):
                bad.append((length(w), w))
    assert count == 2 ** (m - 1) * __import__('math').factorial(m)
    return sorted(bad)


def gern_predicted(m):
    """Gern Theorem 2.3.6: bad elements of D_m are w_k u, k even, 4 <= k <= m, u a commuting
    product of generators neither in nor adjacent to supp w_k = {1,...,k}."""
    out = {}
    for k in range(4, m + 1, 2):
        rest = list(range(k + 2, m + 1))
        for r in range(len(rest) + 1):
            for U in combinations(rest, r):
                if any(abs(a - b) == 1 for a in U for b in U):
                    continue
                name = f'w{k}' + (''.join(f'*s{u}' for u in U))
                out[name] = signed_perm(gern_wn_word(k) + list(U), m)
    return out


failures = []


def expect(cond, what):
    if not cond:
        failures.append(what)
        print('EXPECTATION FAILED:', what, file=sys.stderr)


def part_gern_rows():
    report = {}
    for n in (4, 6, 8, 10):
        word = gern_wn_word(n)
        expect(signed_perm(word, n) == gern_cor_2_2_19(n, n), f'Lemma 2.3.4 word of w_{n} equals Corollary 2.2.19')
        expect(len(word) == (3 * n * n + 2 * n) // 8, f'l(w_{n}) = 3n^2/8 + n/4')
    for n in (6, 7, 8):
        rs = en(n)
        m = n - 1
        rows = {w: rs.word_matrix(parse(w)) for (nn, w, _) in TABLE if nn == n}
        matched, mapped = set(), []
        for name, gw in sorted(gern_predicted(m).items()):
            # gern_predicted works with signed permutations; rebuild the word for the label map
            k = int(name.split('*')[0][1:])
            U = [int(t[1:]) for t in name.split('*')[1:]]
            gword = gern_wn_word(k) + U
            M = rs.word_matrix([gern_to_en(s, n) for s in gword])
            hits = [w for w in rows if rows[w] == M]
            expect(len(hits) == 1, f'E{n}: Gern {name} matches exactly one table row')
            matched.update(hits)
            mapped.append({'gern_element': name, 'gern_word': word_str(gword),
                           'mapped_word': word_str([gern_to_en(s, n) for s in gword]),
                           'table_row': hits[0] if hits else None, 'length': rs.length(M)})
        unmatched = sorted((w for w in rows if w not in matched), key=len)
        report[f'E{n}'] = {'type_D_parabolic': f'D{m} on nodes {sorted({1, n - 1} | set(range(2, n - 1)))}',
                           'gern_bad_elements_mapped': mapped, 'rows_not_from_type_D': unmatched}
    expect(report['E6']['rows_not_from_type_D'] == [], 'E6: every row is Gern\'s')
    expect(report['E7']['rows_not_from_type_D'] == ['1325436210321432543621324356'], 'E7: only w_7 is new')
    expect(report['E8']['rows_not_from_type_D'] == ['1325437210321432543721324357', TABLE[-1][1]], 'E8: only w_7 (in the E7 parabolic) and w_8 are not from type D')
    # The E8 row of length 28 is w_7 relabelled 6 -> 7 (E7 parabolic on {0,...,5,7}).
    rs8 = en(8)
    e7_in_e8 = TABLE[4][1].replace('6', '7')
    expect(rs8.word_matrix(parse(e7_in_e8)) == rs8.word_matrix(parse(TABLE[9][1])), 'E8 row of length 28 equals w_7 under 6 -> 7')
    report['E8_row_28_is_w7_in_parabolic_0_5_7'] = True
    return report


def part_dm_enumeration():
    report = {}
    for m in (5, 6, 7):
        bad = enumerate_bad_dm(m)
        predicted = gern_predicted(m)
        expect({w for _, w in bad} == set(predicted.values()), f'D{m}: bad elements are exactly Gern\'s list')
        report[f'D{m}'] = {'elements': 2 ** (m - 1) * __import__('math').factorial(m), 'noncommuting_terminals': len(bad),
                           'lengths': [l for l, _ in bad], 'gern_names': sorted(predicted),
                           'equals_gern_theorem_2_3_6': {w for _, w in bad} == set(predicted.values())}
    return report


def part_structure():
    report = {}
    for name, (n, wstr, layers) in LAYERS.items():
        rs = en(n)
        W = rs.word_matrix(parse(wstr))
        full = [s for layer in layers for s in layer]
        expect(rs.word_matrix(full) == W, f'{name}: layered word gives the table element')
        expect(len(full) == rs.length(W) == len(wstr), f'{name}: layered word is reduced')
        expect(layers == layers[::-1], f'{name}: layers are palindromic')
        colour = {j: j % 2 for j in range(n - 1)}
        colour[n - 1] = 1  # node n-1 is attached to node 2, so it has the colour of node 3
        classes = [{colour[s] for s in layer} for layer in layers]
        expect(all(len(c) == 1 for c in classes) and all(classes[i] != classes[i + 1] for i in range(len(classes) - 1)),
               f'{name}: layers alternate between the colour classes')
        expect(all(t not in rs.adj[s] for layer in layers for s in layer for t in layer), f'{name}: each layer is an independent set')
        L, R = rs.left_descents(W), rs.right_descents(W)
        table_descents = {frozenset(map(int, d)) for (nn, w, d) in TABLE if w == wstr}
        expect(table_descents == {frozenset(L)} and L == R, f'{name}: L = R = the common descent set of Table 1')
        expect(set(layers[0]) <= L, f'{name}: first layer is contained in L')
        expect(rs.mult(W, W) == rs.identity, f'{name}: involution')
        expect(rs.terminal(W), f'{name}: terminal')
        negated = [r for r in rs.positive_roots() if rs.apply(W, r) == tuple(-c for c in r)]
        fixed = [r for r in rs.positive_roots() if rs.apply(W, r) == r]
        expect([word_str(r) for r in negated] == EXPECTED_NEGATED[name], f'{name}: negated positive roots')
        expect(all(rs.pairing(a, b) == 0 for a in negated for b in negated if a != b), f'{name}: negated roots mutually orthogonal')
        P = rs.identity
        for r in negated:
            P = rs.mult(P, rs.reflection(r))
        expect(P == W, f'{name}: product of the reflections in the negated roots')
        MI = [[W[j][i] + int(i == j) for j in range(n)] for i in range(n)]
        expect(n - rank_of(MI) == len(negated), f'{name}: (-1)-eigenspace dimension equals the number of negated roots')
        # middle layer K and prefix P: W = P i(K) P^-1 with P(alpha_k), k in K, the negated roots
        half = layers[:len(layers) // 2]
        K = layers[len(layers) // 2]
        Pm = rs.word_matrix([s for layer in half for s in layer])
        images = sorted(rs.apply(Pm, tuple(int(i == k) for i in range(n))) for k in K)
        expect(images == sorted(negated), f'{name}: images of the middle layer are the negated roots')
        entry = {'type': f'E{n}', 'word': wstr, 'length': len(wstr), 'layers': layers, 'L_equals_R': sorted(L), 'first_layer': layers[0],
                 'involution': True, 'terminal': True, 'negated_positive_roots': [word_str(r) for r in negated],
                 'orthogonal_reflection_count': len(negated), 'fixed_positive_roots': len(fixed)}
        if fixed:
            span = rank_of([list(r) for r in fixed])
            entry['fixed_root_subsystem'] = {'rank': span, 'positive_roots': len(fixed)}
            if name == 'w8':
                expect(span == 6 and len(fixed) == 30, 'w8: fixed subsystem has rank 6 and 30 positive roots (type D6)')
                entry['fixed_root_subsystem']['type'] = 'D6'
        report[name] = entry
    return report


def part_chain_and_factorizations():
    report = {}
    # Chain in E8 and E7 (right weak order = prefix order; all four are involutions, so prefix = suffix).
    rs8 = en(8)
    words8 = {'w4': '1327213', 'w6': '132543721324357', 'w7': '1325437210321432543721324357', 'w8': TABLE[-1][1]}
    mats8 = {k: rs8.word_matrix(parse(v)) for k, v in words8.items()}
    chain8 = []
    for a, b in [('w4', 'w6'), ('w6', 'w7'), ('w7', 'w8')]:
        pre, suf = rs8.is_prefix(mats8[a], mats8[b]), rs8.is_suffix(mats8[a], mats8[b])
        expect(pre and suf, f'E8: {a} <= {b} in right (and left) weak order')
        chain8.append({'lower': a, 'upper': b, 'prefix': pre, 'suffix': suf})
    rs7 = en(7)
    words7 = {'w4': '1326213', 'w6': '132543621324356', 'w7': '1325436210321432543621324356'}
    mats7 = {k: rs7.word_matrix(parse(v)) for k, v in words7.items()}
    chain7 = []
    for a, b in [('w4', 'w6'), ('w6', 'w7')]:
        pre, suf = rs7.is_prefix(mats7[a], mats7[b]), rs7.is_suffix(mats7[a], mats7[b])
        expect(pre and suf, f'E7: {a} <= {b} in right weak order')
        chain7.append({'lower': a, 'upper': b, 'prefix': pre, 'suffix': suf})
    report['right_weak_order_chain'] = {'E8': {'words': words8, 'relations': chain8}, 'E7': {'words': words7, 'relations': chain7}}
    # Longest-parabolic inner factors and the w_0(J) factorizations.
    factors = {}
    for name, (n, wstr, _) in LAYERS.items():
        rs = en(n)
        W = rs.word_matrix(parse(wstr))
        suff = rs.suffixes(W)
        best, bestJ = 0, None
        for V in suff:
            Lv = tuple(sorted(rs.left_descents(V)))
            if Lv:
                l = rs.length(rs.longest_element(list(Lv)))
                if l > best:
                    best, bestJ = l, Lv
        entry = {'type': f'E{n}', 'right_weak_interval_size': len(suff), 'max_inner_longest_parabolic_length': best, 'achieved_by_J': list(bestJ)}
        I = sorted(rs.left_descents(W))
        J = sorted(set(I) | {0})
        w0J = rs.longest_element(J)
        hit = next((V for V in suff if set(J) <= rs.left_descents(V)), None)
        if hit is not None:
            u = suff[hit]
            v = rs.mult(rs.inverse(w0J), hit)
            additive = rs.length(u) + rs.length(w0J) + rs.length(v) == rs.length(W)
            expect(rs.mult(rs.mult(u, w0J), v) == W and additive, f'{name}: w = u w_0(J) v length-additively')
            entry['w0_factorization'] = {'J': J, 'l_w0_J': rs.length(w0J), 'u': word_str(rs.reduced_word(u)),
                                         'w0_J': word_str(rs.reduced_word(w0J)), 'v': word_str(rs.reduced_word(v)),
                                         'lengths': [rs.length(u), rs.length(w0J), rs.length(v)], 'length_additive': additive,
                                         'J_edges': [[a, b] for a, b in rs.edges if a in J and b in J]}
            if name in EXPECTED_FACTORIZATION:
                eu, ew, ev = EXPECTED_FACTORIZATION[name]
                expect(rs.word_matrix(parse(eu)) == u and rs.word_matrix(parse(ew)) == w0J and rs.word_matrix(parse(ev)) == v,
                       f'{name}: factorization matches the manuscript\'s words')
                expect(rs.length(w0J) == 6 and len(entry['w0_factorization']['J_edges']) == 1, f'{name}: J = I + {{0}} is of type A2 x A1^3 with l(w_0(J)) = 6')
        factors[name] = entry
    expect(factors['w4']['max_inner_longest_parabolic_length'] == 3, 'w4: longest inner parabolic factor has length 3')
    expect(factors['w6']['max_inner_longest_parabolic_length'] == 5, 'w6: longest inner parabolic factor has length 5')
    expect(factors['w7']['max_inner_longest_parabolic_length'] == 6, 'w7: longest inner parabolic factor has length 6')
    expect(factors['w8']['max_inner_longest_parabolic_length'] == 6, 'w8: longest inner parabolic factor has length 6')
    report['inner_longest_parabolic_factors'] = factors
    return report


def main():
    report = {'conventions': 'E_n: chain 0-1-...-(n-2), node n-1 attached to node 2; a string of labels is the product of the simple reflections in the order written; columns of the matrix of w are w(alpha_j). Gern D_m labels: 1,2 attached to 3, chain 3-...-m; label map 1->1, 2->n-1, j->j-1.',
              'gern_plus_two': part_gern_rows(),
              'type_D_enumeration': part_dm_enumeration(),
              'structure': part_structure()}
    report.update(part_chain_and_factorizations())
    report['status'] = 'FAILED' if failures else 'passed'
    report['failed_expectations'] = failures
    out = ROOT / 'results/terminal-structure-certificate.json'
    out.write_text(json.dumps(report, indent=1) + '\n')
    print(json.dumps({'status': report['status'], 'type_D_ranks_enumerated': sorted(report['type_D_enumeration']),
                      'chain_E8': [r['lower'] + '<' + r['upper'] for r in report['right_weak_order_chain']['E8']['relations']],
                      'inner_factor_lengths': {k: v['max_inner_longest_parabolic_length'] for k, v in report['inner_longest_parabolic_factors'].items()}}))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
