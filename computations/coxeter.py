"""Exact small finite Coxeter groups, equal-parameter KL polynomials, intervals.

Python standard library only. Permutations are one-line images; words multiply
on the right. A_r acts on r+1 letters, B_r and D_r on r signed letters.
A: s_i swaps positions i,i+1 (zero based). B: s_0 negates position 0,
s_i swaps positions i-1,i for i>0. D: s_0 maps (a,b) to (-b,-a),
and s_i swaps positions i-1,i for i>0. Thus D has fork at s_2.
"""
from collections import Counter, defaultdict, deque
from functools import cache


def trim(p):
    p = list(p)
    while p and not p[-1]:
        p.pop()
    return tuple(p)


def add(p, q, shift=0, scale=1):
    a = list(p) + [0] * max(0, len(q) + shift - len(p))
    for i, c in enumerate(q):
        a[i + shift] += scale * c
    return trim(a)


def mul(p, q):
    a = ()
    for i, c in enumerate(p):
        a = add(a, q, i, c)
    return a


def compose_permutations(u, v):
    """Composition u after v, for signed one-line permutations."""
    return tuple((1 if a > 0 else -1)*u[abs(a)-1] for a in v)


def permutation_order(w):
    identity = tuple(range(1, len(w)+1))
    v = identity
    for n in range(1, 10001):
        v = compose_permutations(v, w)
        if v == identity:
            return n
    raise ValueError('permutation order unexpectedly large')


class Coxeter:
    def __init__(self, kind, rank):
        assert kind in 'ABD'
        assert rank >= (3 if kind == 'D' else 1)
        self.kind, self.rank = kind, rank
        self.name = f'{kind}{rank}'
        e = tuple(range(1, rank + 2 if kind == 'A' else rank + 1))
        self.elements, self.index = [e], {e: 0}
        self.length, self.words, self.right = [0], [()], []
        for w in self.elements:
            wi = self.index[w]
            row = []
            for s in range(rank):
                v = self.multiply_simple(w, s)
                if v not in self.index:
                    self.index[v] = len(self.elements)
                    self.elements.append(v)
                    self.length.append(self.length[wi] + 1)
                    self.words.append(self.words[wi] + (s,))
                row.append(self.index[v])
            self.right.append(row)
        self.order = len(self.elements)
        self.descents = [tuple(s for s, v in enumerate(row)
                              if self.length[v] < self.length[w])
                         for w, row in enumerate(self.right)]
        self.lower = [set() for _ in self.elements]
        self.lower[0] = {0}
        self.fc = [True] * self.order
        self.braid_witness = [None] * self.order
        self.coxeter_matrix = [[self.generator_order(s, t)
                               for t in range(rank)] for s in range(rank)]
        for w in range(1, self.order):
            s = self.descents[w][0]
            v = self.right[w][s]
            self.lower[w] = self.lower[v] | {self.right[x][s] for x in self.lower[v]}
            self.fc[w] = self.is_fc_recursive(w)
        self.covers_down = [tuple(x for x in sorted(self.lower[w])
                                  if self.length[x] + 1 == self.length[w])
                            for w in range(self.order)]
        self.upper = [set() for _ in self.elements]
        for w, lower in enumerate(self.lower):
            for x in lower:
                self.upper[x].add(w)
        self.inverses = [self.index[self.inverse(w)] for w in self.elements]

    def multiply_simple(self, w, s):
        a = list(w)
        if self.kind == 'A':
            a[s], a[s+1] = a[s+1], a[s]
        elif s > 0:
            a[s-1], a[s] = a[s], a[s-1]
        elif self.kind == 'B':
            a[0] *= -1
        else:
            a[0], a[1] = -a[1], -a[0]
        return tuple(a)

    @staticmethod
    def inverse(w):
        a = [0] * len(w)
        for i, v in enumerate(w):
            a[abs(v)-1] = (i+1) * (1 if v > 0 else -1)
        return tuple(a)

    def generator_order(self, s, t):
        if s == t:
            return 1
        w = 0
        for k in range(1, 7):
            w = self.right[self.right[w][s]][t]
            if not w:
                return k
        raise AssertionError('unexpected Coxeter relation')

    def is_fc_recursive(self, w):
        # Every non-FC word either has a non-FC prefix of length l(w)-1,
        # or has a full noncommuting braid as its suffix.
        for s in self.descents[w]:
            v = self.right[w][s]
            if not self.fc[v]:
                self.braid_witness[w] = self.braid_witness[v] + (s,)
                return False
        for s in range(self.rank):
            for t in range(s+1, self.rank):
                m = self.coxeter_matrix[s][t]
                if m < 3:
                    continue
                suffix = tuple(s if k % 2 == 0 else t for k in range(m))
                v = w
                for a in reversed(suffix):
                    v = self.right[v][a]
                if self.length[v] + m == self.length[w]:
                    self.braid_witness[w] = self.words[v] + suffix
                    return False
        return True

    def word_element(self, word):
        w = 0
        for s in word:
            w = self.right[w][s]
        return w

    def interval(self, x, w):
        return tuple(sorted(self.lower[w] & self.upper[x]))

    def interval_edges(self, nodes):
        nodes = set(nodes)
        return {(x, w) for w in nodes for x in self.covers_down[w] if x in nodes}

    def describe(self, x):
        return dict(id=x, one_line=self.elements[x], length=self.length[x],
                    reduced_word=self.words[x], fully_commutative=self.fc[x],
                    non_fc_reduced_word=self.braid_witness[x])

    @cache
    def kl(self, x, w):
        if x not in self.lower[w]:
            return ()
        if x == w:
            return (1,)
        return self.kl_with_descent(x, w, self.descents[w][0])

    def kl_with_descent(self, x, w, s):
        assert s in self.descents[w]
        v, xs = self.right[w][s], self.right[x][s]
        c = int(self.length[xs] < self.length[x])
        p = add((), self.kl(xs, v), 1-c)
        p = add(p, self.kl(x, v), c)
        for z in self.lower[v]:
            if z == v or x not in self.lower[z] or s not in self.descents[z]:
                continue
            mu = self.mu(z, v)
            if mu:
                p = add(p, self.kl(x, z), (self.length[w]-self.length[z])//2, -mu)
        return p

    def mu(self, x, w):
        d = self.length[w] - self.length[x]
        if d <= 0 or d % 2 == 0:
            return 0
        p = self.kl(x, w)
        return p[(d-1)//2] if (d-1)//2 < len(p) else 0

    @cache
    def rpoly(self, x, w):
        if x not in self.lower[w]:
            return ()
        if x == w:
            return (1,)
        s = self.descents[w][0]
        v, xs = self.right[w][s], self.right[x][s]
        if self.length[xs] < self.length[x]:
            return self.rpoly(xs, v)
        p = add((), self.rpoly(xs, v), 1)
        p = add(p, self.rpoly(x, v), 1)
        return add(p, self.rpoly(x, v), scale=-1)

    def verify_reciprocity(self, x, w):
        d = self.length[w]-self.length[x]
        p = self.kl(x, w)
        lhs = trim([p[d-i] if 0 <= d-i < len(p) else 0 for i in range(d+1)])
        rhs = ()
        for z in self.interval(x, w):
            rhs = add(rhs, mul(self.rpoly(x, z), self.kl(z, w)))
        assert lhs == rhs, (self.name, x, w, lhs, rhs)

    def validate(self):
        import math
        expected = math.factorial(self.rank+1) if self.kind == 'A' else (
            2**(self.rank-(self.kind == 'D')) * math.factorial(self.rank))
        assert self.order == expected
        count, fc_count, mus = 0, 0, Counter()
        w0 = self.length.index(max(self.length))
        for w in range(self.order):
            assert all(abs(self.length[x]-self.length[w]) == 1 for x in self.right[w])
            if not self.fc[w]:
                assert self.word_element(self.braid_witness[w]) == w
                assert len(self.braid_witness[w]) == self.length[w]
            for x in self.lower[w]:
                count += 1
                p, d = self.kl(x, w), self.length[w]-self.length[x]
                assert p[0] == 1 and min(p) >= 0
                assert x == w or 2*(len(p)-1) < d
                assert p == self.kl(self.inverses[x], self.inverses[w])
                if self.fc[x] and x != w:
                    fc_count += 1
                    mus[self.mu(x, w)] += 1
            assert self.kl(w, w0) == (1,)
        return dict(group_order=self.order, max_length=max(self.length),
                    fully_commutative_elements=sum(self.fc),
                    comparable_pairs_including_equal=count,
                    strict_pairs_with_fc_bottom=fc_count,
                    fc_bottom_mu_distribution=dict(sorted(mus.items())),
                    checks=['group order', 'length change under generators',
                            'non-FC reduced braid witnesses', 'KL positivity',
                            'KL degree bound and constant term',
                            'KL invariance under inverse', 'P(x,w0)=1'])


def poset_data(group, x, w):
    nodes = group.interval(x, w)
    edges = group.interval_edges(nodes)
    ranks = {v: group.length[v]-group.length[x] for v in nodes}
    down, up = {v: set() for v in nodes}, {v: set() for v in nodes}
    for a, b in edges:
        down[b].add(a)
        up[a].add(b)
    # Canonical color refinement, using signatures as tuples of integers.
    colors = {v: ranks[v] for v in nodes}
    for _ in range(len(nodes)):
        sig = {v: (colors[v], tuple(sorted(colors[a] for a in down[v])),
                   tuple(sorted(colors[b] for b in up[v]))) for v in nodes}
        ids = {s: i for i, s in enumerate(sorted(set(sig.values())))}
        new = {v: ids[sig[v]] for v in nodes}
        if len(set(new.values())) == len(set(colors.values())):
            colors = new
            break
        colors = new
    fingerprint = (tuple(Counter(ranks.values())[i] for i in range(max(ranks.values())+1)),
                   tuple(sorted(Counter((colors[v], tuple(sorted(colors[a] for a in down[v])),
                                         tuple(sorted(colors[b] for b in up[v])))
                                        for v in nodes).items())))
    return dict(nodes=nodes, edges=edges, ranks=ranks, down=down, up=up,
                colors=colors, fingerprint=fingerprint)


def isomorphism(a, b):
    """Exact directed-Hasse graph isomorphism after deterministic refinement."""
    if a['fingerprint'] != b['fingerprint']:
        return None
    classes = defaultdict(list)
    for v in b['nodes']:
        classes[b['colors'][v]].append(v)
    mapping, used = {}, set()

    def candidates(u):
        return [v for v in classes[a['colors'][u]] if v not in used
                and all(((i in a['down'][u]) == (j in b['down'][v]) and
                         (i in a['up'][u]) == (j in b['up'][v]))
                        for i, j in mapping.items())]

    def visit():
        if len(mapping) == len(a['nodes']):
            return dict(mapping)
        choices = [(candidates(u), u) for u in a['nodes'] if u not in mapping]
        choices.sort(key=lambda cu: (len(cu[0]), -len((a['down'][cu[1]] | a['up'][cu[1]]) & mapping.keys()), cu[1]))
        cand, u = choices[0]
        for v in cand:
            mapping[u] = v
            used.add(v)
            result = visit()
            if result is not None:
                return result
            used.remove(v)
            del mapping[u]
        return None
    return visit()


def certify(source, x, w, target, y, v, mapping):
    """Verify bijection/order directly, independently of isomorphism search."""
    a, b = source.interval(x, w), target.interval(y, v)
    assert set(mapping) == set(a) and set(mapping.values()) == set(b)
    for i in a:
        for j in a:
            assert (i in source.lower[j]) == (mapping[i] in target.lower[mapping[j]])
            if i in source.lower[j]:
                assert source.kl(i, j) == target.kl(mapping[i], mapping[j])
                source.verify_reciprocity(i, j)
                target.verify_reciprocity(mapping[i], mapping[j])
                for s in source.descents[j]:
                    if i != j:
                        assert source.kl(i, j) == source.kl_with_descent(i, j, s)
                for s in target.descents[mapping[j]]:
                    if i != j:
                        assert target.kl(mapping[i], mapping[j]) == target.kl_with_descent(mapping[i], mapping[j], s)
    return dict(source=dict(type=source.name, bottom=source.describe(x), top=source.describe(w)),
                target=dict(type=target.name, bottom=target.describe(y), top=target.describe(v)),
                rank=source.length[w]-source.length[x],
                rank_vector=poset_data(source, x, w)['fingerprint'][0],
                polynomial_coefficients=source.kl(x, w), mu=source.mu(x, w),
                map=[dict(source=source.describe(i), target=target.describe(mapping[i])) for i in a],
                source_cover_edges=sorted(source.interval_edges(a)),
                target_cover_edges=sorted(target.interval_edges(b)),
                checks=['bijection', 'all pairwise order relations',
                        'all subinterval KL polynomials agree',
                        'all subinterval R-polynomial reciprocity identities',
                        'all right-descent choices in the KL recurrence'])


def non_type_a_coset_certificate(group, x, w):
    """A pair of cover reflections of product order 4 excludes an A parabolic coset.

    If all interval elements lie in gW_J, each right cover label u^-1 v lies
    in W_J. Any two reflections in a type-A Coxeter group have product order
    1, 2 or 3. The analogous argument uses vu^-1 for right parabolic cosets.
    """
    edges = sorted(group.interval_edges(group.interval(x, w)))
    result = {}
    for side in ['right', 'left']:
        labels = {}
        for u, v in edges:
            ui = group.inverse(group.elements[u])
            t = compose_permutations(ui, group.elements[v]) if side == 'right' else compose_permutations(group.elements[v], ui)
            assert permutation_order(t) == 2
            labels.setdefault(t, (u, v))
        for t in sorted(labels):
            for u in sorted(labels):
                if permutation_order(compose_permutations(t, u)) == 4:
                    result[side+'_cover_labels'] = dict(
                        first=dict(reflection=t, cover=labels[t]),
                        second=dict(reflection=u, cover=labels[u]), product_order=4)
                    break
            if side+'_cover_labels' in result:
                break
    return result
