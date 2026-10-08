"""Targeted signed-permutation KL computation without enumerating the group.

Independent length formula and lifting-property Bruhat comparison provide
a second check on the finite BFS/subword implementation in coxeter.py.
"""
from functools import cache
from pathlib import Path
from coxeter import add, mul, trim


class SparseCoxeter:
    def __init__(self, kind, rank):
        self.kind, self.rank = kind, rank
        self.kl_values = {}

    @cache
    def length(self, w):
        inversions = sum(w[i] > w[j] for i in range(len(w)) for j in range(i+1, len(w)))
        if self.kind == 'A':
            return inversions
        negative_sum_pairs = sum(-w[i] > w[j] for i in range(len(w)) for j in range(i+1, len(w)))
        return inversions+negative_sum_pairs+(sum(x < 0 for x in w) if self.kind == 'B' else 0)

    @cache
    def right(self, w, s):
        a = list(w)
        if self.kind == 'A':
            a[s], a[s+1] = a[s+1], a[s]
        elif s:
            a[s-1], a[s] = a[s], a[s-1]
        elif self.kind == 'B':
            a[0] *= -1
        else:
            a[0], a[1] = -a[1], -a[0]
        return tuple(a)

    @cache
    def descents(self, w):
        return tuple(s for s in range(self.rank) if self.length(self.right(w, s)) < self.length(w))

    @cache
    def leq(self, x, w):
        if x == w:
            return True
        if self.length(x) >= self.length(w):
            return False
        s = self.descents(w)[0]
        xs, ws = self.right(x, s), self.right(w, s)
        return self.leq(xs if self.length(xs) < self.length(x) else x, ws)

    @cache
    def lower(self, w):
        if not self.length(w):
            return frozenset([w])
        s = self.descents(w)[0]
        base = self.lower(self.right(w, s))
        return base | {self.right(x, s) for x in base}

    @cache
    def corrections(self, w, s):
        v = self.right(w, s)
        result = []
        for z in self.lower(v):
            delta = self.length(v)-self.length(z)
            if delta <= 0 or delta % 2 == 0 or s not in self.descents(z):
                continue
            p = self.kl(z, v)
            mu = p[(delta-1)//2] if len(p) > (delta-1)//2 else 0
            if mu:
                result.append((z, (delta+1)//2, mu))
        return tuple(result)

    @cache
    def kl(self, x, w):
        if x == w:
            p = (1,)
        elif not self.leq(x, w):
            p = ()
        else:
            p = self.with_descent(x, w, self.descents(w)[0])
        self.kl_values[x, w] = p
        return p

    def with_descent(self, x, w, s):
        v, xs = self.right(w, s), self.right(x, s)
        c = int(self.length(xs) < self.length(x))
        p = add((), self.kl(xs, v), 1-c)
        p = add(p, self.kl(x, v), c)
        for z, exponent, mu in self.corrections(w, s):
            if self.leq(x, z):
                p = add(p, self.kl(x, z), exponent, -mu)
        return p

    @cache
    def rpoly(self, x, w):
        if x == w:
            return (1,)
        if not self.leq(x, w):
            return ()
        s = self.descents(w)[0]
        v, xs = self.right(w, s), self.right(x, s)
        if self.length(xs) < self.length(x):
            return self.rpoly(xs, v)
        p = add((), self.rpoly(xs, v), 1)
        p = add(p, self.rpoly(x, v), 1)
        return add(p, self.rpoly(x, v), scale=-1)

    def verify_reciprocity(self, x, w):
        delta = self.length(w)-self.length(x)
        p = self.kl(x, w)
        lhs = trim([p[delta-i] if 0 <= delta-i < len(p) else 0 for i in range(delta+1)])
        rhs = ()
        for z in self.lower(w):
            if self.leq(x, z):
                rhs = add(rhs, mul(self.rpoly(x, z), self.kl(z, w)))
        assert lhs == rhs, (lhs, rhs)


def bad_d6(certificate_path=None):
    from time import perf_counter
    start = perf_counter()
    g = SparseCoxeter('D', 6)
    x = (-1, -2, 4, 3, 6, 5)
    w = (-1, -6, 3, -4, 5, -2)
    p = g.kl(x, w)
    print('D6 KL', p, 'seconds', round(perf_counter()-start, 3), flush=True)
    for s in g.descents(w):
        assert g.with_descent(x, w, s) == p
    g.verify_reciprocity(x, w)
    interval = [z for z in g.lower(w) if g.leq(x, z)]
    from collections import Counter
    ranks = Counter(g.length(z)-g.length(x) for z in interval)
    if certificate_path is not None:
        write_dependency_certificate(g, x, w, certificate_path)
    return dict(type='D6', bottom=x, top=w, bottom_length=g.length(x),
                top_length=g.length(w), interval_rank=g.length(w)-g.length(x),
                lower_ideal_size=len(g.lower(w)), interval_size=len(interval),
                rank_vector=[ranks[i] for i in range(g.length(w)-g.length(x)+1)],
                polynomial_coefficients=p, mu=p[5] if len(p) > 5 else 0,
                checks=['independent signed inversion length formula',
                        'lifting-property Bruhat comparisons',
                        'KL recurrence with every right descent of top',
                        'independent R-polynomial reciprocity identity'],
                elapsed_seconds=round(perf_counter()-start, 3),
                kl_cache_entries=len(g.kl_values),
                recurrence_certificate=relative_label(certificate_path) if certificate_path else None)


def relative_label(path):
    """Repository-relative label for a results file (no machine-specific paths in JSON)."""
    path = Path(path).resolve()
    root = Path(__file__).resolve().parents[1]
    return str(path.relative_to(root)) if root in path.parents else path.name


def write_dependency_certificate(g, x, w, filename):
    """Plain JSON records for all evaluated KL dependencies, including zeros.

    The verifier checks missing/extra correction terms itself, so it does not
    rely on the listed correction terms being exhaustive.
    """
    import json
    from pathlib import Path
    elements = sorted({z for pair in g.kl_values for z in pair} | set(g.lower(w)))
    indices = {z: i for i, z in enumerate(elements)}
    rows = []
    for (a, b), p in sorted(g.kl_values.items()):
        row = dict(x=indices[a], w=indices[b], polynomial=p)
        if a != b and g.leq(a, b):
            s = g.descents(b)[0]
            row['right_descent'] = s
            row['correction_terms'] = sorted([
                [indices[z], exponent, mu]
                for z, exponent, mu in g.corrections(b, s) if g.leq(a, z)])
        rows.append(row)
    root_recurrences = []
    for s in g.descents(w):
        root_recurrences.append(dict(right_descent=s, correction_terms=sorted([
            [indices[z], exponent, mu]
            for z, exponent, mu in g.corrections(w, s) if g.leq(x, z)])))
    data = dict(format_version=1, group=dict(type=g.kind, rank=g.rank),
                conventions='As documented in computations/coxeter.py. Polynomials use ascending powers of q.',
                elements=elements, root=[indices[x], indices[w]],
                records=rows, root_recurrences=root_recurrences,
                number_records=len(rows))
    Path(filename).write_text(json.dumps(data, separators=(',', ':'))+'\n')


if __name__ == '__main__':
    import json
    print(json.dumps(bad_d6(), indent=2))
