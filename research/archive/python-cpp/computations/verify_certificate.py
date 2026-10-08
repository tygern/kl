"""Verify saved data without invoking the KL computation routine.

Usage: python3 computations/verify_certificate.py [certificate.json]
Default: results/d6-recurrence-certificate.json

For every supplied KL dependency, checks the standard recurrence from the
polynomial records, recalculating Bruhat comparisons, the entire correction
index set, degree/positivity constraints, and all required references.
Separately checks q^d P(q^-1)=sum R_xz P_zw at the root.
"""
import json
from pathlib import Path
import sys
from coxeter import add, mul, trim
from sparse_kl import SparseCoxeter


def verify(filename):
    data = json.loads(Path(filename).read_text())
    g = SparseCoxeter(data['group']['type'], data['group']['rank'])
    elements = [tuple(w) for w in data['elements']]
    assert len(set(elements)) == len(elements)
    for w in elements:
        assert sorted(abs(a) for a in w) == list(range(1, len(w)+1))
        assert sum(a < 0 for a in w) % 2 == 0
    records = {(elements[r['x']], elements[r['w']]): r for r in data['records']}
    assert len(records) == len(data['records']) == data['number_records']

    def p(x, w):
        return tuple(records[x, w]['polynomial'])

    def recurrence(x, w, s):
        assert s in g.descents(w)
        v, xs = g.right(w, s), g.right(x, s)
        c = int(g.length(xs) < g.length(x))
        value = add((), p(xs, v), 1-c)
        value = add(value, p(x, v), c)
        terms = []
        for z in g.lower(v):
            delta = g.length(v)-g.length(z)
            if delta <= 0 or delta % 2 == 0 or s not in g.descents(z):
                continue
            polynomial = p(z, v)
            mu = polynomial[(delta-1)//2] if len(polynomial) > (delta-1)//2 else 0
            if mu and g.leq(x, z):
                exponent = (delta+1)//2
                terms.append((z, exponent, mu))
                value = add(value, p(x, z), exponent, -mu)
        return value, sorted(terms)

    comparable, strict, zero = 0, 0, 0
    for (x, w), row in records.items():
        polynomial = p(x, w)
        if x == w:
            assert polynomial == (1,)
            comparable += 1
        elif not g.leq(x, w):
            assert polynomial == ()
            zero += 1
        else:
            comparable += 1
            strict += 1
            assert polynomial[0] == 1 and min(polynomial) >= 0
            assert 2*(len(polynomial)-1) < g.length(w)-g.length(x)
            value, terms = recurrence(x, w, row['right_descent'])
            listed = sorted((elements[z], e, m) for z, e, m in row['correction_terms'])
            assert terms == listed
            assert polynomial == value
    x, w = (elements[i] for i in data['root'])
    for row in data['root_recurrences']:
        value, terms = recurrence(x, w, row['right_descent'])
        assert value == p(x, w)
        assert terms == sorted((elements[z], e, m) for z, e, m in row['correction_terms'])
    assert {row['right_descent'] for row in data['root_recurrences']} == set(g.descents(w))
    polynomial, delta = p(x, w), g.length(w)-g.length(x)
    lhs = trim([polynomial[delta-i] if 0 <= delta-i < len(polynomial) else 0 for i in range(delta+1)])
    rhs = ()
    for z in g.lower(w):
        if g.leq(x, z):
            rhs = add(rhs, mul(g.rpoly(x, z), p(z, w)))
    assert lhs == rhs
    # The verifier's KL evaluator must never have run.
    assert not g.kl_values
    result = dict(records=len(records), comparable_records=comparable,
                  strict_recurrence_checks=strict, incomparable_zero_checks=zero,
                  root_polynomial=polynomial, root_right_descent_checks=len(g.descents(w)),
                  root_R_reciprocity=True, evaluator_called=False)
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    default = Path(__file__).resolve().parents[1] / 'results' / 'd6-recurrence-certificate.json'
    verify(sys.argv[1] if len(sys.argv) > 1 else default)
