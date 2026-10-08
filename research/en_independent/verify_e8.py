"""Independent integral-matrix verification of the two E8 coset certificates.

No imports from the enumeration engines. Checks root inversion lengths,
word equality across D7/E7 decompositions, complete FC-catalogue closure,
all eligible bottoms, and exact parabolic/coset Poincare products.
Run: python3 research/en_independent/verify_e8.py
"""
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
N = 8
EDGES = tuple((s, s + 1) for s in range(6)) + ((2, 7),)
ADJ = [tuple(t if s == i else s for s, t in EDGES if i in (s, t)) for i in range(N)]
IDENTITY = tuple(tuple(int(i == j) for i in range(N)) for j in range(N))

def right(a, s):
    cols = list(a)
    cols[s] = tuple(-v for v in a[s])
    for t in ADJ[s]:
        cols[t] = tuple(u + v for u, v in zip(a[t], a[s]))
    return tuple(cols)

def element(word):
    a = IDENTITY
    for s in word:
        a = right(a, s)
    return a

def positive(a):
    assert any(a) and (all(v >= 0 for v in a) or all(v <= 0 for v in a)), a
    return any(v > 0 for v in a)

def desc(a):
    return frozenset(s for s, col in enumerate(a) if not positive(col))

def mask(ds):
    return sum(1 << s for s in ds)

def independent(ds):
    return not any(s in ds and t in ds for s, t in EDGES)

def weak(a):
    # Equivalent to no suffix ts with s,t adjacent, tested without words.
    return all(positive(tuple(u + v for u, v in zip(a[s], a[t])))
               for s in desc(a) for t in ADJ[s])

def image(a, root):
    return tuple(sum(root[s] * a[s][k] for s in range(N)) for k in range(N))

def roots():
    found = set(IDENTITY)
    todo = list(IDENTITY)
    for a in todo:
        for s in range(N):
            b = list(a)
            b[s] = -a[s] + sum(a[t] for t in ADJ[s])
            b = tuple(b)
            if b not in found:
                found.add(b)
                todo.append(b)
    assert len(found) == 240
    assert all(positive(a) or positive(tuple(-v for v in a)) for a in found)
    return frozenset(found)

def leq(x, lx, w, lw):
    while x != w:
        if lx >= lw:
            return False
        s = min(desc(w))
        if s in desc(x):
            x, lx = right(x, s), lx - 1
        w, lw = right(w, s), lw - 1
    return True

VOLATILE_KEYS = {'seconds', 'elapsed_seconds'}

def strip_volatile(value):
    if isinstance(value, dict):
        return {k: strip_volatile(v) for k, v in value.items() if k not in VOLATILE_KEYS}
    if isinstance(value, list):
        return [strip_volatile(v) for v in value]
    return value

def content_hash(path):
    # Generated JSON certificates carry wall-clock timing fields; hash a canonical
    # copy without them so the audit is reproducible across machines.
    if path.suffix == '.json':
        data = json.dumps(strip_volatile(json.loads(path.read_text())), sort_keys=True, separators=(',', ':'))
        return hashlib.sha256(data.encode()).hexdigest()
    return hashlib.sha256(path.read_bytes()).hexdigest()

def convolution(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x * y
    return c

def poincare(exponents):
    a = [1]
    for m in exponents:
        a = convolution(a, [1] * (m + 1))
    return a

def main():
    paths = [ROOT / 'research/en_independent/e8-d7-terminals.json',
             ROOT / 'research/en_e8/e8-terminals.json',
             ROOT / 'research/en_independent/e8-fc.json']
    d7, e7, catalogue = [json.loads(p.read_text()) for p in paths]
    assert d7['complete'] and e7['complete'] and catalogue['complete']
    assert d7['parabolic_order'] * d7['cosets'] == 696729600
    assert e7['parabolic_order'] * e7['coset_count'] == 696729600
    assert d7['parabolic_histogram'] == poincare([1,3,5,7,9,11,6])
    assert convolution(d7['parabolic_histogram'], d7['coset_histogram']) == poincare([1,7,11,13,17,19,23,29])
    assert d7['full_right_terminals'] == e7['right_terminal_count'] == 2160
    assert d7['commuting_terminals'] == e7['commuting_terminals'] == 58
    assert sum(independent(set(s for s in range(N) if bits & (1 << s)))
               for bits in range(1 << N)) == 58
    all_roots = roots()
    positive_roots = [r for r in all_roots if positive(r)]
    def length(a):
        return sum(not positive(image(a, r)) for r in positive_roots)
    fc = {}
    for row in catalogue['elements']:
        a = element(row['word'])
        assert a not in fc and independent(desc(a))
        assert len(row['word']) == row['length']
        assert mask(desc(a)) == row['Rmask']
        assert mask(desc(element(reversed(row['word'])))) == row['Lmask']
        fc[a] = row
    assert len(fc) == 10846
    # Recursive FC criterion, together with a minimal-missing-element argument,
    # proves that the finite catalogue is both sound and complete.
    ascents = 0
    for a, row in fc.items():
        ds = desc(a)
        for s in ds:
            assert right(a, s) in fc
            assert fc[right(a, s)]['length'] == row['length'] - 1
        for s in range(N):
            if s in ds:
                continue
            ascents += 1
            b = right(a, s)
            valid = independent(desc(b)) and all(right(b, t) in fc for t in desc(b))
            assert valid == (b in fc)
    assert ascents == catalogue['ascents_tested']
    canonical = []
    for source in (d7, e7):
        records = {}
        for row in source['bad']:
            a = element(row['word'])
            ai = element(reversed(row['word']))
            assert a not in records and weak(a) and weak(ai) and a not in fc
            assert length(a) == row['length'] == len(row['word'])
            assert {image(a, r) for r in all_roots} == all_roots
            rd, ld = desc(a), desc(ai)
            assert mask(rd) == row['Rmask'] and mask(ld) == row['Lmask']
            eligible = []
            for x, info in fc.items():
                if info['Rmask'] & row['Rmask'] != row['Rmask'] or info['Lmask'] & row['Lmask'] != row['Lmask']:
                    continue
                if leq(x, info['length'], a, row['length']):
                    eligible.append((x, info['length']))
            supplied = [(element(b['word']), b['length']) for b in row['bottoms']]
            assert set(eligible) == set(supplied)
            assert len(eligible) == 1
            records[a] = (row['length'], frozenset(rd), frozenset(ld), tuple(eligible))
        canonical.append(records)
    assert canonical[0] == canonical[1]
    summary = {'status': 'passed', 'complete_E8_terminal_classification': True,
               'distinct_parabolic_methods': ['E7, 240 cosets', 'D7, 2160 cosets'],
               'matrix_words_and_inversion_lengths_agree': True,
               'terminal_count': 64, 'commuting_terminals': 58,
               'bad_terminal_lengths': sorted(v[0] for v in canonical[0].values()),
               'eligible_ranks': sorted(v[0] - v[3][0][1] for v in canonical[0].values()),
               'complete_FC_catalogue': len(fc), 'FC_ascents_verified': ascents,
               'parabolic_and_coset_poincare_products_match': True,
               'source_sha256': {str(p.relative_to(ROOT)): content_hash(p)
                                 for p in paths + [Path(__file__), ROOT / 'research/en_independent/e8_d7_cosets.cpp',
                                                  ROOT / 'research/en_independent/fc_catalogue.cpp',
                                                  ROOT / 'research/en_e8/parabolic_terminals.cpp']},
               'hashing': 'JSON inputs are hashed after removing the volatile keys ' + ', '.join(sorted(VOLATILE_KEYS)) + ' and re-serialising with sorted keys; other files are hashed byte for byte.',
               'priority': 'Not established by this computational audit.'}
    output = ROOT / 'results/e8-independent-audit.json'
    output.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()
