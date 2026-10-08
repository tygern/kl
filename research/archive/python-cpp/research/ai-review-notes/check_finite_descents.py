"""Check all eleven finite table rows; no imports from enumeration engines.

python3 research/ai-review-notes/check_finite_descents.py
The matching is a short certificate for alpha(support) <= |support|-|matching|.
"""
import json
from pathlib import Path

ROWS = [
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


def check(n, word_string, lower_string):
    word = list(map(int, word_string))
    support = set(word)
    edges = [(i, i+1) for i in range(n-2)] + [(2, n-1)]
    support_edges = [e for e in edges if set(e) <= support]
    independent_set = set(map(int, lower_string))
    assert not any(set(e) <= independent_set for e in edges)

    def multiply(seq):
        columns = [[int(i == j) for i in range(n)] for j in range(n)]
        for s in seq:
            old = columns[s][:]
            columns[s] = [-v for v in old]
            for a, b in edges:
                if s in (a, b):
                    t = b if s == a else a
                    columns[t] = [u+v for u, v in zip(columns[t], old)]
        return columns

    def desc(columns):
        assert all(any(c) and (all(v >= 0 for v in c) or all(v <= 0 for v in c)) for c in columns)
        return {j for j, c in enumerate(columns) if any(v < 0 for v in c)}

    assert desc(multiply(word)) == desc(multiply(reversed(word))) == independent_set
    matchings = []
    for bits in range(1 << len(support_edges)):
        chosen = [e for i, e in enumerate(support_edges) if bits & (1 << i)]
        vertices = [s for e in chosen for s in e]
        if len(vertices) == len(set(vertices)):
            matchings.append(chosen)
    matching = max(matchings, key=len)
    assert len(support) - len(matching) == len(independent_set)
    return dict(type=f'E{n}', word=word_string, length=len(word),
                support=sorted(support), common_descents=sorted(independent_set),
                matching=matching, independence_number=len(independent_set))


if __name__ == '__main__':
    result = [check(*row) for row in ROWS]
    target = Path(__file__).with_name('finite-descents.json')
    target.write_text(json.dumps(result, indent=2) + '\n')
    for row in result:
        print(row['type'], row['length'], row['support'], row['common_descents'], row['matching'])
