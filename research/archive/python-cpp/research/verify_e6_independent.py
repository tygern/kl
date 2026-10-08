"""Independent exact E6 terminal certificate; standard library only.

This does not import any project implementation. Elements are integral matrices
(stored by columns). Enumeration uses the faithful geometric representation.
Run: python3 research/verify_e6_independent.py
"""
import hashlib
import json
from collections import Counter
from pathlib import Path

N = 6
EDGES = ((0, 1), (1, 2), (2, 3), (3, 4), (2, 5))
ADJ = [tuple(t if s == i else s for s, t in EDGES if i in (s, t)) for i in range(N)]
IDENTITY = tuple(tuple(int(i == j) for i in range(N)) for j in range(N))

def right(a, s):
    cols = list(a)
    cols[s] = tuple(-v for v in a[s])
    for t in ADJ[s]:
        cols[t] = tuple(u + v for u, v in zip(a[t], a[s]))
    return tuple(cols)

def negative(v):
    assert all(c >= 0 for c in v) or all(c <= 0 for c in v), v
    return any(c < 0 for c in v)

def main():
    elements = [IDENTITY]
    index = {IDENTITY: 0}
    words = [()]
    length = [0]
    action = []
    for k, a in enumerate(elements):
        row = []
        for s in range(N):
            b = right(a, s)
            if b not in index:
                index[b] = len(elements)
                elements.append(b)
                words.append(words[k] + (s,))
                length.append(length[k] + 1)
            row.append(index[b])
        action.append(tuple(row))
    assert len(elements) == 51840
    inv = []
    for word in words:
        k = 0
        for s in reversed(word):
            k = action[k][s]
        inv.append(k)
    desc = [set(s for s in range(N) if negative(a[s])) for a in elements]
    left_desc = [desc[inv[k]] for k in range(len(elements))]
    left = [tuple(inv[action[inv[k]][s]] for s in range(N)) for k in range(len(elements))]
    fc = []
    for k in range(len(elements)):
        # A reduced braid suffix is detected after stripping descents. Matsumoto
        # gives the equivalence with full commutativity in simply laced type.
        fc.append(all(fc[action[k][s]] for s in desc[k]) and
                  not any(s in desc[k] and t in desc[k] for s, t in EDGES))
    assert sum(fc) == 662
    def leq(x, w):
        while x != w:
            if length[x] >= length[w]:
                return False
            s = min(desc[w])
            if s in desc[x]:
                x = action[x][s]
            w = action[w][s]
        return True
    terminal = []
    star_checks = 0
    commuting_terminals = 0
    for k in range(len(elements)):
        weak_bad = True
        for table, ds in ((action, desc), (left, left_desc)):
            for s in ds[k]:
                if any(t in ds[table[k][s]] for t in ADJ[s]):
                    weak_bad = False
            if fc[k]:
                # Check the complete length-three star domain, including stars
                # that raise length. There is exactly one domain-preserving move.
                for s, t in EDGES:
                    if len(ds[k] & {s, t}) == 1:
                        moves = [table[k][u] for u in (s, t)
                                 if len(ds[table[k][u]] & {s, t}) == 1]
                        assert len(moves) == 1
                        assert fc[moves[0]]
                        star_checks += 1
        if weak_bad:
            if fc[k]:
                commuting_terminals += 1
                assert len(words[k]) == len(set(words[k]))
                assert not any(s in words[k] and t in words[k] for s, t in EDGES)
            else:
                compatible = [x for x in range(len(elements)) if fc[x]
                              and desc[k] <= desc[x] and left_desc[k] <= left_desc[x]
                              and leq(x, k)]
                terminal.append({"word": words[k], "length": length[k],
                                 "right_descents": sorted(desc[k]),
                                 "left_descents": sorted(left_desc[k]),
                                 "candidates": [{"word": words[x], "length": length[x],
                                                 "interval_rank": length[k]-length[x]}
                                                for x in compatible]})
    assert len(terminal) == 1
    assert all(c["interval_rank"] % 2 == 0 for w in terminal for c in w["candidates"])
    result = {"type": "E6", "edges": EDGES, "element_count": len(elements),
              "fc_count": sum(fc), "length_distribution": dict(sorted(Counter(length).items())),
              "fc_star_checks": star_checks, "commuting_terminals": commuting_terminals,
              "non_fc_weak_bad_terminals": terminal,
              "status": "Independent finite terminal enumeration; all compatible bad ranks even.",
              "code_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    output = Path(__file__).resolve().parents[1] / "results/e6-independent-certificate.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
