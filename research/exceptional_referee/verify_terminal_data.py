"""Referee checks independent of either exceptional-group implementation.

Checks full Poincare polynomials, commuting terminal counts, terminal word
transcriptions, and the D6 identification. Uses only Python's standard library.
Run: python3 research/exceptional_referee/verify_terminal_data.py
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def poincare(degrees):
    result = [1]
    for d in degrees:
        next_result = [0] * (len(result) + d - 1)
        for i, coefficient in enumerate(result):
            for j in range(d):
                next_result[i + j] += coefficient
        result = next_result
    return result


def independent_sets(n, edges):
    return sum(not any(mask & (1 << s) and mask & (1 << t) for s, t in edges)
               for mask in range(1 << n))


def root_bounds(n, edges):
    adjacency = [set() for _ in range(n)]
    for s, t in edges:
        adjacency[s].add(t)
        adjacency[t].add(s)
    roots = {tuple(int(i == j) for i in range(n)) for j in range(n)}
    todo = list(roots)
    for root in todo:
        for s in range(n):
            image = list(root)
            image[s] = sum(root[t] for t in adjacency[s]) - root[s]
            image = tuple(image)
            if image not in roots:
                roots.add(image)
                todo.append(image)
    assert all(all(c >= 0 for c in root) or all(c <= 0 for c in root) for root in roots)
    return len(roots), max(abs(c) for root in roots for c in root)


def d6_product(word):
    # Gern's generators use one-based labels; right multiplication.
    w = list(range(1, 7))
    for s in word:
        if s == 1:
            w[0], w[1] = -w[1], -w[0]
        else:
            w[s - 2], w[s - 1] = w[s - 1], w[s - 2]
    return tuple(w)


def d6_length(w):
    return sum(w[i] > w[j] for i in range(6) for j in range(i + 1, 6)) + sum(
        -w[i] > w[j] for i in range(6) for j in range(i + 1, 6))


def main():
    certs = {n: json.loads((ROOT / f"results/e{n}-independent-certificate.json").read_text())
             for n in (6, 7)}
    for n, degrees in [(6, [2, 5, 6, 8, 9, 12]),
                       (7, [2, 6, 8, 10, 12, 14, 18])]:
        hist = certs[n]["length_distribution"]
        if isinstance(hist, dict):
            hist = [hist[str(i)] for i in range(len(hist))]
        assert hist == poincare(degrees)
        edges = [(i, i + 1) for i in range(n - 2)] + [(2, n - 1)]
        assert independent_sets(n, edges) == {6: 22, 7: 36}[n]
        assert root_bounds(n, edges) == {6: (72, 3), 7: (126, 4)}[n]

    expected_words = ["1326213", "13256213", "132543621324356",
                      "1325436210321432543621324356"]
    assert ["".join(map(str, b["word"])) for b in certs[7]["bad"]] == expected_words
    assert [len(b["word"]) for b in certs[7]["bad"]] == [7, 8, 15, 28]
    assert [b["bottoms"][0]["rank"] for b in certs[7]["bad"]] == [4, 4, 11, 24]
    assert all(len(b["bottoms"]) == 1 for b in certs[7]["bad"])
    assert set(certs[7]["bad"][3]["word"]) == set(range(7))

    # In Gern, the D6 edges are 1--3, 2--3, 3--4--5--6.
    d_to_e = {1: 1, 2: 6, 3: 2, 4: 3, 5: 4, 6: 5}
    e_to_d = {e: d for d, e in d_to_e.items()}
    d_edges = [(1, 3), (2, 3), (3, 4), (4, 5), (5, 6)]
    mapped_edges = {tuple(sorted((d_to_e[s], d_to_e[t]))) for s, t in d_edges}
    assert mapped_edges == {(1, 2), (2, 3), (3, 4), (4, 5), (2, 6)}
    bad = certs[7]["bad"][2]
    top = d6_product([e_to_d[s] for s in bad["word"]])
    bottom = d6_product([e_to_d[s] for s in bad["bottoms"][0]["word"]])
    assert top == (-1, -6, 3, -4, 5, -2)
    assert bottom == (-1, -2, 4, 3, 6, 5)
    assert d6_length(top) == 15 and d6_length(bottom) == 4
    output = {"status": "All assertions passed",
              "poincare_polynomials": "E6 and E7 match coefficient by coefficient",
              "commuting_terminal_counts": [22, 36],
              "maximum_absolute_root_coefficients": [3, 4],
              "E7_terminal_ranks": [4, 4, 11, 24],
              "D6_map_Gern_to_E7": d_to_e,
              "D6_bottom": bottom, "D6_top": top}
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
