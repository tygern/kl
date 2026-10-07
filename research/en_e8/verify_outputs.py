"""Cross-check all terminal outputs using unbounded integer matrices.

This does not enumerate a Coxeter group. It verifies reduced word lengths,
descent masks, both terminal conditions, and equivalence of outputs obtained
from flat E7 cosets, recursive cosets, and the independent D7 implementation.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def check_certificate(path, n):
    data = json.loads((ROOT / path).read_text())
    edges = [(i, i + 1) for i in range(n - 2)] + [(2, n - 1)]
    adj = [set() for _ in range(n)]
    for s, t in edges:
        adj[s].add(t)
        adj[t].add(s)
    identity = tuple(tuple(int(i == j) for i in range(n)) for j in range(n))

    def positive(v):
        assert any(v)
        assert all(c >= 0 for c in v) or all(c <= 0 for c in v)
        return all(c >= 0 for c in v)

    def element(word):
        a = identity
        for s in word:
            assert positive(a[s]), (path, "nonreduced word", word)
            b = list(a)
            b[s] = tuple(-c for c in a[s])
            for t in adj[s]:
                b[t] = tuple(c + d for c, d in zip(a[t], a[s]))
            a = tuple(b)
        return a

    def desc(a):
        return sum(1 << s for s in range(n) if not positive(a[s]))

    def weak(a):
        for s in range(n):
            if not positive(a[s]):
                for t in adj[s]:
                    assert positive(tuple(c + d for c, d in zip(a[s], a[t])))

    result = {}
    bads = data.get("bad", data.get("non_fc_weak_bad_terminals"))
    for b in bads:
        w = element(b["word"])
        wi = element(list(reversed(b["word"])))
        assert len(b["word"]) == b["length"]
        weak(w)
        weak(wi)
        rd = b.get("Rmask", sum(1 << s for s in b.get("right_descents", [])))
        ld = b.get("Lmask", sum(1 << s for s in b.get("left_descents", [])))
        assert rd == desc(w) and ld == desc(wi)
        bottoms = set()
        for x in b.get("bottoms", b.get("candidates", [])):
            a = element(x["word"])
            ai = element(list(reversed(x["word"])))
            assert len(x["word"]) == x["length"]
            assert rd & desc(a) == rd and ld & desc(ai) == ld
            support = set(x["word"])
            assert len(support) == len(x["word"])
            assert not any(s in support and t in support for s, t in edges)
            # A product of distinct commuting generators is below w whenever
            # its support is contained in w: select one of each as a subword.
            assert support <= set(b["word"])
            gap = b["length"] - x["length"]
            assert gap == x.get("rank", x.get("interval_rank"))
            bottoms.add((a, x["length"], gap))
        assert w not in result
        result[w] = (b["length"], rd, ld, frozenset(bottoms))
    return data, result


def convolve(a, b):
    result = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i + j] += x * y
    return result


def poincare(degrees):
    result = [1]
    for d in degrees:
        result = convolve(result, [1] * d)
    return result


def main():
    for n in (6, 7):
        _, reference = check_certificate(f"results/e{n}-independent-certificate.json", n)
        _, flat = check_certificate(f"research/en_e8/e{n}-validation.json", n)
        _, recursive = check_certificate(f"research/en_e8/e{n}-recursive.json", n)
        assert reference == flat == recursive

    flat_data, flat = check_certificate("research/en_e8/e8-terminals.json", 8)
    recursive_data, recursive = check_certificate("research/en_e8/e8-recursive.json", 8)
    d7_data, d7 = check_certificate("research/en_independent/e8-d7-terminals.json", 8)
    assert flat == recursive == d7
    assert flat_data["right_terminal_count"] == recursive_data["right_terminal_count"] == d7_data["full_right_terminals"] == 2160
    assert flat_data["commuting_terminals"] == recursive_data["commuting_terminals"] == d7_data["commuting_terminals"] == 58
    assert len(flat) == 6
    assert d7_data["parabolic_histogram"] == poincare([2, 4, 6, 7, 8, 10, 12])
    assert convolve(d7_data["parabolic_histogram"], d7_data["coset_histogram"]) == poincare([2, 8, 12, 14, 18, 20, 24, 30])
    print(json.dumps({
        "status": "All assertions passed",
        "E6_E7": "Flat and recursive methods match independent complete enumerations",
        "E8": "Flat E7, recursive parabolic, and independent D7 methods agree as integer matrices",
        "right_terminal_count": 2160,
        "commuting_terminal_count": 58,
        "noncommuting_terminal_count": 6,
        "eligible_gaps": sorted(next(iter(entry[3]))[2] for entry in flat.values()),
        "D7_and_E8_Poincare_histograms": "verified coefficient by coefficient",
    }, indent=2))


if __name__ == "__main__":
    main()
