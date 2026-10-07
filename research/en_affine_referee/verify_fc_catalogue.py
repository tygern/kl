"""Independently certify the complete E9 FC catalogue and base cover exclusion.

Uses arbitrary-precision row matrices from verify_families.py. Completeness
is checked by FC-ascent closure, not assumed from the input's status field.
"""
from collections import Counter
import json
from verify_families import ROOT, HERE, N, E, ADJ, EDGES, col, positive, right, word_matrix, reduced_word


def desc(a):
    return sum(1 << s for s in range(N) if not positive(col(a, s)))


def independent(mask):
    return not any(mask & (1 << s) and mask & (1 << t) for s, t in EDGES)


def left(a, s):
    out = list(a)
    out[s] = tuple(-a[s][j] + sum(a[t][j] for t in ADJ[s]) for j in range(N))
    return tuple(out)


def main():
    data = json.loads((ROOT / "research/en_independent/e9-fc.json").read_text())
    states = {}
    by_word = {}
    histogram = Counter()
    previous_length = -1
    for row in data["elements"]:
        word = tuple(row["word"])
        assert len(word) == row["length"] >= previous_length
        previous_length = row["length"]
        if not word:
            a = ai = E
        else:
            parent, inverse_parent = by_word[word[:-1]]
            assert positive(col(parent, word[-1]))
            a = right(parent, word[-1])
            ai = left(inverse_parent, word[-1])
        assert a not in states
        d = desc(a)
        assert d == row["Rmask"] and desc(ai) == row["Lmask"]
        assert independent(d)
        # Inductively certify each listed element FC, using all predecessors.
        for s in range(N):
            if d & (1 << s):
                assert states[right(a, s)] == len(word) - 1
        states[a] = len(word)
        by_word[word] = (a, ai)
        histogram[len(word)] += 1
    assert len(states) == 44199 and max(histogram) == 44
    ascents = missing_noncommuting = missing_bad_predecessor = 0
    for a, length in states.items():
        d = desc(a)
        for s in range(N):
            if d & (1 << s):
                continue
            ascents += 1
            b = right(a, s)
            if b in states:
                assert states[b] == length + 1
                continue
            bd = desc(b)
            if not independent(bd):
                missing_noncommuting += 1
                continue
            assert any(right(b, t) not in states for t in range(N) if bd & (1 << t)), "missing FC ascent"
            missing_bad_predecessor += 1
    # A shortest omitted FC element would have all shorter FC predecessors
    # present and would violate the last assertion. Therefore closure is complete.
    assert [histogram[i] for i in range(45)] == data["length_distribution"]
    assert ascents == data["ascents_tested"]
    inputs = json.loads((ROOT / "research/en_families/e9_full_support_eligible.json").read_text())
    word = next(row["word"] for row in inputs if row["length"] == 27)
    covers = {}
    for i in range(len(word)):
        a = word_matrix(word[:i] + word[i + 1:], reduced=False)
        rw = reduced_word(a)
        if len(rw) == 26:
            covers[a] = rw
    assert len(covers) == 21
    assert all(a not in states for a in covers)
    output = {"status": "All assertions passed", "FC_catalogue_independently_complete": True,
              "FC_count": len(states), "maximum_FC_length": 44,
              "FC_ascent_tests": ascents, "missing_with_noncommuting_descents": missing_noncommuting,
              "missing_with_nonFC_predecessor": missing_bad_predecessor,
              "base_Bruhat_cover_count": len(covers), "FC_base_covers": 0,
              "base_cover_reduced_words": sorted(covers.values()),
              "conclusion": "For every k>=0 and every fully commutative x, ordinary mu(x,b_k)=0 for the length 27+92k family"}
    (HERE / "fc-cover-certificate.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: v for k, v in output.items() if k != "base_cover_reduced_words"}, indent=2))


if __name__ == "__main__":
    main()
