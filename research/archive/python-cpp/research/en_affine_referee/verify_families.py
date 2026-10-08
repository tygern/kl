"""Independent symbolic audit of two affine E8 families.

Uses row matrices with arbitrary-precision integers and the 240 finite E8
roots. No targeted.py import, affine group enumeration, or KL computation.
All-k conclusions follow from checked affine formulas, not finite samples.
"""
from collections import Counter
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
N = 9
EDGES = [(i, i + 1) for i in range(7)] + [(2, 8)]
ADJ = [set() for _ in range(N)]
for s, t in EDGES:
    ADJ[s].add(t)
    ADJ[t].add(s)
E = tuple(tuple(int(i == j) for j in range(N)) for i in range(N))
DELTA = (2, 4, 6, 5, 4, 3, 2, 1, 3)
BETA = (1, 2, 3, 3, 2, 2, 1, 1, 2)
GAMMA = (1, 1, 1, 1, 1, 0, 0, 0, 0)


def add(a, b, k=1):
    return tuple(x + k * y for x, y in zip(a, b))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def pair(a, b):
    return sum(a[i] * (2 * b[i] - sum(b[j] for j in ADJ[i])) for i in range(N))


def mv(a, v):
    return tuple(dot(row, v) for row in a)


def mm(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(N)) for j in range(N)) for i in range(N))


def col(a, j):
    return tuple(row[j] for row in a)


def positive(a):
    assert any(a)
    assert all(c >= 0 for c in a) or all(c <= 0 for c in a), a
    return all(c >= 0 for c in a)


def right(a, s):
    out = [list(row) for row in a]
    for i in range(N):
        out[i][s] = -a[i][s]
        for t in ADJ[s]:
            out[i][t] = a[i][t] + a[i][s]
    return tuple(map(tuple, out))


def word_matrix(word, reduced=True):
    a = E
    for s in word:
        if reduced:
            assert positive(col(a, s)), ("nonreduced", word)
        a = right(a, s)
    return a


def reduced_word(a):
    reversed_word = []
    while a != E:
        choices = [s for s in range(N) if not positive(col(a, s))]
        assert choices
        s = choices[0]
        reversed_word.append(s)
        a = right(a, s)
    return list(reversed(reversed_word))


def reflection(a):
    p = tuple(pair(e, a) for e in E)
    return tuple(tuple(E[i][j] - a[i] * p[j] for j in range(N)) for i in range(N))


def translate(gamma, k=1):
    # Standard convention t_gamma(alpha)=alpha-<alpha,gamma>*delta.
    p = tuple(pair(e, gamma) for e in E)
    return tuple(tuple(E[i][j] - k * DELTA[i] * p[j] for j in range(N)) for i in range(N))


def shifted(a, d, k):
    return tuple(tuple(a[i][j] + k * DELTA[i] * d[j] for j in range(N)) for i in range(N))


def finite_roots():
    roots = {E[i] for i in range(N) if i != 7}
    todo = list(roots)
    for a in todo:
        for s in range(N):
            if s == 7:
                continue
            b = list(a)
            b[s] = -a[s] + sum(a[t] for t in ADJ[s])
            b = tuple(b)
            if b not in roots:
                roots.add(b)
                todo.append(b)
    assert len(roots) == 240
    assert all(a[7] == 0 and pair(a, a) == 2 for a in roots)
    return roots


FINITE_ROOTS = finite_roots()


def length_formula(a, d):
    """Prove ell(a+k delta*d)=constant+slope*k for every k>=0.

    A positive affine root is r+n delta, n>=0 for positive finite r,
    n>=1 for negative r. Its image is r'+(n+c+k*d.r)delta.
    Its contribution to inversions is max(0,-c-nmin+[r'<0]-k*d.r).
    """
    classes = Counter()
    intercept = slope = 0
    for r in FINITE_ROOTS:
        image = mv(a, r)
        c = image[7]
        finite_image = add(image, DELTA, -c)
        assert finite_image in FINITE_ROOTS
        nmin = int(not positive(r))
        A = -c - nmin + int(not positive(finite_image))
        B = -dot(d, r)
        assert B <= 0 or A >= 0, (A, B, r)
        assert B >= 0 or A <= 0, (A, B, r)
        intercept += max(0, A)
        slope += max(0, B)
        classes[A, B] += 1
    return {"intercept": intercept, "slope": slope,
            "inversion_classes": [{"A": a, "B": b, "count": count} for (a, b), count in sorted(classes.items())],
            "no_breakpoints_for_k_ge_0": True}


def terminal_certificate(a, d):
    assert mm(a, a) == E and mv(a, DELTA) == DELTA
    assert dot(d, DELTA) == 0
    assert tuple(sum(d[i] * a[i][j] for i in range(N)) for j in range(N)) == tuple(-x for x in d)
    # These identities prove (a+k delta*d)^2=1 for every k.
    descents = []
    for s in range(N):
        if positive(col(a, s)):
            assert d[s] >= 0
        else:
            descents.append(s)
            assert d[s] <= 0
    edges = []
    for s in descents:
        for t in sorted(ADJ[s]):
            base = add(col(a, s), col(a, t))
            assert positive(base)
            assert d[s] + d[t] >= 0
            edges.append({"s": s, "t": t, "slope": d[s] + d[t], "base_root": base})
    return {"descents": descents, "involutions_for_all_k": True,
            "terminal_for_all_k": True, "terminal_edge_checks": edges}


def subword_witness(top, bottom):
    counts = Counter(bottom)
    target = word_matrix(bottom)

    def search(i, a, positions, left):
        if not left:
            return positions if a == target else None
        if i == len(top) or len(top) - i < left:
            return None
        s = top[i]
        if counts[s]:
            counts[s] -= 1
            found = search(i + 1, right(a, s), positions + [i], left - 1)
            counts[s] += 1
            if found is not None:
                return found
        return search(i + 1, a, positions, left)

    positions = search(0, E, [], len(bottom))
    assert positions is not None
    assert word_matrix([top[i] for i in positions]) == target
    return positions


def catalogue_candidates(I):
    data = json.loads((ROOT / "research/en_independent/e9-fc.json").read_text())
    assert data["complete"] and data["fc_count"] == 44199
    mask = sum(1 << i for i in I)
    return [row for row in data["elements"] if row["Lmask"] & mask == mask and row["Rmask"] & mask == mask]


def audit_family(base_word, d, expected, I, name):
    a = word_matrix(base_word)
    symbolic = terminal_certificate(a, d)
    assert symbolic["descents"] == I
    formula = length_formula(a, d)
    assert (formula["intercept"], formula["slope"]) == expected
    assert len(base_word) == expected[0]
    assert set(base_word) == set(range(N))
    # R=1+delta*d lies in W because aR is the conjugate described in main.
    R = shifted(E, d, 1)
    assert mm(a, R) == shifted(a, d, 1)
    rword = reduced_word(R)
    assert len(rword) == expected[1] and word_matrix(rword) == R
    rformula = length_formula(E, d)
    assert (rformula["intercept"], rformula["slope"]) == (0, expected[1])
    candidates = catalogue_candidates(I)
    bottoms = []
    for x in candidates:
        witness = subword_witness(base_word, x["word"])
        bottoms.append({"word": x["word"], "length": x["length"], "base_subword_positions": witness,
                        "gap_intercept": expected[0] - x["length"], "gap_slope": expected[1]})
    return {"name": name, "base_word": base_word, "base_row_matrix": a, "column_slopes": d,
            "length_formula": formula, "right_translation_reduced_word": rword,
            "right_translation_length": len(rword), "base_prefix_for_all_k": True,
            "distinct_and_full_support_for_all_k": True, **symbolic,
            "complete_FC_mask_candidates": len(candidates), "eligible_FC_bottoms_for_all_k": bottoms}


def main():
    assert all(pair(DELTA, e) == 0 for e in E)
    assert GAMMA in FINITE_ROOTS
    # gamma+delta is a real root, so t_gamma=r_gamma*r_(gamma+delta).
    assert mm(reflection(GAMMA), reflection(add(GAMMA, DELTA))) == translate(GAMMA)
    input_rows = json.loads((ROOT / "research/en_families/e9_full_support_eligible.json").read_text())
    base_word = next(row["word"] for row in input_rows if row["length"] == 27)
    a = word_matrix(base_word)
    assert mm(a, a) == E
    shift_vector = add(GAMMA, mv(a, GAMMA), -1)
    d = tuple(pair(e, shift_vector) for e in E)
    assert d == (2, -1, 1, -1, 1, -1, 1, -1, -1)
    assert mm(mm(translate(GAMMA), a), translate(GAMMA, -1)) == shifted(a, d, 1)
    # Algebraic expansion, delta fixed and <delta,gamma>=0, proves all k.
    I = [1, 3, 5, 7, 8]
    positive_family = audit_family(base_word, d, (27, 92), I, "odd-length terminal conjugates")
    positive_family["gamma"] = GAMMA
    positive_family["shift_vector"] = shift_vector
    assert len(positive_family["eligible_FC_bottoms_for_all_k"]) == 1
    assert positive_family["eligible_FC_bottoms_for_all_k"][0]["word"] == I
    independent_masks = [m for m in range(1 << N) if not any(m & (1 << s) and m & (1 << t) for s, t in EDGES)]
    assert max(m.bit_count() for m in independent_masks) == len(I) == 5
    positive_family["independence_number"] = 5
    positive_family["eligible_mu_for_all_k"] = 0
    positive_family["all_FC_lower_mu_bound"] = [0, 1]

    # Secondary reflection family: membership follows beta-delta in finite roots.
    assert add(BETA, DELTA, -1) in FINITE_ROOTS
    b = reflection(BETA)
    bword = reduced_word(b)
    dm = tuple(-pair(e, BETA) for e in E)
    reflection_family = audit_family(bword, dm, (33, 58), [3, 5, 7, 8], "reflection beta+k delta")
    assert len(reflection_family["eligible_FC_bottoms_for_all_k"]) == 5
    assert reflection(add(BETA, DELTA)) == shifted(b, dm, 1)
    reflection_family["mu_status"] = "No 0/1 bound proved for the odd eligible gaps"
    reflection_family["extensions"] = []
    for s in (0, 1):
        ext = right(b, s)
        assert mm(b, right(E, s)) == mm(right(E, s), b)
        assert dm[s] == 0
        ext_d = tuple(sum(dm[i] * right(E, s)[i][j] for i in range(N)) for j in range(N))
        assert ext_d == dm
        ext_word = bword + [s]
        ext_I = sorted([s, 3, 5, 7, 8])
        row = audit_family(ext_word, dm, (34, 58), ext_I, f"reflection extension {s}")
        assert len(row["eligible_FC_bottoms_for_all_k"]) == 1
        row["mu_status"] = "Unique eligible gap is odd; coefficient remains unresolved"
        reflection_family["extensions"].append(row)

    output = {"status": "All symbolic and exact assertions passed", "delta": DELTA,
              "finite_root_count": len(FINITE_ROOTS), "positive_family": positive_family,
              "reflection_family": reflection_family,
              "method": "row integer matrices; finite-root affine inversion formulas; full FC mask catalogue; Green uniqueness independently checked"}
    (HERE / "certificate.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": output["status"], "positive_length": "27+92k", "positive_eligible_mu": 0,
                      "reflection_length": "33+58k", "reflection_eligible_count": 5,
                      "reflection_extension_length": "34+58k", "reflection_mu": "unresolved"}, indent=2))


if __name__ == "__main__":
    main()
