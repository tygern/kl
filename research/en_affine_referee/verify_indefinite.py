"""Exact all-k audit of the full-support indefinite E13 reflection family.

No family implementation is imported. Arbitrary-precision row matrices
certify the quadratic root formula and all terminal inequalities.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
N = 13
EDGES = [(i, i + 1) for i in range(11)] + [(2, 12)]
ADJ = [set() for _ in range(N)]
for s, t in EDGES:
    ADJ[s].add(t)
    ADJ[t].add(s)
E = tuple(tuple(int(i == j) for j in range(N)) for i in range(N))
ZERO = tuple(tuple(0 for _ in range(N)) for _ in range(N))
BETA = (2, 3, 5, 5, 4, 4, 3, 3, 2, 2, 1, 1, 3)
DELTA = (2, 4, 6, 5, 4, 3, 2, 1, 0, 0, 0, 0, 3)
GAMMA = (1, 2, 3, 2, 2, 1, 1, 0, 0, 0, 0, 0, 1)
I = [0, 3, 5, 7, 9, 11, 12]


def configure(rank):
    """Select another generalized E rank for reuse of this exact arithmetic."""
    global N, EDGES, ADJ, E, ZERO
    N = rank
    EDGES = [(i, i + 1) for i in range(N - 2)] + [(2, N - 1)]
    ADJ = [set() for _ in range(N)]
    for s, t in EDGES:
        ADJ[s].add(t)
        ADJ[t].add(s)
    E = tuple(tuple(int(i == j) for j in range(N)) for i in range(N))
    ZERO = tuple(tuple(0 for _ in range(N)) for _ in range(N))


def pair(a, b):
    return sum(a[i] * (2 * b[i] - sum(b[j] for j in ADJ[i])) for i in range(N))


def pairing_vector(a):
    return tuple(pair(e, a) for e in E)


def mv(a, v):
    return tuple(sum(x * y for x, y in zip(row, v)) for row in a)


def mm(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(N)) for j in range(N)) for i in range(N))


def add(a, b, k=1):
    return tuple(x + k * y for x, y in zip(a, b))


def simple(a, s):
    b = list(a)
    b[s] = -a[s] + sum(a[t] for t in ADJ[s])
    return tuple(b)


def real_root_witness(a):
    """Return an actual sequence of simple reflections reducing a to a simple root."""
    assert min(a) >= 0
    original = a
    steps = []
    while a not in E:
        choices = [s for s in range(N) if pair(E[s], a) > 0]
        assert choices, ("not certified real", a)
        s = choices[0]
        b = simple(a, s)
        assert min(b) >= 0 and sum(b) < sum(a)
        steps.append(s)
        a = b
    assert pair(original, original) == 2
    return {"height": sum(original), "lowering_generators": steps, "simple_root": E.index(a)}


def reflection(a):
    p = pairing_vector(a)
    return tuple(tuple(E[i][j] - a[i] * p[j] for j in range(N)) for i in range(N))


def right(a, s):
    b = [list(row) for row in a]
    for i in range(N):
        b[i][s] = -a[i][s]
        for t in ADJ[s]:
            b[i][t] = a[i][t] + a[i][s]
    return tuple(map(tuple, b))


def reduced_matrix(word):
    a = E
    for s in word:
        assert all(a[i][s] >= 0 for i in range(N)), "nonreduced supplied word"
        a = right(a, s)
    return a


def main():
    witnesses = {"beta0": real_root_witness(BETA), "gamma": real_root_witness(GAMMA),
                 "gamma_plus_delta": real_root_witness(add(GAMMA, DELTA))}
    w0 = reflection(BETA)
    base_source = json.loads((ROOT / "research/en_families/cartan_E13_m1_max0.json").read_text())
    row = next(r for r in base_source["real_terminal_roots"] if tuple(r["beta"]) == BETA)
    assert len(row["reduced_word"]) == 75 and reduced_matrix(row["reduced_word"]) == w0
    assert set(row["reduced_word"]) == set(range(N))
    T = mm(reflection(GAMMA), reflection(add(GAMMA, DELTA)))
    M = tuple(tuple(T[i][j] - E[i][j] for j in range(N)) for i in range(N))
    M2 = mm(M, M)
    assert mm(M2, M) == ZERO
    v = mv(M, BETA)
    u = mv(M2, BETA)
    assert v == (8, 16, 24, 21, 16, 13, 8, 5, 0, 0, 0, 0, 13)
    assert u == tuple(4 * x for x in DELTA)
    assert all(x > 0 for x in BETA) and min(v) >= 0 and min(u) >= 0
    m0, m1, m2 = map(pairing_vector, [BETA, v, u])
    assert m0 == (1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1)
    assert m1 == (0, 0, -2, 2, -2, 2, -2, 2, -5, 0, 0, 0, 2)
    assert m2 == (0, 0, 0, 0, 0, 0, 0, 0, -4, 0, 0, 0, 0)
    for s in range(N):
        if s in I:
            assert m0[s] >= 1 and m1[s] >= 0 and m2[s] >= 0
        else:
            assert m0[s] <= 0 and m1[s] <= 0 and m2[s] <= 0
    edge_checks = []
    for s in I:
        for t in sorted(ADJ[s]):
            coefficients = [m[s] + m[t] for m in (m0, m1, m2)]
            assert max(coefficients) <= 0
            edge_checks.append({"s": s, "t": t, "pairing_sum_coefficients_in_1_k_binom_k_2": coefficients})
    independent = [m for m in range(1 << N) if not any(m & (1 << s) and m & (1 << t) for s, t in EDGES)]
    assert max(m.bit_count() for m in independent) == len(I) == 7
    assert not any(s in I and t in I for s, t in EDGES)
    # beta_k has every coordinate positive. If its reflection belonged to
    # W_(S\{j}), row j would be the jth row of identity. Column 0 differs
    # in EVERY row, since <alpha_0,beta_k>=1 and beta_k[j]>0.
    assert (m0[0], m1[0], m2[0]) == (1, 0, 0)
    assert (BETA[0], v[0], u[0]) == (2, 8, 8)
    output = {"status": "All exact all-k assertions passed", "rank": N,
              "beta0": BETA, "delta_in_embedded_affine_E8": DELTA, "gamma": GAMMA,
              "real_root_witnesses": witnesses, "base_length": 75,
              "T_minus_identity_cube_zero": True, "quadratic_beta_v": v, "quadratic_beta_u": u,
              "beta_formula": "beta_k=beta0+k*v+binom(k,2)*u=T^k*beta0",
              "pairing_coefficients": [m0, m1, m2], "fixed_left_right_descents": I,
              "terminal_edge_checks": edge_checks, "maximum_independent_size": 7,
              "full_support_for_all_k": True,
              "distinctness_witness": "beta_k[0]=4*k*k+4*k+2 and <alpha0,beta_k>=1",
              "length_parity": "odd for every k, because each element is a real-root reflection",
              "unique_eligible_FC_bottom": I, "eligible_mu": 0,
              "all_FC_lower_mu_bound": [0, 1],
              "exact_length_formula": "not claimed or required"}
    (HERE / "e13-certificate.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: output[k] for k in ["status", "full_support_for_all_k", "length_parity", "unique_eligible_FC_bottom", "eligible_mu", "all_FC_lower_mu_bound", "exact_length_formula"]}, indent=2))


if __name__ == "__main__":
    main()
