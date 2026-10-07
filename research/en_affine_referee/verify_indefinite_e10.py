"""Independent all-k E10 reflection-family audit, using referee arithmetic."""
import json
import verify_indefinite as q

q.configure(10)
BETA = (3, 7, 10, 9, 7, 6, 4, 3, 1, 6)
DELTA = (2, 4, 6, 5, 4, 3, 2, 1, 0, 3)
GAMMA = (1, 2, 3, 2, 2, 1, 1, 0, 0, 1)
I = [1, 3, 5, 7, 9]


def main():
    witnesses = {"beta0": q.real_root_witness(BETA), "gamma": q.real_root_witness(GAMMA),
                 "gamma_plus_delta": q.real_root_witness(q.add(GAMMA, DELTA))}
    source = json.loads((q.ROOT / "research/en_families/cartan_E10_m2_max1.json").read_text())
    row = next(r for r in source["real_terminal_roots"] if tuple(r["beta"]) == BETA)
    assert len(row["reduced_word"]) == 101
    assert q.reduced_matrix(row["reduced_word"]) == q.reflection(BETA)
    assert set(row["reduced_word"]) == set(range(q.N))
    T = q.mm(q.reflection(GAMMA), q.reflection(q.add(GAMMA, DELTA)))
    M = tuple(tuple(T[i][j] - q.E[i][j] for j in range(q.N)) for i in range(q.N))
    M2 = q.mm(M, M)
    assert q.mm(M2, M) == q.ZERO
    v = q.mv(M, BETA)
    u = q.mv(M2, BETA)
    assert v == (7, 14, 21, 18, 14, 11, 7, 4, 0, 11)
    assert u == tuple(2 * x for x in DELTA)
    assert q.mv(M, u) == (0,) * q.N
    assert min(BETA) > 0 and min(v) >= 0 and min(u) >= 0
    m0, m1, m2 = map(q.pairing_vector, [BETA, v, u])
    assert m0 == (-1, 1, -2, 1, -1, 1, -1, 1, -1, 2)
    assert m1 == (0, 0, -1, 1, -1, 1, -1, 1, -4, 1)
    assert m2 == (0, 0, 0, 0, 0, 0, 0, 0, -2, 0)
    for s in range(q.N):
        coefficients = (m0[s], m1[s], m2[s])
        if s in I:
            assert coefficients[0] >= 1 and min(coefficients[1:]) >= 0
        else:
            assert max(coefficients) <= 0
    edge_checks = []
    for s in I:
        for t in sorted(q.ADJ[s]):
            coefficients = [m[s] + m[t] for m in (m0, m1, m2)]
            assert max(coefficients) <= 0
            edge_checks.append({"s": s, "t": t, "pairing_sum_coefficients": coefficients})
    independent = [m for m in range(1 << q.N) if not any(m & (1 << s) and m & (1 << t) for s, t in q.EDGES)]
    assert max(m.bit_count() for m in independent) == len(I) == 5
    assert not any(s in I and t in I for s, t in q.EDGES)
    assert (m0[1], m1[1], m2[1]) == (1, 0, 0)
    assert (BETA[0], v[0], u[0]) == (3, 7, 4)
    output = {"status": "All exact all-k assertions passed", "rank": q.N,
              "beta0": BETA, "embedded_affine_delta": DELTA, "gamma": GAMMA,
              "real_root_witnesses": witnesses, "base_length": 101,
              "T_minus_identity_cube_zero": True, "N_beta0": v, "N_squared_beta0": u, "N_u_zero": True,
              "beta_formula": "beta_k=beta0+k*v+binom(k,2)*u=T^k*beta0",
              "pairing_coefficients": [m0, m1, m2], "fixed_left_right_descents": I,
              "terminal_edge_checks": edge_checks, "maximum_independent_size": 5,
              "full_support_for_all_k": True,
              "full_support_witness": "<alpha1,beta_k>=1, so column1 of r_beta_k differs from identity in every row",
              "distinctness_witness": "beta_k[0]=2*k*k+5*k+3",
              "length_parity": "odd for every k, because each element is a real-root reflection",
              "unique_eligible_FC_bottom": I, "eligible_mu": 0,
              "all_FC_lower_mu_bound": [0, 1], "exact_length_formula": "not claimed or required"}
    (q.HERE / "e10-certificate.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({k: output[k] for k in ["status", "rank", "full_support_for_all_k", "unique_eligible_FC_bottom", "eligible_mu", "all_FC_lower_mu_bound"]}, indent=2))


if __name__ == "__main__":
    main()
