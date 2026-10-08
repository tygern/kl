"""Checks on b = b_{r,k} = r_{beta_{r,k}} in E_{4r+1}:
  * b is a reflection with L(b)=R(b)=I_r, odd length; l(b) vs 2dp(beta)-1
  * for each s in I_r: x = bs has s in L(x), L(x)=I_r, R(x)=I_r\{s}, l(x)=l(b)-1, x not FC;
    sbs = r_{s beta}, s beta positive non-simple, l(sbs)=l(b)-2, sbs not FC
  * exhaustive covers (single deletions from a reduced word of b): every cover is bs, sb
    (s in I_r) or has L,R containing I_r; no cover is FC.
"""
import sys
from elib import *

pairs = [(3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (5, 0), (5, 1), (6, 0), (7, 0)]
if len(sys.argv) > 1:
    pairs = [tuple(map(int, p.split(','))) for p in sys.argv[1:]]

for (r, k) in pairs:
    n = 4 * r + 1
    A = cartan(n)
    beta = uniform_root(r, k)
    Ir = I_r(r)
    assert pairing(A, beta, beta) == 2, "beta not norm 2"
    assert all(x >= 1 for x in beta)
    b = reflection_matrix(A, beta)
    wb = reduced_word(A, b)
    lb = len(wb)
    assert word_to_matrix(A, wb) == b
    assert word_to_matrix(A, list(reversed(wb))) == b  # involution
    Rb = right_descents(b)
    dp = depth(A, beta)
    print(f"r={r} k={k} n={n}: l(b)={lb} (odd={lb % 2 == 1}), 2dp-1={2*dp-1}, R(b)=I_r: {Rb == Ir}, "
          f"8r^2+3+116k={8*r*r+3+116*k}, FC(b)={is_fc_word(A, wb)}")

    # --- bs for s in I_r
    for s in Ir:
        x = mult_right_simple(A, b, s)
        wx = reduced_word(A, x)
        xinv = word_to_matrix(A, list(reversed(wx)))
        Lx = right_descents(xinv)
        Rx = right_descents(x)
        sbeta = simple_refl(A, s, beta)
        m_s = sum(A[s][j] * beta[j] for j in range(n))
        sbs = mult_left_simple(A, s, x)
        assert sbs == reflection_matrix(A, sbeta), "sbs != r_{s beta}"
        lsbs = len(reduced_word(A, sbs))
        ok = (len(wx) == lb - 1 and s in Lx and Lx == Ir and Rx == [t for t in Ir if t != s]
              and is_positive(sbeta) and sum(sbeta) >= 2 and lsbs == lb - 2
              and not is_fc_word(A, wx) and not is_fc(A, sbs) and m_s >= 1)
        zeros = [j for j in range(n) if sbeta[j] == 0]
        print(f"   s={s:2d}: m_s={m_s}, l(bs)={len(wx)}, s in L(bs)={s in Lx}, L(bs)=I_r={Lx == Ir}, "
              f"R(bs)=I_r-s={Rx == [t for t in Ir if t != s]}, s.beta>0 nonsimple={is_positive(sbeta) and sum(sbeta) >= 2} "
              f"(zero coords {zeros}), l(sbs)={lsbs}, FC(bs)={is_fc_word(A, wx)}, FC(sbs)={is_fc(A, sbs)}  -> {'OK' if ok else 'FAIL'}")

    # --- exhaustive covers via single deletions
    covers = {}
    for i in range(lb):
        w = wb[:i] + wb[i + 1:]
        M = word_to_matrix(A, w)
        if M in covers:
            continue
        if len(reduced_word(A, M)) == lb - 1:
            covers[M] = w
    sb_set = {mult_left_simple(A, s, b) for s in Ir}
    bs_set = {mult_right_simple(A, b, s) for s in Ir}
    n_fc = 0
    n_other = 0
    bad_other = 0
    for M, w in covers.items():
        if is_fc_word(A, reduced_word(A, M)):
            n_fc += 1
        if M not in sb_set and M not in bs_set:
            n_other += 1
            Linv = word_to_matrix(A, list(reversed(reduced_word(A, M))))
            if not (set(Ir) <= set(right_descents(Linv)) and set(Ir) <= set(right_descents(M))):
                bad_other += 1
    print(f"   covers: {len(covers)} total, {len(sb_set | bs_set)} of form sb/bs, {n_other} others "
          f"(of which {bad_other} fail L,R ⊇ I_r), FC covers: {n_fc}")
    sys.stdout.flush()
