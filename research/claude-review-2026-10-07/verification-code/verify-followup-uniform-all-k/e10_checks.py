"""E10 appendix family: beta_k = beta_0 - k gamma + (k^2+3k) delta, b_k = r_{beta_k}.
Check no cover of b_k is FC (k=0,1) and the bs covers are non-FC via the reflection lemma."""
from elib import *
import sys

n = 10
A = cartan(n)
delta = (2, 4, 6, 5, 4, 3, 2, 1, 0, 3)
gamma = (1, 2, 3, 2, 2, 1, 1, 0, 0, 1)
beta0 = (3, 7, 10, 9, 7, 6, 4, 3, 1, 6)
I = [1, 3, 5, 7, 9]
for k in (0, 1):
    beta = tuple(beta0[i] - k * gamma[i] + (k * k + 3 * k) * delta[i] for i in range(n))
    assert pairing(A, beta, beta) == 2
    b = reflection_matrix(A, beta)
    wb = reduced_word(A, b)
    lb = len(wb)
    print(f"E10 k={k}: l(b)={lb}, 2dp-1={2*depth(A, beta)-1}, R(b)={right_descents(b)}, FC(b)={is_fc_word(A, wb)}")
    for s in I:
        x = mult_right_simple(A, b, s)
        sbeta = simple_refl(A, s, beta)
        print(f"   s={s}: l(bs)={len(reduced_word(A, x))}, s.beta nonsimple={is_positive(sbeta) and sum(sbeta) >= 2}, FC(bs)={is_fc(A, x)}")
    covers = {}
    for i in range(lb):
        w = wb[:i] + wb[i + 1:]
        M = word_to_matrix(A, w)
        if M not in covers and len(reduced_word(A, M)) == lb - 1:
            covers[M] = w
    n_fc = sum(1 for M in covers if is_fc(A, M))
    print(f"   covers: {len(covers)}, FC covers: {n_fc}")
    sys.stdout.flush()
