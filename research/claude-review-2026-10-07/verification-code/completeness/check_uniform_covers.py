"""Uniform family b_{r,k} = r_{beta_{r,k}} in E_{4r+1}: build the reflection matrix,
check length, descents, terminality, then enumerate all Bruhat covers (single-letter
deletions of a reduced word, keep length l-1) and test each for full commutativity.
If no cover is FC, Theorem uniform's 'eventual vanishing' already holds at that k."""
import sys
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/completeness')
from en import En

def beta_rk(r, k):
    n = 4 * r + 1
    beta = [0] * n
    beta[0], beta[1], beta[2], beta[4 * r] = r - 1, r, 2 * r - 1, r
    for j in range(3, 4 * r):
        beta[j] = 2 * r - j // 2
    delta = [0] * n; gamma = [0] * n
    for j, v in enumerate((2, 4, 6, 5, 4, 3, 2, 1)):
        delta[j] = v
    delta[4 * r] = 3
    for j, v in enumerate((1, 2, 3, 2, 2, 1, 1, 0)):
        gamma[j] = v
    gamma[4 * r] = 1
    a = 2 * r - 4
    return tuple(beta[j] - a * k * gamma[j] + (a * k * k + r * k) * delta[j] for j in range(n))

def reflection(G, beta):
    n = G.n
    cols = []
    for j in range(n):
        m = G.pair(beta, j)
        col = tuple((1 if i == j else 0) - m * beta[i] for i in range(n))
        cols.append(col)
    return tuple(cols)

def covers(G, w):
    word = G.word(w)
    l = len(word)
    seen = {}
    for p in range(l):
        sub = word[:p] + word[p + 1:]
        v = G.from_word(sub)
        if v in seen:
            continue
        lv = G.length(v)
        if lv == l - 1:
            seen[v] = sub
    return word, seen

if __name__ == '__main__':
    for r, k in [(3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (5, 0)]:
        n = 4 * r + 1
        G = En(n)
        beta = beta_rk(r, k)
        assert G.pair(beta, 0) == r - 2 and sum(beta[i] * G.pair(beta, i) for i in range(n)) == 2
        b = reflection(G, beta)
        word = G.word(b)
        L = G.L(b); R = G.R(b)
        Ir = sorted({0, 4 * r} | set(range(3, 4 * r, 2)))
        print(f"r={r} k={k} E{n}: len={len(word)} (8r^2+3+116k={8*r*r+3+116*k}) R={R} L==R={L==R} R==I_r={R==Ir} "
              f"terminal={G.terminal(b)} fullsupp={len(set(word))==n} FC={G.is_fc_word(word)}")
        word, cov = covers(G, b)
        fc_cov = [v for v, sub in cov.items() if G.is_fc_word(G.word(v))]
        print(f"   distinct Bruhat covers: {len(cov)}; FC covers: {len(fc_cov)}")
        sys.stdout.flush()
