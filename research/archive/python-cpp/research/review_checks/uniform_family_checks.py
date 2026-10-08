"""Finite checks of the ingredients of the manuscript's Lemma on covers of
terminal reflections, for the rank-uniform family b_{r,k} = r_{beta_{r,k}}
in W(E_{4r+1}), and the observed length formula l(b_{r,k}) = 8r^2 + 3 + 116k.

For each (r,k) checked:
  * beta_{r,k} has norm 2, b = r_beta is a terminal involution of odd length
    8r^2 + 3 + 116k with L(b) = R(b) = I_r, full support, not fully
    commutative (Stembridge's heap criterion on a reduced word);
  * for every s in I_r the element b s has L(b s) = I_r, R(b s) = I_r minus {s},
    length l(b) - 1, and is not fully commutative; s b is its inverse;
  * s b s = r_{s beta} with s beta a positive non-simple root, and
    (alpha_s, s beta) = -(alpha_s, beta) < 0 (the pairing used in Lemma A);
  * for selected (r,k): every Bruhat cover of b (one letter deleted from a
    reduced word, keeping length l(b) - 1; exhaustive by the strong exchange
    condition) is not fully commutative.

Nothing is imported from the repository.  Elements are integer matrices
stored as tuples of columns (column j = w(alpha_j) in the simple-root basis);
exact integer arithmetic only.  E_n labelling: chain 0-1-...-(n-2), node n-1
attached to node 2.  A string of labels denotes the product of the simple
reflections in the order written.

Run: python3 research/review_checks/uniform_family_checks.py [--full]
Default: r = 3..8, k = 0..3 (descent, length and non-FC checks) and covers for
(r,k) in {(3,0), (3,1), (4,0)}.  --full adds r = 9, 10, 12 and covers for
(3,2), (4,1), (5,0), (6,0), (7,0).  Writes results/uniform-family-certificate.json.
Python 3.10+, standard library only.
"""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]


class En:
    def __init__(self, n):
        self.n = n
        adj = [[] for _ in range(n)]
        for i in range(n - 2):
            adj[i].append(i + 1)
            adj[i + 1].append(i)
        adj[n - 1].append(2)
        adj[2].append(n - 1)
        self.adj = adj
        self.adjset = [set(a) for a in adj]
        self.e = tuple(tuple(1 if i == j else 0 for i in range(n)) for j in range(n))

    def pair(self, v, i):
        """(v, alpha_i) for the symmetric bilinear form with (alpha_i, alpha_i) = 2."""
        return 2 * v[i] - sum(v[j] for j in self.adj[i])

    def form(self, v, w):
        return sum(v[i] * self.pair(w, i) for i in range(self.n))

    def refl_vec(self, v, i):
        c = self.pair(v, i)
        if c == 0:
            return v
        v = list(v)
        v[i] -= c
        return tuple(v)

    def rmul(self, w, i):
        cols = list(w)
        ai = cols[i]
        for j in self.adj[i]:
            cols[j] = tuple(x + y for x, y in zip(cols[j], ai))
        cols[i] = tuple(-x for x in ai)
        return tuple(cols)

    def lmul(self, i, w):
        return tuple(self.refl_vec(col, i) for col in w)

    def mul(self, u, v):
        n = self.n
        out = []
        for j in range(n):
            col = [0] * n
            for i in range(n):
                c = v[j][i]
                if c:
                    ui = u[i]
                    for k in range(n):
                        col[k] += c * ui[k]
            out.append(tuple(col))
        return tuple(out)

    @staticmethod
    def neg(col):
        return any(x < 0 for x in col) and all(x <= 0 for x in col)

    @staticmethod
    def pos(col):
        return any(x > 0 for x in col) and all(x >= 0 for x in col)

    def R(self, w):
        return [i for i in range(self.n) if self.neg(w[i])]

    def word(self, w):
        """A reduced word, by stripping the smallest right descent repeatedly."""
        out = []
        while True:
            d = next((i for i in range(self.n) if self.neg(w[i])), -1)
            if d < 0:
                break
            out.append(d)
            w = self.rmul(w, d)
        out.reverse()
        return out

    def length(self, w):
        return len(self.word(w))

    def inverse(self, w):
        v = self.e
        for i in reversed(self.word(w)):
            v = self.rmul(v, i)
        return v

    def L(self, w):
        return self.R(self.inverse(w))

    def from_word(self, word):
        w = self.e
        for i in word:
            w = self.rmul(w, i)
        return w

    def reflection(self, beta):
        n = self.n
        cols = []
        for j in range(n):
            m = self.pair(beta, j)
            cols.append(tuple((1 if i == j else 0) - m * beta[i] for i in range(n)))
        return tuple(cols)

    def right_terminal(self, w):
        for s in self.R(w):
            for t in self.adj[s]:
                if not self.pos(tuple(a + b for a, b in zip(w[s], w[t]))):
                    return False
        return True

    def terminal(self, w):
        return self.right_terminal(w) and self.right_terminal(self.inverse(w))

    def is_fc_word(self, word):
        """Stembridge's heap criterion for a reduced word in a simply laced group: fully
        commutative iff between two consecutive occurrences of a generator s there are at
        least two occurrences of neighbours of s."""
        last = {}
        for p, s in enumerate(word):
            if s in last:
                if sum(1 for t in word[last[s] + 1:p] if t in self.adjset[s]) < 2:
                    return False
            last[s] = p
        return True

    def covers(self, w):
        """All Bruhat covers of w: delete one letter of a reduced word, keep length l(w)-1."""
        word = self.word(w)
        seen = {}
        for p in range(len(word)):
            sub = word[:p] + word[p + 1:]
            v = self.from_word(sub)
            if v not in seen and self.length(v) == len(word) - 1:
                seen[v] = sub
        return seen


def uniform_data(r):
    """beta_r, gamma, delta of the manuscript's Section 4.1 in E_{4r+1}; a = 2r - 4."""
    n = 4 * r + 1
    beta = [0] * n
    beta[0], beta[1], beta[2], beta[4 * r] = r - 1, r, 2 * r - 1, r
    for j in range(3, 4 * r):
        beta[j] = 2 * r - j // 2
    delta = [0] * n
    gamma = [0] * n
    for j, v in enumerate((2, 4, 6, 5, 4, 3, 2, 1)):
        delta[j] = v
    delta[4 * r] = 3
    for j, v in enumerate((1, 2, 3, 2, 2, 1, 1, 0)):
        gamma[j] = v
    gamma[4 * r] = 1
    return n, tuple(beta), tuple(gamma), tuple(delta), 2 * r - 4


def beta_rk(r, k):
    n, beta, gamma, delta, a = uniform_data(r)
    return tuple(beta[j] - a * k * gamma[j] + (a * k * k + r * k) * delta[j] for j in range(n))


def I_r(r):
    return sorted({0, 4 * r} | set(range(3, 4 * r, 2)))


failures = []


def expect(cond, what):
    if not cond:
        failures.append(what)
        print('EXPECTATION FAILED:', what, file=sys.stderr)


def check_pair(r, k):
    n = 4 * r + 1
    G = En(n)
    beta = beta_rk(r, k)
    I = I_r(r)
    expect(G.form(beta, beta) == 2, f'r={r} k={k}: norm of beta is 2')
    expect(all(c > 0 for c in beta), f'r={r} k={k}: beta is positive')
    b = G.reflection(beta)
    word = G.word(b)
    predicted = 8 * r * r + 3 + 116 * k
    expect(len(word) == predicted, f'r={r} k={k}: l(b) = 8r^2+3+116k')
    expect(len(word) % 2 == 1, f'r={r} k={k}: odd length')
    expect(G.R(b) == I and G.L(b) == I, f'r={r} k={k}: L(b) = R(b) = I_r')
    expect(len(I) == 2 * r + 1, f'r={r}: |I_r| = 2r+1')
    expect(not any(t in I for s in I for t in G.adj[s]), f'r={r}: I_r is independent')
    expect(G.terminal(b), f'r={r} k={k}: terminal')
    expect(len(set(word)) == n, f'r={r} k={k}: full support')
    expect(not G.is_fc_word(word), f'r={r} k={k}: b is not fully commutative')
    expect(G.mul(b, b) == G.e, f'r={r} k={k}: involution')
    bs_rows = []
    for s in I:
        x = G.rmul(b, s)
        xw = G.word(x)
        m = G.pair(beta, s)
        sbeta = G.refl_vec(beta, s)
        ok = (len(xw) == len(word) - 1 and G.L(x) == I and G.R(x) == sorted(set(I) - {s}) and not G.is_fc_word(xw))
        expect(ok, f'r={r} k={k} s={s}: bs has L=I_r, R=I_r-{{s}}, length l(b)-1 and is not FC')
        expect(G.inverse(x) == G.lmul(s, b), f'r={r} k={k} s={s}: (bs)^-1 = sb')
        expect(m >= 1, f'r={r} k={k} s={s}: (alpha_s, beta) >= 1 since s is a descent')
        expect(all(c >= 0 for c in sbeta) and sum(sbeta) > 1, f'r={r} k={k} s={s}: s beta is a positive non-simple root')
        expect(G.pair(sbeta, s) == -m, f'r={r} k={k} s={s}: (alpha_s, s beta) = -(alpha_s, beta)')
        expect(G.lmul(s, x) == G.reflection(sbeta), f'r={r} k={k} s={s}: s b s = r_(s beta)')
        bs_rows.append({'s': s, 'pairing_m': m, 'bs_length': len(xw), 'bs_L': G.L(x), 'bs_R': G.R(x), 'bs_fully_commutative': G.is_fc_word(xw)})
    return {'r': r, 'k': k, 'rank': n, 'beta': list(beta), 'I_r': I, 'length': len(word), 'length_formula_8r2_plus_3_plus_116k': predicted,
            'L_equals_R_equals_I_r': True, 'terminal': True, 'full_support': True, 'fully_commutative': False, 'involution': True,
            'bs_elements': bs_rows}


def check_covers(r, k):
    n = 4 * r + 1
    G = En(n)
    b = G.reflection(beta_rk(r, k))
    cov = G.covers(b)
    fc = 0
    min_support = n
    for v, sub in cov.items():
        w = G.word(v)
        min_support = min(min_support, len(set(w)))
        if G.is_fc_word(w):
            fc += 1
    expect(fc == 0, f'r={r} k={k}: no Bruhat cover of b is fully commutative')
    return {'r': r, 'k': k, 'rank': n, 'length': 8 * r * r + 3 + 116 * k, 'distinct_covers': len(cov), 'fully_commutative_covers': fc,
            'minimum_cover_support': min_support}


def main():
    full = '--full' in sys.argv[1:]
    pairs = [(r, k) for r in range(3, 9) for k in range(4)]
    cover_pairs = [(3, 0), (3, 1), (4, 0)]
    if full:
        pairs += [(9, 0), (9, 3), (10, 0), (10, 2), (12, 0), (12, 1)]
        cover_pairs += [(3, 2), (4, 1), (5, 0), (6, 0), (7, 0)]
    rows = [check_pair(r, k) for r, k in pairs]
    covers = [check_covers(r, k) for r, k in cover_pairs]
    expected_cover_counts = {(3, 0): 74, (3, 1): 162, (3, 2): 230, (4, 0): 130, (4, 1): 218, (5, 0): 202, (6, 0): 290, (7, 0): 394}
    for c in covers:
        key = (c['r'], c['k'])
        if key in expected_cover_counts:
            expect(c['distinct_covers'] == expected_cover_counts[key], f'r={key[0]} k={key[1]}: {expected_cover_counts[key]} distinct covers')
    report = {'family': 'b_{r,k} = r_{beta_{r,k}} in W(E_{4r+1}), beta_{r,k} = beta_r - a k gamma + (a k^2 + r k) delta, a = 2r-4',
              'conventions': 'E_n: chain 0-1-...-(n-2), node n-1 attached to node 2; words act on the right; columns of w are w(alpha_j)',
              'checked_pairs': [[r, k] for r, k in pairs], 'length_formula_holds_for_all_checked_pairs': all(r['length'] == r['length_formula_8r2_plus_3_plus_116k'] for r in rows),
              'bs_checks': 'for every s in I_r: L(bs) = I_r, R(bs) = I_r - {s}, l(bs) = l(b) - 1, bs not FC, sbs = r_{s beta} with s beta positive non-simple and (alpha_s, s beta) = -(alpha_s, beta) < 0',
              'pairs': rows, 'bruhat_covers': covers,
              'status': 'FAILED' if failures else 'passed', 'failed_expectations': failures}
    out = ROOT / 'results/uniform-family-certificate.json'
    out.write_text(json.dumps(report, indent=1) + '\n')
    print(json.dumps({'status': report['status'], 'pairs_checked': len(rows), 'lengths': {f'({r["r"]},{r["k"]})': r['length'] for r in rows},
                      'covers': {f'({c["r"]},{c["k"]})': [c['distinct_covers'], c['fully_commutative_covers']] for c in covers}}))
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
