"""Independent check of finding followup-uniform-all-k-D1.

(1) exact length of b_{r,k} = r_{beta_{r,k}} in E_{4r+1} (descent stripping on the
    reflection matrix, exact integers), compared with 8r^2+3+116k and with 2*ht-1;
(2) the explicit word family of length 8r^2+7r-2 in E_{4r+1}: reduced + FC, using the
    parabolic-strip criterion (ws non-FC iff w^J has a reduced word ending in s t,
    t adjacent to s, J = commutant of s), validated against Catalan/D4/E6/E7 counts.
"""
import sys
from itertools import product

def cartan(n):
    # chain 0-1-...-(n-2), node n-1 attached to 2
    A = [[0]*n for _ in range(n)]
    for i in range(n): A[i][i] = 2
    def e(i, j): A[i][j] = A[j][i] = -1
    for i in range(n-2): e(i, i+1)
    e(2, n-1)
    return A

def pair(A, u, v):
    n = len(A)
    Av = [sum(A[i][j]*v[j] for j in range(n)) for i in range(n)]
    return sum(u[i]*Av[i] for i in range(n))

def beta_rk(r, k):
    n = 4*r+1
    b = [0]*n
    b[0], b[1], b[2], b[4*r] = r-1, r, 2*r-1, r
    for j in range(3, 4*r): b[j] = 2*r - j//2
    delta = [0]*n; gamma = [0]*n
    for i, v in enumerate((2,4,6,5,4,3,2,1)): delta[i] = v
    delta[4*r] = 3
    for i, v in enumerate((1,2,3,2,2,1,1,0)): gamma[i] = v
    gamma[4*r] = 1
    a = 2*r-4
    return [b[i] - a*k*gamma[i] + (a*k*k + r*k)*delta[i] for i in range(n)]

def refl_cols(A, beta):
    n = len(A)
    cols = []
    for j in range(n):
        p = sum(A[j][i]*beta[i] for i in range(n))  # (alpha_j, beta)
        cols.append([(1 if i == j else 0) - p*beta[i] for i in range(n)])
    return cols

def length_of_matrix(A, cols):
    """cols[j] = w(alpha_j) in root coordinates. Strip right descents."""
    n = len(A)
    cols = [c[:] for c in cols]
    L = 0
    while True:
        s = next((j for j in range(n) if all(x <= 0 for x in cols[j]) and any(x < 0 for x in cols[j])), None)
        if s is None:
            return L
        # right-multiply by s: (ws)(alpha_j) = w(alpha_j) - A[j][s] w(alpha_s)
        cs = cols[s]
        for j in range(n):
            if j != s and A[j][s] != 0:
                cols[j] = [cols[j][i] - A[j][s]*cs[i] for i in range(n)]
        cols[s] = [-x for x in cs]
        L += 1

# ---------- height-vector representation and FC test ----------
class Group:
    def __init__(self, A):
        self.A = A; self.n = len(A)
        self.nb = [[j for j in range(self.n) if j != i and A[i][j] != 0] for i in range(self.n)]
        self.J = [[j for j in range(self.n) if j != i and A[i][j] == 0] for i in range(self.n)]
    def mul(self, p, s):
        # p_i = ht(w(alpha_i)); (ws)(alpha_i) = w(alpha_i) - A[i][s] w(alpha_s)
        q = list(p); ps = p[s]
        for i in range(self.n):
            if self.A[i][s]:
                q[i] = p[i] - self.A[i][s]*ps
        return q
    def extend_is_fc(self, p, s):
        """w FC with height vector p, s not a right descent: is ws FC?"""
        q = list(p)
        J = self.J[s]
        while True:
            j = next((j for j in J if q[j] < 0), None)
            if j is None: break
            q = self.mul(q, j)
        for t in self.nb[s]:
            if q[t] < 0:
                qt = self.mul(q, t)
                if qt[s] < 0:
                    return False
        return True
    def word_is_reduced_fc(self, word):
        p = [1]*self.n
        for s in word:
            if p[s] <= 0: return (False, 'not reduced')
            if not self.extend_is_fc(p, s): return (False, 'not FC')
            p = self.mul(p, s)
        return (True, 'ok')
    def enumerate_fc(self):
        layer = {tuple([1]*self.n)}; total = 0; maxlen = 0
        while layer:
            total += len(layer); nxt = set()
            for p in layer:
                for s in range(self.n):
                    if p[s] > 0 and self.extend_is_fc(p, s):
                        nxt.add(tuple(self.mul(p, s)))
            if nxt: maxlen += 1
            layer = nxt
        return total, maxlen

def chain_cartan(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n): A[i][i] = 2
    for i in range(n-1): A[i][i+1] = A[i+1][i] = -1
    return A

def D_cartan(n):
    A = chain_cartan(n-1) ; A = [row+[0] for row in A] + [[0]*n]
    A[n-1][n-1] = 2; A[n-1][n-3] = A[n-3][n-1] = -1  # fork at node n-3
    return A

def family_word(n):
    p = n - 1
    w = list(range(n-2, 1, -1)) + [p]
    w += [1, 2, 3, 0, 1, 2, p]
    use0 = False
    for t in range(5, n-1, 2):
        for i in range(t-1, 1, -1):
            w += [i, i+1]
        last = (t == n-2)
        w += ([0, 1, 2] if (use0 and not last) else [1, 2]) + [p]
        use0 = not use0
    w += list(range(0, n-1))
    return w

if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if mode in ('validate', 'all'):
        print('validation of FC test (counts, max length):')
        for name, A, expect in [('A3', chain_cartan(3), 14), ('A4', chain_cartan(4), 42), ('A5', chain_cartan(5), 132),
                                ('D4', D_cartan(4), 48), ('D5', D_cartan(5), 167), ('E6', cartan(6), 662), ('E7', cartan(7), 2670)]:
            tot, ml = Group(A).enumerate_fc()
            print(f'  {name}: count={tot} expected={expect} ok={tot==expect} maxlen={ml}')
        sys.stdout.flush()
    if mode in ('lengths', 'all'):
        print('lengths of b_{r,k}:')
        for r in range(3, 8):
            for k in range(0, 3):
                if r >= 6 and k > 1: continue
                n = 4*r+1; A = cartan(n); beta = beta_rk(r, k)
                assert pair(A, beta, beta) == 2
                assert all(x > 0 for x in beta)
                L = length_of_matrix(A, refl_cols(A, beta))
                ht = sum(beta)
                print(f'  r={r} k={k} E{n}: l(b)={L}  8r^2+3+116k={8*r*r+3+116*k}  2ht-1={2*ht-1}  l(b)-1={L-1}  family=8r^2+7r-2={8*r*r+7*r-2}')
                sys.stdout.flush()
    if mode in ('family', 'all'):
        print('explicit family words:')
        for n in list(range(11, 30, 2)) + [33, 41, 73]:
            G = Group(cartan(n)); w = family_word(n)
            ok, why = G.word_is_reduced_fc(w)
            extra = ''
            if (n-1) % 4 == 0:
                r = (n-1)//4
                extra = f'  r={r}: 8r^2+7r-2={8*r*r+7*r-2}'
            print(f'  E{n}: len={len(w)} reduced&FC={ok} ({why}) support={sorted(set(w))==list(range(n))}{extra}')
            sys.stdout.flush()
