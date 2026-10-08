"""Part (a): are the non-new rows of Table 1 exactly Gern's D-type bad elements?

1. Sanity-check my reading of Gern: bracket form (Lemma 2.3.4) as signed
   permutation equals Corollary 2.2.19 and has length 3n^2/8+n/4.
2. Map Gern's w4, w4*u, w6 into E_n via 1->1, 2->n-1, j->j-1 and compare
   matrices with the table words.
3. Exhaustively enumerate all bad (noncommuting terminal) elements of D7 via
   signed permutations; restrict to D5, D6 parabolics; compare with Gern's
   Theorem 2.3.6 and with the table.
"""
from itertools import permutations, product
import sys
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-new-terminals-structure')
from en_tools import *

TABLE = [
    (6, '1325213', '135'),
    (7, '1326213', '136'),
    (7, '13256213', '1356'),
    (7, '132543621324356', '1356'),
    (7, '1325436210321432543621324356', '1356'),
    (8, '1327213', '137'),
    (8, '13257213', '1357'),
    (8, '61327213', '1367'),
    (8, '132543721324357', '1357'),
    (8, '1325437210321432543721324357', '1357'),
    (8, '7534231270123456210321432' + '5437210321432543721324357', '1357'),
]

# ---------- 1. Gern sanity ----------
def gern_cor_2_2_19(n, m):
    """Signed permutation of w_n inside D_m (m >= n), Corollary 2.2.19."""
    w = []
    for i in range(1, m + 1):
        if i > n:
            w.append(i)
        elif n % 2 == 0:
            if i == 1:
                w.append((-1) ** (n // 2))
            elif i % 2 == 1:
                w.append(i)
            else:
                w.append(-(n + 2 - i))
        else:
            if i == 1:
                w.append((-1) ** ((n - 1) // 2))
            elif i % 2 == 1:
                w.append(i)
            else:
                w.append(-(n + 1 - i))
    return tuple(w)

print("== Gern bracket form vs Corollary 2.2.19 ==")
for n in range(4, 11):
    word = gern_wn_word(n)
    sp = signed_perm_of_gern_word(word, n)
    cor = gern_cor_2_2_19(n, n)
    ne = n if n % 2 == 0 else n - 1
    expected_len = (3 * ne * ne + 2 * ne) // 8
    print(f"w_{n}: word={word_str(word)} len(word)={len(word)} expected ell={expected_len} "
          f"signed perm={sp} matches Cor 2.2.19: {sp == cor}")
    assert sp == cor and len(word) == expected_len

# ---------- 2. Map into E_n and compare with table ----------
print("\n== Gern elements mapped into E_n (1->1, 2->n-1, j->j-1) ==")
def mapped_matrix(rs, n_en, gern_word):
    return rs.word_matrix([gern_to_en(s, n_en) for s in gern_word])

results = {}
for n_en in (6, 7, 8):
    rs = RootSystem(n_en, en_edges(n_en))
    m = n_en - 1  # D_m parabolic on nodes 1..n_en-1
    gern_list = []  # (description, gern word)
    # w4 * u, u commuting product inside {s6,...,s_m}
    rest = list(range(6, m + 1))
    for r in range(len(rest) + 1):
        for U in combinations(rest, r):
            if any(abs(a - b) == 1 for a in U for b in U):
                continue
            gern_list.append((f"w4*{''.join(map(str,U)) or 'e'}", gern_wn_word(4) + list(U)))
    if m >= 6:
        rest = list(range(8, m + 1))
        for r in range(len(rest) + 1):
            for U in combinations(rest, r):
                if any(abs(a - b) == 1 for a in U for b in U):
                    continue
                gern_list.append((f"w6*{''.join(map(str,U)) or 'e'}", gern_wn_word(6) + list(U)))
    if m >= 8:
        gern_list.append(("w8", gern_wn_word(8)))
    table_rows = [(w, x) for (nn, w, x) in TABLE if nn == n_en]
    table_mats = {w: rs.word_matrix(parse(w)) for w, x in table_rows}
    print(f"\nE_{n_en} (D_{m} parabolic): Gern predicts {len(gern_list)} bad elements")
    matched = set()
    for desc, gw in gern_list:
        M = mapped_matrix(rs, n_en, gw)
        mapped_word = word_str([gern_to_en(s, n_en) for s in gw])
        hits = [w for w in table_mats if table_mats[w] == M]
        ell = rs.length(M)
        print(f"  {desc:10s} Gern word {word_str(gw):20s} -> E_n word {mapped_word:30s} "
              f"ell={ell} terminal={rs.terminal(M)} L={sorted(rs.left_descents(M))} "
              f"R={sorted(rs.right_descents(M))} matches table row: {hits}")
        assert len(hits) == 1
        matched.update(hits)
    unmatched = [w for w in table_mats if w not in matched]
    print(f"  Table rows NOT accounted for by Gern: {unmatched}")
    results[n_en] = unmatched

# ---------- 3. Exhaustive D7 enumeration via signed permutations ----------
print("\n== Exhaustive enumeration of bad elements in D7 (signed permutations, Gern labels) ==")
m = 7
# Realize D7 roots: alpha_1 = e1+e2, alpha_i = e_i - e_{i-1} (i>=2), in Gern's convention
# s1 = (1,-2)(-1,2): reflection in e1+e2.  s_i = (i-1,i): reflection in e_i - e_{i-1}.
# Positive roots: e_j - e_i (i<j), e_i + e_j (i<j).  w(e_i) = sign * e_{|w(i)|}.
def image_root(w, root):
    # root as dict {index: coeff}
    out = {}
    for i, c in root.items():
        j = w[i - 1]
        out[abs(j)] = out.get(abs(j), 0) + (c if j > 0 else -c)
    return out

pos_roots = []
for i in range(1, m + 1):
    for j in range(i + 1, m + 1):
        pos_roots.append({j: 1, i: -1})
        pos_roots.append({i: 1, j: 1})

def root_is_positive(r):
    # positive iff in the span with nonneg simple-root coefficients:
    # e_j - e_i (i<j) positive; e_i + e_j positive; negatives otherwise.
    items = [(k, v) for k, v in r.items() if v != 0]
    assert len(items) == 2
    (a, ca), (b, cb) = items
    if ca > 0 and cb > 0:
        return True
    if ca < 0 and cb < 0:
        return False
    # one +1, one -1: e_hi - e_lo positive
    hi = a if ca > 0 else b
    lo = b if ca > 0 else a
    return hi > lo

simple_roots = {1: {1: 1, 2: 1}}
for i in range(2, m + 1):
    simple_roots[i] = {i: 1, i - 1: -1}

def descents_right(w):
    return {s for s in range(1, m + 1) if not root_is_positive(image_root(w, simple_roots[s]))}

def inv(w):
    out = [0] * m
    for i, j in enumerate(w, 1):
        out[abs(j) - 1] = i if j > 0 else -i
    return tuple(out)

adj = {1: {3}, 2: {3}}
for i in range(3, m + 1):
    adj[i] = {i - 1, i + 1} if i < m else {i - 1}
adj[3] = {1, 2, 4}

def right_terminal(w):
    R = descents_right(w)
    for s in R:
        for t in adj[s]:
            r = {k: simple_roots[s].get(k, 0) + simple_roots[t].get(k, 0) for k in set(simple_roots[s]) | set(simple_roots[t])}
            r = {k: v for k, v in r.items() if v}
            if not root_is_positive(image_root(w, r)):
                return False
    return True

def length(w):
    return sum(1 for r in pos_roots if not root_is_positive(image_root(w, r)))

def is_commuting_product(w):
    R = descents_right(w)
    if any(t in R for s in R for t in adj[s]):
        return False
    return length(w) == len(R) and signed_perm_of_gern_word(sorted(R), m) == w

count = 0
bad = []
for perm in permutations(range(1, m + 1)):
    for signs in product((1, -1), repeat=m):
        if sum(1 for s in signs if s < 0) % 2:
            continue
        w = tuple(s * p for s, p in zip(signs, perm))
        count += 1
        if right_terminal(w) and right_terminal(inv(w)) and not is_commuting_product(w):
            bad.append(w)
print(f"elements enumerated: {count} (expected {2**(m-1)*720*7})")
print(f"bad elements in D7: {len(bad)}")
gern_pred = {
    'w4': gern_cor_2_2_19(4, m),
    'w4*s6': signed_perm_of_gern_word(gern_wn_word(4) + [6], m),
    'w4*s7': signed_perm_of_gern_word(gern_wn_word(4) + [7], m),
    'w6=w7': gern_cor_2_2_19(6, m),
}
for w in sorted(bad, key=length):
    names = [k for k, v in gern_pred.items() if v == w]
    # support: generators appearing in a reduced word
    supp = set()
    ww = w
    while descents_right(ww):
        s = min(descents_right(ww))
        supp.add(s)
        g = signed_perm_of_gern_word([s], m)
        ww = tuple((1 if x > 0 else -1) * ww[abs(x) - 1] for x in g)  # ww * s
    print(f"  {w} ell={length(w)} supp={sorted(supp)} Gern name: {names}")
assert set(bad) == set(gern_pred.values())
print("D7 bad elements == {w4, w4 s6, w4 s7, w6}: True")
print("restricted to D6 (fixes 7):", sum(1 for w in bad if w[6] == 7))
print("restricted to D5 (fixes 6,7):", sum(1 for w in bad if w[5] == 6 and w[6] == 7))
