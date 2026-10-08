import sys
from fc_check import adjacency_E, is_fc_reduced, length_by_matrix

def Z2(T, p):   # width-2 zigzag 2 1 3 2 ... T T-1, closed by p
    w = []
    for j in range(2, T+1): w += [j, j-1]
    return w + [p]
def Z2z(T, p):  # width-2 zigzag with an extra 0 at the bottom: 2 1 0 3 2 4 3 ... T T-1 p
    w = [2, 1, 0]
    for j in range(3, T+1): w += [j, j-1]
    return w + [p]
def Z3(T, p):   # width-3 zigzag 2 1 0 3 2 1 ... T T-1 T-2, closed by p
    w = []
    for j in range(2, T+1): w += [j, j-1, j-2]
    return w + [p]

def best_block_word(n):
    """DFS over block sequences A | Z2(n-2) | blocks with decreasing tops | final run 2..n-2,
    keeping only FC words; returns the longest found."""
    p = n-1; adj = adjacency_E(n)
    A = list(range(n-2, -1, -1)) + [p]
    final = list(range(2, n-1))
    best = [None]
    def dfs(word, maxtop):
        cand = word + final
        if is_fc_reduced(cand, adj) and (best[0] is None or len(cand) > len(best[0])):
            best[0] = cand
        for T in range(maxtop, 2, -1):
            for blk in (Z3(T,p), Z2(T,p), Z2z(T,p)):
                w2 = word + blk
                if is_fc_reduced(w2, adj):
                    dfs(w2, T-1)
    dfs(A + Z2(n-2, p), n-3)
    return best[0]

for n in [int(x) for x in sys.argv[1:]]:
    w = best_block_word(n)
    r = (n-1)/4
    print(f"n={n}: block-structured FC length {len(w)}   (8r^2+2={8*r*r+2:.0f}, 8r^2+7r-2={8*r*r+7*r-2:.0f})")
    print("   word:", ' '.join(map(str, w)))
