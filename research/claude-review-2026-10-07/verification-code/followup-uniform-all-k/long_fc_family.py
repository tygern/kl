"""Explicit long FC elements of E_n (n odd), generalizing the computed maximum-length
elements of E_11, E_13:  A = (n-2)(n-3)...(3)(2) p ;  B_1 = 1 2 3 0 1 2 p ;
B(t) for t = 5,7,...,n-2: runs (t-1 t)(t-2 t-1)...(2 3) then (1 2) or (0 1 2), then p,
with the bottom alternating 12 / 012 / 12 / ...;  final run 0 1 2 ... (n-2).
Each word is checked to be reduced (length via root matrices) and FC (heap criterion),
and compared with l(b_{r,0}) - 1 = 8r^2 + 2 and l(b_{r,k}) - 1 at n = 4r + 1."""
import sys
sys.path.insert(0, '/private/tmp/claude-501/-Users-ic-workspace-kl/45cb4eb2-517e-46a6-a21d-77de64370e52/scratchpad/followup-uniform-all-k')
from uni import En


def family_word(n):
    p = n - 1
    w = list(range(n - 2, 1, -1)) + [p]          # A
    w += [1, 2, 3, 0, 1, 2, p]                     # B_1 (top 3)
    use0 = False
    for t in range(5, n - 1, 2):                   # tops 5,7,...,n-2
        for i in range(t - 1, 1, -1):              # runs (i i+1) for i = t-1 .. 2
            w += [i, i + 1]
        last = (t == n - 2)
        w += ([0, 1, 2] if (use0 and not last) else [1, 2]) + [p]
        use0 = not use0
    w += list(range(0, n - 1))                     # final run
    return w


if __name__ == '__main__':
    ns = [int(a) for a in sys.argv[1:]] or list(range(11, 62, 2))
    for n in ns:
        G = En(n)
        w = family_word(n)
        x = G.from_word(w)
        red = (G.length(x) == len(w))
        fc = G.is_fc_word(w)
        line = f"E{n}: family length {len(w)} reduced={red} FC={fc}"
        if (n - 1) % 4 == 0:
            r = (n - 1) // 4
            k0 = 0
            while 8 * r * r + 2 + 116 * k0 < len(w):
                k0 += 1
            line += (f"  [r={r}: l(b_r0)-1={8*r*r+2}, excess={len(w)-(8*r*r+2)}; "
                     f"smallest k with l(b_rk)-1 >= family length: {k0}]")
        print(line)
        sys.stdout.flush()
