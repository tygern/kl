# Gern Example 1.1.6 / Prop 2.2.3 (right action): s1: (w1,w2)->(-w2,-w1); s_j (j>=2): swap positions j-1,j
def rmul(w, j):
    w = list(w)
    if j == 1:
        w[0], w[1] = -w[1], -w[0]
    else:
        w[j-2], w[j-1] = w[j-1], w[j-2]
    return tuple(w)
def length(w):
    n = len(w)
    return sum(1 for i in range(n) for j in range(i+1, n) if w[i] > w[j]) + \
           sum(1 for i in range(n) for j in range(i+1, n) if w[i] + w[j] < 0)
def word_to_perm(word, n):
    w = tuple(range(1, n+1))
    for j in word:
        w = rmul(w, j)
    return w
# Example 2.2.1: s2 s3 s1 s2 s4 in D4 should be (-2,-3,4,1)
print("Ex 2.2.1:", word_to_perm([2,3,1,2,4], 4))
# Lemma 2.3.4 n=8: [2,0][4,0][6,0][8,0][6,4][7,6][8,8]
def br(j, i):
    # [j,i] = [i,j]^{-1} with [0,j]=s1 s2 ... sj, [i,j]=s_i...s_j (i>=2)
    if i == 0:
        return list(range(j, 0, -1))
    return list(range(j, i-1, -1))
word8 = br(2,0)+br(4,0)+br(6,0)+br(8,0)+br(6,4)+br(7,6)+br(8,8)
print("word8 =", word8, "len", len(word8))
w8 = word_to_perm(word8, 8)
print("w8 =", w8, "length", length(w8))
print("target (1,-8,3,-6,5,-4,7,-2) length", length((1,-8,3,-6,5,-4,7,-2)))
x8 = word_to_perm([1,2,4,6,8], 8)
print("x8 =", x8, "length", length(x8))
word6 = br(2,0)+br(4,0)+br(6,0)+br(5,4)+br(6,6)
w6 = word_to_perm(word6, 6); print("w6 =", w6, length(w6), "x6 =", word_to_perm([1,2,4,6],6))
# right descents by Prop 2.2.4
def rdes(w):
    return [j for j in range(1, len(w)+1) if (j==1 and w[0]+w[1]<0) or (j>=2 and w[j-2]>w[j-1])]
print("R(w8) =", rdes(w8), "R(x8) =", rdes(x8))
inv = [0]*8
for i,v in enumerate(w8): inv[abs(v)-1] = (i+1)*(1 if v>0 else -1)
print("w8^{-1} =", tuple(inv))
