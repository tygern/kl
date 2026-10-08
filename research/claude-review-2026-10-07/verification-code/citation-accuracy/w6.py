# Signed permutations on {±1..±n}; w as tuple (w(1),...,w(n)); composition (uv)(i)=u(v(i))
n=6
def gen_gern(i):
    # Gern: s1 -> (1,-2)(-1,2); s_i -> (i-1,i)(-(i-1),-i)
    w=list(range(1,n+1))
    if i==1:
        w[0]=-2; w[1]=-1
    else:
        w[i-2],w[i-1]=w[i-1],w[i-2]
    return tuple(w)
def ap(w,i):
    return w[i-1] if i>0 else -w[-i-1]
def mul(u,v):
    return tuple(ap(u,v[i]) for i in range(n))
def word(ws):
    w=tuple(range(1,n+1))
    for s in ws: w=mul(w,gen_gern(s))
    return w
# Gern Lemma 2.3.4 for n=6, k=1: w6=[2,0][4,0][6,0][5,4][6,6]
# [j,0]=s_j...s_2 s_1 ; [5,4]=s5 s4 ; [6,6]=s6
w6 = word([2,1, 4,3,2,1, 6,5,4,3,2,1, 5,4, 6])
print("Gern w6 =",w6)
x6 = word([1,2,4,6]); print("x6 =",x6)
# manuscript E7 word 132543621324356 with E7->Gern map 1->1,6->2,2->3,3->4,4->5,5->6
m={1:1,6:2,2:3,3:4,4:5,5:6}
e7="132543621324356"
print("E7 word in Gern labels:", [m[int(c)] for c in e7])
print("manuscript w =", word([m[int(c)] for c in e7]))
# E8 word 132543721324357 with E8->Gern: 1->1,7->2,2->3,3->4,4->5,5->6
m8={1:1,7:2,2:3,3:4,4:5,5:6}
print("E8 word ->", word([m8[int(c)] for c in "132543721324357"]))
# length via inversions in D_n: #{i<j: w(i)>w(j)} + #{i<=j... } use Gern Prop 2.2.2 formula: sum_i |{j>i: a_i>a_j}| + |{j>i: -a_i>a_j}|
def length(w):
    L=0
    for i in range(n):
        for j in range(i+1,n):
            if w[i]>w[j]: L+=1
            if -w[i]>w[j]: L+=1
    return L
print("len w6 =",length(w6),"len x6 =",length(x6))
