# E9 numbering: chain 0-1-...-7, node 8 attached to 2. Exact integer arithmetic.
n=9
edges={(i,i+1) for i in range(7)}|{(2,8)}
adj=lambda i,j: (i,j) in edges or (j,i) in edges
def pair(a,b):
    s=0
    for i in range(n):
        for j in range(n):
            if i==j: c=2
            elif adj(i,j): c=-1
            else: c=0
            s+=a[i]*c*b[j]
    return s
def simple(i): return tuple(1 if j==i else 0 for j in range(n))
def refl_s(i,v):  # s_i(v)
    c=pair(simple(i),v)
    return tuple(v[j]-(c if j==i else 0) for j in range(n))
def act_word(word,v):   # w = s_{word[0]} s_{word[1]} ... ; w(v) applies rightmost first
    for i in reversed(word): v=refl_s(i,v)
    return v
def refl_root(beta,v):
    c=pair(beta,v); assert pair(beta,beta)==2
    return tuple(v[j]-c*beta[j] for j in range(n))
beta=(1,2,3,3,2,2,1,1,2)
delta=(2,4,6,5,4,3,2,1,3)
assert pair(beta,beta)==2, "beta not a norm-2 real root"
print("beta norm:",pair(beta,beta), " delta null:", all(pair(delta,simple(i))==0 for i in range(n)))
word=[8,7,5,6,3,4,5,2,3,4,1,2,3,0,1,2,8,2,3,4,1,2,3,0,1,2,8,5,6,7,4,5,3]
M=[act_word(word,simple(j)) for j in range(n)]
Mr=[refl_root(beta,simple(j)) for j in range(n)]
print("word == r_beta:", M==Mr)
# descents: columns negative
def neg(v): return all(c<=0 for c in v) and any(c<0 for c in v)
def pos(v): return all(c>=0 for c in v) and any(c>0 for c in v)
R=[j for j in range(n) if neg(M[j])]
print("right descents:",R, " all columns signed:", all(neg(M[j]) or pos(M[j]) for j in range(n)))
# terminal on right: for t in R, s adjacent to t: w(alpha_s+alpha_t) positive
term=all(pos(tuple(M[s][i]+M[t][i] for i in range(n))) for t in R for s in range(n) if adj(s,t))
print("right-terminal:",term)
# involution -> left descents = right descents, left terminal also
I=[simple(j) for j in range(n)]
M2=[act_word(word+word,simple(j)) for j in range(n)]
print("involution:", M2==I)
# full support: row j != e_j for all j
print("full support (row test):", all(tuple(M[i][j] for i in range(n))!=simple(j) for j in range(n)))
# length via affine inversion count using finite E8 roots (omit node 7 = affine node)
import itertools
fin=[i for i in range(n) if i!=7]
# generate finite positive roots of E8 parabolic by closure
pos_roots=set(simple(i) for i in fin)
frontier=list(pos_roots)
while frontier:
    new=[]
    for r in frontier:
        for i in fin:
            c=pair(r,simple(i))
            if c<0:
                r2=tuple(r[j]+(1 if j==i else 0) for j in range(n))
                if r2 not in pos_roots: pos_roots.add(r2); new.append(r2)
    frontier=new
print("finite positive roots:",len(pos_roots))
allfin=list(pos_roots)+[tuple(-c for c in r) for r in pos_roots]
def decompose(v):
    # v = eta' + h*delta with eta' finite (coordinate 7 == 0)
    h=v[7]//delta[7]  # delta[7]=1
    etap=tuple(v[j]-h*delta[j] for j in range(n))
    assert etap[7]==0
    return etap,h
L=0
for eta in allfin:
    eps=0 if eta in pos_roots else 1
    img=act_word(word,eta)
    etap,h=decompose(img)
    epsp=0 if etap in pos_roots else 1
    L+=max(0,epsp-h-eps)
print("length by inversion count:",L, " word length:",len(word), " reduced:",L==len(word))
# Independence number of E9 diagram
best=0
for mask in range(1<<n):
    S=[i for i in range(n) if mask>>i&1]
    if all(not adj(a,b) for a in S for b in S if a<b): best=max(best,len(S))
print("independence number:",best," |R|=",len(R))
# bottoms: x=s3s5s7s8 etc; check Bruhat subword presence trivially and gaps
for bw in ([3,5,7,8],[0,3,5,7,8],[1,3,5,7,8],[0,1,3,5,7,8],[1,0,3,5,7,8]):
    print("bottom",bw,"length",len(bw),"gap",L-len(bw),"odd" if (L-len(bw))%2 else "even")
# extensions by s0 / s1
for g in (0,1):
    w2=word+[g]
    M3=[act_word(w2,simple(j)) for j in range(n)]
    R2=[j for j in range(n) if neg(M3[j])]
    L2=0
    for eta in allfin:
        eps=0 if eta in pos_roots else 1
        etap,h=decompose(act_word(w2,eta)); epsp=0 if etap in pos_roots else 1
        L2+=max(0,epsp-h-eps)
    t2=all(pos(tuple(M3[s][i]+M3[t][i] for i in range(n))) for t in R2 for s in range(n) if adj(s,t))
    # left descents via inverse word
    Minv=[act_word(list(reversed(w2)),simple(j)) for j in range(n)]
    Lset=[j for j in range(n) if neg(Minv[j])]
    tl=all(pos(tuple(Minv[s][i]+Minv[t][i] for i in range(n))) for t in Lset for s in range(n) if adj(s,t))
    print(f"extension s{g}: length {L2}, R={R2}, L={Lset}, right-terminal={t2}, left-terminal={tl}, gap to length-5 bottom = {L2-5}")
