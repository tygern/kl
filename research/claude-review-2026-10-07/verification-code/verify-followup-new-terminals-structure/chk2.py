from fractions import Fraction
n=8
edges={(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)}
def adj(i,j): return (i,j) in edges or (j,i) in edges
A=[[2 if i==j else (-1 if adj(i,j) else 0) for j in range(n)] for i in range(n)]
def I(): return [[1 if i==j else 0 for j in range(n)] for i in range(n)]
def mul(X,Y): return [[sum(X[i][k]*Y[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def mv(X,v): return [sum(X[i][k]*v[k] for k in range(n)) for i in range(n)]
def simple(i):
    M=I()
    for j in range(n): M[i][j]-=A[i][j]
    return M
S=[simple(i) for i in range(n)]
word="7534231270123456210321432"+"5437210321432543721324357"
print("len word",len(word))
W=I()
for c in word: W=mul(W,S[int(c)])
print("w^2==I:", mul(W,W)==I())
def rank(M):
    M=[[Fraction(x) for x in row] for row in M]; r=0; rows=len(M); cols=len(M[0])
    for c in range(cols):
        p=next((i for i in range(r,rows) if M[i][c]!=0),None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]
        for i in range(rows):
            if i!=r and M[i][c]!=0:
                f=M[i][c]/M[r][c]; M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        r+=1
    return r
Wp=[[W[i][j]+(1 if i==j else 0) for j in range(n)] for i in range(n)]
Wm=[[W[i][j]-(1 if i==j else 0) for j in range(n)] for i in range(n)]
print("dim(-1)-eig:", n-rank(Wp), " dim(+1)-eig:", n-rank(Wm))
# positive roots
fr=[tuple(1 if j==i else 0 for j in range(n)) for i in range(n)]
roots=set(fr); frontier=list(fr)
while frontier:
    new=[]
    for r in frontier:
        for i in range(n):
            t=tuple(mv(S[i],list(r)))
            if all(c>=0 for c in t) and t not in roots: roots.add(t); new.append(t)
    frontier=new
print("num pos roots",len(roots))
print("length =", sum(1 for r in roots if all(c<=0 for c in mv(W,list(r)))))
cand=[r for r in sorted(roots,key=sum) if mv(W,list(r))==[-c for c in r]]
print("roots negated by w:", [(r,sum(r)) for r in cand])
def refl(r):
    rA=[sum(r[k]*A[k][j] for k in range(n)) for j in range(n)]
    return [[(1 if i==j else 0)-r[i]*rA[j] for j in range(n)] for i in range(n)]
for a in cand:
    for b in cand:
        if a<b and mul(refl(a),refl(b))==W:
            ip=sum(a[i]*A[i][j]*b[j] for i in range(n) for j in range(n))
            print("W = r_a r_b with a=",''.join(map(str,a))," b=",''.join(map(str,b)),"(a,b)=",ip)
# descents
print("right descents:", [i for i in range(n) if any(c<0 for c in mv(W,[1 if j==i else 0 for j in range(n)]))])
# char poly via Faddeev-LeVerrier with Fractions
def charpoly(M):
    Mf=[[Fraction(x) for x in row] for row in M]; N=I(); N=[[Fraction(x) for x in row] for row in N]
    coeffs=[Fraction(1)]; Mk=[[Fraction(0)]*n for _ in range(n)]
    for k in range(1,n+1):
        Mk=mul(Mf,[[Mk[i][j]+coeffs[-1]*(1 if i==j else 0) for j in range(n)] for i in range(n)]) if k>1 else [row[:] for row in Mf]
        c=-sum(Mk[i][i] for i in range(n))/k
        coeffs.append(c)
    return coeffs
print("charpoly coeffs (x^8 ... x^0):", charpoly(W))
