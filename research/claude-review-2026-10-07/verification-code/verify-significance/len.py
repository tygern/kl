import sys
from fractions import Fraction
def En(n):
    # nodes 0..n-1; chain 0-1-...-(n-2), node n-1 attached to 2
    A=[[0]*n for _ in range(n)]
    for i in range(n): A[i][i]=2
    for i in range(n-2):
        A[i][i+1]=A[i+1][i]=-1
    A[2][n-1]=A[n-1][2]=-1
    return A
def pair(A,u,v):
    n=len(A); return sum(u[i]*A[i][j]*v[j] for i in range(n) for j in range(n))
def refl_matrix(A,beta):
    n=len(A)
    # columns: image of alpha_j = alpha_j - (alpha_j,beta) beta
    cols=[]
    for j in range(n):
        m=sum(A[j][i]*beta[i] for i in range(n))
        col=[(1 if i==j else 0)-m*beta[i] for i in range(n)]
        cols.append(col)
    return cols
def length(A,cols):
    n=len(A); cols=[c[:] for c in cols]; L=0; word=[]
    while True:
        ds=[j for j in range(n) if all(x<=0 for x in cols[j]) and any(x<0 for x in cols[j])]
        if not ds:
            assert all(all(x>=0 for x in c) for c in cols)
            # must be identity
            assert all(cols[j][i]==(1 if i==j else 0) for j in range(n) for i in range(n))
            return L,word
        s=ds[0]
        # w <- w s : new columns: (ws)(alpha_j) = w(alpha_j - A[s][j] alpha_s) = w alpha_j - A[s][j] w alpha_s
        ws=cols[s]
        new=[]
        for j in range(n):
            new.append([cols[j][i]-A[s][j]*ws[i] for i in range(n)])
        cols=new; L+=1; word.append(s)
def beta_rk(r,k):
    n=4*r+1; b=[0]*n
    b[0]=r-1;b[1]=r;b[2]=2*r-1;b[n-1]=r
    for j in range(3,4*r): b[j]=2*r-j//2
    a=2*r-4
    delta=[0]*n; gamma=[0]*n
    for j,v in enumerate((2,4,6,5,4,3,2,1)): delta[j]=v
    delta[n-1]=3
    for j,v in enumerate((1,2,3,2,2,1,1,0)): gamma[j]=v
    gamma[n-1]=1
    return [b[i]-a*k*gamma[i]+(a*k*k+r*k)*delta[i] for i in range(n)]
for r in (3,4):
    A=En(4*r+1)
    for k in range(0,4):
        beta=beta_rk(r,k)
        assert pair(A,beta,beta)==2
        L,word=length(A,refl_matrix(A,beta))
        n=4*r+1
        Rdes=sorted(set([word[-1]]))
        print(f"r={r} k={k} ht={sum(beta)} len={L} odd={L%2==1} supp={sorted(set(word))==list(range(n))}")
