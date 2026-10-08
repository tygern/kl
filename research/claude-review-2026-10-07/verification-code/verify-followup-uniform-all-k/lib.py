def cartan(r):
    n=4*r+1
    A=[[0]*n for _ in range(n)]
    for i in range(n): A[i][i]=2
    def edge(i,j): A[i][j]=-1; A[j][i]=-1
    for i in range(4*r-1): edge(i,i+1)
    edge(2,4*r)
    return A

def pair(A,u,v):
    n=len(A); return sum(u[i]*A[i][j]*v[j] for i in range(n) for j in range(n))

def refl_matrix(A,beta):
    n=len(A)
    # columns r_beta(alpha_j) = alpha_j - (alpha_j,beta) beta
    cols=[]
    for j in range(n):
        p=sum(A[j][i]*beta[i] for i in range(n))
        col=[(1 if i==j else 0)-p*beta[i] for i in range(n)]
        cols.append(col)
    return cols  # list of columns

def matmul(A_,X,Y):  # X,Y list of columns; (XY)(alpha_j) = X(Y alpha_j)
    n=len(X)
    out=[]
    for j in range(n):
        v=Y[j]
        col=[sum(X[i][a]*v[i] for i in range(n)) for a in range(n)]
        out.append(col)
    return out

def length(A,W):
    n=len(A); W=[c[:] for c in W]; L=0
    while True:
        j=next((j for j in range(n) if all(x<=0 for x in W[j]) and any(x<0 for x in W[j])),None)
        if j is None:
            assert all(W[i][a]==(1 if i==a else 0) for i in range(n) for a in range(n)), "not identity"
            return L
        # W <- W s_j : col_i <- col_i - A[j][i] col_j
        cj=W[j]
        W=[[W[i][a]-A[j][i]*cj[a] for a in range(n)] for i in range(n)]
        L+=1
        if L>100000: raise Exception("too long")

def data(r):
    n=4*r+1; A=cartan(r)
    beta=[0]*n
    beta[0]=r-1; beta[1]=r; beta[2]=2*r-1; beta[4*r]=r
    for j in range(3,4*r): beta[j]=2*r-j//2
    delta=[0]*n; gamma=[0]*n
    for j,v in enumerate([2,4,6,5,4,3,2,1]): delta[j]=v
    delta[4*r]=3
    for j,v in enumerate([1,2,3,2,2,1,1,0]): gamma[j]=v
    gamma[4*r]=1
    return n,A,beta,delta,gamma

def apply(W,v):
    n=len(W); return [sum(W[i][a]*v[i] for i in range(n)) for a in range(n)]

