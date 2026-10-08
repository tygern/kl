from lib import data, pair, refl_matrix, matmul, length, apply
from collections import deque
from fractions import Fraction as Fr

# (1) coset-representative count #{alpha>0 in D_{4r} : (alpha,F omega)<0}, and check beta_r = s_0 F alpha_0 in the realization
for r in range(3,7):
    N=4*r
    s=[0]+[1 if i%2==1 else -1 for i in range(1,N+1)]  # sign pattern of 2*F*omega: +1 odd, -1 even
    cnt=0
    for i in range(1,N+1):
        for j in range(i+1,N+1):
            if s[j]-s[i]<0: cnt+=1    # e_j - e_i
            if s[i]+s[j]<0: cnt+=1    # e_i + e_j
    print(f"r={r}: #{{alpha>0: (alpha,F omega)<0}} = {cnt}, 4r^2={4*r*r}, r(2r+1)+r(2r-1)={r*(2*r+1)+r*(2*r-1)}")

# (2) lowering distance from beta_r to a simple root, BFS restricted to positive roots of height <= ht(beta_r)
for r in (3,4):
    n,A,beta,delta,gamma=data(r)
    ht=sum(beta)
    start=tuple(beta)
    dist={start:0}; dq=deque([start]); found=None
    while dq:
        v=dq.popleft()
        if sum(v)==1: found=dist[v]; break
        for j in range(n):
            p=sum(A[j][i]*v[i] for i in range(n))
            if p==0: continue
            w=list(v); w[j]-=p; w=tuple(w)
            if any(x<0 for x in w) or sum(w)>ht: continue
            if w not in dist:
                dist[w]=dist[v]+1; dq.append(w)
    l=8*r*r+3
    print(f"r={r}: lowering distance beta_r -> simple root (height-bounded BFS) = {found}; 4r^2+1={4*r*r+1}; "
          f"Brink-Howlett dp=min{{l(w):w beta<0}} <= {found+1}; l(b_r0)={l} = 2*{(l+1)//2}-1 so dp >= {(l+1)//2}")
