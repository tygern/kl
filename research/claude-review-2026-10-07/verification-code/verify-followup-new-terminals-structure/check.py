import sys
from collections import deque
# E_n manuscript labelling: chain 0-1-2-...-(n-2), node n-1 attached to node 2
def edges(n):
    E=set()
    for i in range(n-2): E.add((i,i+1)); E.add((i+1,i))
    E.add((2,n-1)); E.add((n-1,2))
    return E
def cartan(n):
    E=edges(n)
    return [[2 if i==j else (-1 if (i,j) in E else 0) for j in range(n)] for i in range(n)]
def act(A,s,v):  # s_s(v) = v - <v,alpha_s^vee> alpha_s ; coefficient = sum_j v_j A[s][j]
    c=sum(v[j]*A[s][j] for j in range(len(v)))
    w=list(v); w[s]-=c; return tuple(w)
def word_matrix(A,n,word):
    # columns: images of simple roots under w = s_{word[0]} s_{word[1]} ... ; apply rightmost first
    cols=[]
    for j in range(n):
        v=tuple(1 if k==j else 0 for k in range(n))
        for s in reversed(word): v=act(A,s,v)
        cols.append(v)
    return tuple(cols)
def length(M):  # number of positive roots sent negative; compute via root enumeration
    return None
def pos_roots(A,n):
    roots=set(); start=[tuple(1 if k==j else 0 for k in range(n)) for j in range(n)]
    dq=deque(start); roots.update(start)
    while dq:
        v=dq.popleft()
        for s in range(n):
            w=act(A,s,v)
            if all(x>=0 for x in w) and w not in roots: roots.add(w); dq.append(w)
    return roots
def apply(M,v):
    n=len(v); return tuple(sum(v[j]*M[j][k] for j in range(n)) for k in range(n))
def rdesc(M,A,n):
    return {s for s in range(n) if any(x<0 for x in apply(M,tuple(1 if k==s else 0 for k in range(n))))}
def is_right_terminal(M,A,n,E):
    for s in rdesc(M,A,n):
        for t in range(n):
            if (s,t) in E:
                v=tuple((1 if k in (s,t) else 0) for k in range(n))
                if any(x<0 for x in apply(M,v)): return False
    return True
def inverse(M,n):
    # M columns are images of simple roots; as matrix M[j][k]: image of alpha_j has coordinate k. inverse via solving: brute force use word? simpler: compute via sympy-free integer inversion using adjugate is annoying; use the fact that M is in a finite group: inverse = M^(order-1). Instead compute left descents via transpose trick is wrong for non-simply-laced; E is simply laced so left descents of w = right descents of w^{-1}; w^{-1} acts... use matrix inverse via fractions.
    from fractions import Fraction
    N=[[Fraction(M[j][k]) for k in range(n)]+[Fraction(1 if j==k else 0) for k in range(n)] for j in range(n)]
    for c in range(n):
        p=next(r for r in range(c,n) if N[r][c]!=0); N[c],N[p]=N[p],N[c]
        pv=N[c][c]; N[c]=[x/pv for x in N[c]]
        for r in range(n):
            if r!=c and N[r][c]!=0:
                f=N[r][c]; N[r]=[a-f*b for a,b in zip(N[r],N[c])]
    return tuple(tuple(int(N[j][n+k]) for k in range(n)) for j in range(n))
def is_terminal(M,A,n,E):
    return is_right_terminal(M,A,n,E) and is_right_terminal(inverse(M,n),A,n,E)
def is_commuting_product(M,A,n,E):
    D=rdesc(M,A,n)
    if any((s,t) in E for s in D for t in D): return False
    return M==word_matrix(A,n,sorted(D))
def support_of_word(word): return set(word)

rows={6:["1325213"],7:["1326213","13256213","132543621324356","1325436210321432543621324356"],
      8:["1327213","13257213","61327213","132543721324357","1325437210321432543721324357",
         "7534231270123456210321432"+"5437210321432543721324357"]}
# Gern words in Gern labels (1..n), D_n: 1,2 attached to 3, chain 3..n
gern={"w4":"2143214","w6":"214321654321546"}
def gern_to_E(word,n):  # Gern 1->1, 2->n-1, j->j-1
    return [1 if c=='1' else (n-1 if c=='2' else int(c)-1) for c in word]
for n in (6,7,8):
    A=cartan(n); E=edges(n); PR=pos_roots(A,n); assert len(PR)=={6:36,7:63,8:120}[n]
    def ell(M): return sum(1 for r in PR if any(x<0 for x in apply(M,r)))
    print("== E%d"%n)
    rowM={}
    for r in rows[n]:
        word=[int(c) for c in r]; M=word_matrix(A,n,word); rowM[r]=M
        print(" row",r,"len",ell(M),"reduced",ell(M)==len(word),"supp",sorted(set(word)),"terminal",is_terminal(M,A,n,E),"commprod",is_commuting_product(M,A,n,E))
    # Gern elements mapped into the D_{n-1} parabolic {1..n-2, n-1}
    cands={}
    for name,gw in gern.items():
        if name=="w6" and n==6: continue
        cands[name]=word_matrix(A,n,gern_to_E(gw,n))
    # w4 * u with u commuting generators outside supp(w4)={1,2,3,n-1}, in the D parabolic (nodes 4..n-2) and also not node 0 (outside D)
    for u in range(4,n-1):
        cands["w4*s%d"%u]=word_matrix(A,n,gern_to_E(gern["w4"],n)+[u])
    for name,M in cands.items():
        match=[r for r,RM in rowM.items() if RM==M]
        print(" Gern",name,"-> row",match, "len",ell(M),"terminal",is_terminal(M,A,n,E))
    # exhaustive: noncommuting terminals in the D_{n-1} parabolic (omit node 0)
    gens=[g for g in range(n) if g!=0]
    I=tuple(tuple(1 if j==k else 0 for k in range(n)) for j in range(n))
    seen={I}; dq=deque([I]); cnt=0
    while dq:
        M=dq.popleft()
        for s in gens:
            # right multiply by s: columns of Ms = images under M of s(alpha_j)
            Ms=tuple(apply(M,act(A,s,tuple(1 if k==j else 0 for k in range(n)))) for j in range(n))
            if Ms not in seen: seen.add(Ms); dq.append(Ms)
    print(" D%d parabolic size"%(n-1),len(seen))
    bad=[M for M in seen if is_terminal(M,A,n,E) and not is_commuting_product(M,A,n,E)]
    print(" noncommuting terminals in D parabolic:",len(bad))
    for M in bad:
        match=[r for r,RM in rowM.items() if RM==M]
        print("   len",ell(M),"row",match, "fixes node 0 & n-2?", )
    # also: which table rows lie in the D parabolic
    for r,RM in rowM.items():
        print(" row",r,"in D parabolic:",RM in seen)
