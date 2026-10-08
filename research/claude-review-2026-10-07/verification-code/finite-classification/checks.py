import itertools
def graph(n):
    E=set()
    for i in range(n-2): E.add((i,i+1))
    E.add((2,n-1))
    return E
def adj(E,i,j): return (i,j) in E or (j,i) in E
def indep_sets(nodes,E):
    cnt=0; mx=0
    for r in range(len(nodes)+1):
        for S in itertools.combinations(nodes,r):
            if all(not adj(E,a,b) for a,b in itertools.combinations(S,2)):
                cnt+=1; mx=max(mx,r)
    return cnt,mx
for n in (6,7,8):
    E=graph(n); print("E%d: #independent sets (incl empty) = %d, independence number = %d (ceil(n/2)=%d)"%((n,)+indep_sets(list(range(n)),E)+(-(-n//2),)))
for n in range(6,13):
    E=graph(n); c,m=indep_sets(list(range(n)),E); print("  E%d independence number %d vs ceil(n/2)=%d"%(n,m,-(-n//2)))
# induced support graphs in the table
rows={ 'E6 len7':(6,'1325213'), 'E7 len7':(7,'1326213'),'E7 len8':(7,'13256213'),'E7 len15':(7,'132543621324356'),'E7 len28':(7,'1325436210321432543621324356'),
 'E8 len7':(8,'1327213'),'E8 len8a':(8,'13257213'),'E8 len8b':(8,'61327213'),'E8 len15':(8,'132543721324357'),'E8 len28':(8,'1325437210321432543721324357'),'E8 len50':(8,'75342312701234562103214325437210321432543721324357')}
for k,(n,w) in rows.items():
    E=graph(n); supp=sorted(set(int(c) for c in w)); c,m=indep_sets(supp,E)
    print(k, "supp",supp,"independence number",m, "#edges in induced graph", sum(1 for a,b in itertools.combinations(supp,2) if adj(E,a,b)))
# Cartan matrix E8 and the weight (4,7,10,8,6,4,2,5)
n=8;E=graph(n)
C=[[2 if i==j else (-1 if adj(E,i,j) else 0) for j in range(n)] for i in range(n)]
v=(4,7,10,8,6,4,2,5)
print("C*(4,7,10,8,6,4,2,5) =", [sum(C[i][j]*v[j] for j in range(n)) for i in range(n)])
import fractions
# number of cosets E8/E7, E8/D7, E8/A7
print("E8/E7 index", 696729600//2903040, " E8/D7 index", 696729600//(2**6*5040), " E8/A7 index", 696729600//40320)
# Gern's D6 signed permutation check. Gern labels: s1=(1,-2)(-1,2), s_i=(i-1,i) for i>=2; right mult by s_i acts on positions.
def gern_perm(word,n=6):
    w=list(range(1,n+1))
    for i in word:
        if i==1:
            w[0],w[1]=-w[1],-w[0]
        else:
            w[i-2],w[i-1]=w[i-1],w[i-2]
    return tuple(w)
# paper mapping: Gern labels 1,2,3,4,5,6 -> paper 1,n-1,2,3,4,5 ; inverse: paper 1->1, n-1->2, 2->3, 3->4, 4->5, 5->6
def to_gern(word,n):
    inv={1:1,n-1:2,2:3,3:4,4:5,5:6}
    return [inv[int(c)] for c in word]
for n,w,x in ((7,'132543621324356','1356'),(8,'132543721324357','1357')):
    gw=to_gern(w,n); gx=to_gern(x,n)
    print("E%d: Gern word for b ="%n, gw, "-> signed perm", gern_perm(gw), " x ->", gx, gern_perm(gx))
print("Gern Cor 2.2.19 w_6 =", tuple((-1) if i==1 else (i if i%2==1 else -(6+2-i)) for i in range(1,7)))
# length of w_6 via Gern Prop 2.2.2
def dlen(w):
    n=len(w); return sum(1 for i in range(n) for j in range(i+1,n) if w[i]>w[j]) + sum(1 for i in range(n) for j in range(i+1,n) if -w[i]>w[j])
print("length of (-1,-6,3,-4,5,-2) by inversion count:", dlen((-1,-6,3,-4,5,-2)), " 3n^2/8+n/4 =", 3*36/8+6/4)
# parity of table gaps
for g in (4,4,4,11,24,4,4,4,11,24,46): pass
print("odd gaps among table rows:", [g for g in (4,4,4,11,24,4,4,4,11,24,46) if g%2])
