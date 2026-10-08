# E8 parabolic on nodes J=(0..7, b) of E_{4r+1}, b=4r attached to node 2.
# Local labelling: indices 0..7 = chain nodes 0..7, index 8 = branch node b.
A=[[2 if i==j else 0 for j in range(9)] for i in range(9)]
def e(i,j): A[i][j]=-1; A[j][i]=-1
for i in range(7): e(i,i+1)
e(2,8)
def pr(u,v): return sum(u[i]*A[i][j]*v[j] for i in range(9) for j in range(9))
# generate E8 roots by closure under simple reflections from simple roots
FIN=[0,1,2,3,4,5,6,8]  # finite E8: omit affine node 7
simple=[tuple(1 if i==j else 0 for i in range(9)) for j in FIN]
roots=set(simple); frontier=list(simple)
while frontier:
    new=[]
    for v in frontier:
        for j in FIN:
            p=sum(A[j][i]*v[i] for i in range(9))
            w=tuple(v[i]-(p if i==j else 0) for i in range(9))
            if w not in roots: roots.add(w); new.append(w)
    frontier=new
assert len(roots)==240, len(roots)
gamma=(1,2,3,2,2,1,1,0,1)
assert gamma in roots
# pairing with beta_r: (eta,beta_r) = sum_j eta_j m_j, m=(r-2,2-r,-1,1,-1,1,-1,1, m_b=1)
# Also (delta,beta_r) = -a = -(2r-4).  Condition for alpha=eta+m*delta in N(T^k):
# (alpha,beta_r) = (eta,beta_r) - a*m <= -1 for all m>=eps(eta) (eps=0 if eta>0 else 1), a>=2.
# worst case m=eps:  eta>0: (eta,beta_r)<=-1 ; eta<0: (eta,beta_r) <= a-1 = 2r-5.
def pairing_coeffs(eta):
    # (eta,beta_r) = (r-2)*(eta0-eta1) + f(eta)
    c=eta[0]-eta[1]
    f=-eta[2]+eta[3]-eta[4]+eta[5]-eta[6]+eta[7]+eta[8]
    return c,f
pos=[v for v in roots if all(x>=0 for x in v)]
assert len(pos)==120
viol=[]; cnt_pos=0; cnt_neg=0; tot=0
for eta in roots:
    pg=pr(eta,gamma)
    if pg<1: continue
    tot+=pg
    c,f=pairing_coeffs(eta)
    ispos=all(x>=0 for x in eta)
    if ispos:
        cnt_pos+=1
        # need (r-2)c+f <= -1 for all r>=3: c<=0 and c+f<=-1 (at r=3)
        ok = (c<=0) and (c+f<=-1)
    else:
        cnt_neg+=1
        # need (r-2)c+f <= 2r-5 = 2(r-2)-1  <=> (r-2)(c-2) + f+1 <= 0 for all r>=3: c<=2 and (c-2)+f+1<=0
        ok = (c<=2) and ((c-2)+f+1<=0)
    if not ok: viol.append((eta,c,f,ispos))
print("roots with (eta,gamma)>=1:",cnt_pos,"positive,",cnt_neg,"negative; sum of pairings (=l(T)) =",tot)
print("violations:",viol)
# exact verification of the pairing decomposition against the full E_{4r+1} form for r=3..6
import sys
sys.path.insert(0,'.')
from lib import data, pair
for r in range(3,7):
    n,Afull,beta,delta,gam=data(r)
    J=list(range(8))+[4*r]
    bad=0
    for eta in roots:
        full=[0]*n
        for loc,g in enumerate(J): full[g]=eta[loc]
        c,f=pairing_coeffs(eta)
        if pair(Afull,full,beta)!=(r-2)*c+f: bad+=1
    print(f"r={r}: pairing decomposition mismatches: {bad}; (delta,beta_r)={pair(Afull,delta,beta)} = -(2r-4)={-(2*r-4)}")
