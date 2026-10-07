"""Search affine-parabolic root translations preserving terminal-root cones."""
from targeted import en
from pathlib import Path
import json,argparse
p=argparse.ArgumentParser();p.add_argument('rank',type=int,default=13);p.add_argument('--max-entry',type=int,default=1);a=p.parse_args();n=a.rank;g=en(n)
delta=tuple([2,4,6,5,4,3,2,1]+[0]*(n-9)+[3])
def pair(x,y):return sum(x[i]*(2*y[i]-sum(y[j] for j in g.adj[i])) for i in range(n))
def reflroot(root,x):return tuple(v-pair(x,root)*r for v,r in zip(x,root))
def translate(gamma,x):
 gd=tuple(x+y for x,y in zip(gamma,delta))
 return reflroot(gamma,reflroot(gd,x))
def mvec(x):return tuple(pair(g.e[i],x) for i in range(n))
def finite_roots():
 roots=set(g.e[i] for i in list(range(7))+[n-1]);todo=list(roots)
 for root in todo:
  for i in list(range(7))+[n-1]:
   z=reflroot(g.e[i],root)
   if z not in roots:roots.add(z);todo.append(z)
 return sorted(roots)
roots=finite_roots();assert len(roots)==240
inputs=json.load(open(f'research/en_families/cartan_E{n}_m{a.max_entry}_max{int(a.max_entry>1)}.json'))['real_terminal_roots'];out=[]
for inp in inputs:
 beta=tuple(inp['beta']);m0=tuple(inp['pairings']);I=tuple(inp['I']);dirs=[]
 if not inp['maximum_independent']:continue
 for gamma in roots:
  b1=translate(gamma,beta);b2=translate(gamma,b1);b3=translate(gamma,b2)
  v=tuple(y-x for x,y in zip(beta,b1));u=tuple(z-2*y+x for x,y,z in zip(beta,b1,b2))
  if not any(v):continue
  assert tuple(z-2*y+x for x,y,z in zip(b1,b2,b3))==u
  # Integer k>=0 expansion beta+k*v+binom(k,2)*u; coefficientwise
  # nonnegative is sufficient for strictly positive roots for all k.
  if min(v)<0 or min(u)<0:continue
  mv,mu=mvec(v),mvec(u)
  if any(mv[i]<0 or mu[i]<0 for i in I):continue
  if any(mv[i]>0 or mu[i]>0 for i in range(n) if i not in I):continue
  if any(mv[s]+mv[t]>0 or mu[s]+mu[t]>0 for s in I for t in g.adj[s]):continue
  checks=[]
  for k in range(4):
   bk=tuple(x+k*y+(k*(k-1)//2)*z for x,y,z in zip(beta,v,u));mk=mvec(bk)
   w=tuple(tuple(int(i==j)-mk[j]*bk[i] for i in range(n)) for j in range(n))
   assert g.terminal(w) and g.desc(w)==I
   checks.append(dict(k=k,height=sum(bk),length=g.length(w)))
  dirs.append(dict(gamma=gamma,v=v,u=u,pairing_v=mv,pairing_u=mu,checks=checks))
 row=dict(rank=n,beta=beta,m=m0,I=I,delta=delta,directions=dirs);out.append(row);print(json.dumps(row),flush=True)
Path(__file__).with_name(f'indefinite_E{n}_translations.json').write_text(json.dumps(out,indent=2)+'\n')
