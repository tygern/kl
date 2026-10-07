"""Exact translation directions preserving specified affine terminal conditions."""
from targeted import en
from pathlib import Path
import json
G=en(9);delta=(2,4,6,5,4,3,2,1,3)
def pair(a,b):return sum(a[i]*(2*b[i]-sum(b[j] for j in G.adj[i])) for i in range(9))
def apply(w,a):return tuple(sum(w[j][i]*a[j] for j in range(9)) for i in range(9))
def reflect(a,s):
 out=list(a);out[s]-=pair(G.e[s],a);return tuple(out)
roots=set(G.e[i] for i in range(9) if i!=7);todo=list(roots)
for a in todo:
 for s in range(9):
  if s==7:continue
  b=reflect(a,s)
  if b not in roots:roots.add(b);todo.append(b)
rows=json.load(open('research/en_families/e9_full_support_eligible.json'));out=[]
for row in rows:
 w=G.elt(row['word']);I=G.desc(w);dirs={}
 for gamma in roots:
  wg=apply(w,gamma);v=tuple(x-y for x,y in zip(gamma,wg));d=tuple(pair(G.e[j],v) for j in range(9))
  if not any(d):continue
  if any(d[s]>0 for s in I) or any(d[s]<0 for s in range(9) if s not in I):continue
  if any(d[s]+d[t]<0 for s in I for t in G.adj[s]):continue
  if d in dirs:continue
  lengths=[];checks=[]
  U=tuple(tuple(int(i==j)+d[j]*delta[i] for i in range(9)) for j in range(9))
  lu=G.length(U)
  for k in range(4):
   wk=tuple(tuple(w[j][i]+k*d[j]*delta[i] for i in range(9)) for j in range(9))
   assert G.terminal(wk) and G.desc(wk)==I
   lengths.append(G.length(wk));checks.append(G.length(wk)==G.length(w)+k*lu)
  dirs[d]=dict(gamma=gamma,shift_vector=v,column_slopes=d,translation_length=lu,lengths=lengths,length_additive=all(checks))
 result=dict(base_length=row['length'],word=row['word'],directions=list(dirs.values()))
 out.append(result);print(json.dumps(result),flush=True)
Path(__file__).with_name('translation_cones.json').write_text(json.dumps(out,indent=2)+'\n')
