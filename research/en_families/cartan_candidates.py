"""Exhaust specified small root pairings for terminal reflections in E_n.
All arithmetic exact; candidate root reality is certified by simple reflections.
"""
from targeted import en
from fractions import Fraction
from itertools import product
from pathlib import Path
import argparse,json,time,sys
sys.setrecursionlimit(20000)

def inverse_cartan(g):
 n=g.n;A=[[Fraction(2 if i==j else -1 if j in g.adj[i] else 0) for j in range(n)]+[Fraction(int(i==j)) for j in range(n)] for i in range(n)]
 for j in range(n):
  pivot=next(i for i in range(j,n) if A[i][j]);A[j],A[pivot]=A[pivot],A[j];v=A[j][j];A[j]=[x/v for x in A[j]]
  for i in range(n):
   if i==j:continue
   v=A[i][j]
   if v:A[i]=[x-v*y for x,y in zip(A[i],A[j])]
 return [row[n:] for row in A]

def independent_sets(g):
 out=[]
 def visit(i,chosen):
  if i==g.n:out.append(tuple(chosen));return
  visit(i+1,chosen)
  if all(j not in chosen for j in g.adj[i]):visit(i+1,chosen+[i])
 visit(0,[]);return out

def root_reduction(g,b):
 steps=[];a=b
 while sum(a)>1:
  pairing=[2*a[i]-sum(a[j] for j in g.adj[i]) for i in range(g.n)]
  positives=[i for i in range(g.n) if pairing[i]>0]
  if not positives:return None
  s=positives[0];v=list(a);v[s]-=pairing[s]
  if min(v)<0:return None
  a=tuple(v);steps.append(s)
 if a not in g.e:return None
 return dict(conjugator=steps,simple_root=g.e.index(a))

def search(n,max_entry,max_only):
 g=en(n);start=time.monotonic();inv=inverse_cartan(g)
 denom=1
 from math import lcm
 for row in inv:
  for v in row:denom=lcm(denom,v.denominator)
 N=[[int(v*denom) for v in row] for row in inv]
 Is=independent_sets(g);alpha=max(map(len,Is));Is=[I for I in Is if I and (not max_only or len(I)==alpha)]
 checked=0;positive_integer=0;norm_two=0;out=[]
 for I in Is:
  for pos in product(range(1,max_entry+1),repeat=len(I)):
   m=[None]*n
   for s,v in zip(I,pos):m[s]=v
   rest=[j for j in range(n) if j not in I];options=[]
   for j in rest:
    demand=max([m[s] for s in g.adj[j] if s in I] or [0])
    options.append(tuple(range(-max_entry,-demand+1)))
   for values in product(*options):
    checked+=1
    for j,v in zip(rest,values):m[j]=v
    bs=[sum(x*y for x,y in zip(row,m)) for row in N]
    if any(x<=0 or x%denom for x in bs):continue
    positive_integer+=1;b=tuple(x//denom for x in bs)
    if sum(x*y for x,y in zip(b,m))!=2:continue
    norm_two+=1;cert=root_reduction(g,b)
    if cert is None:continue
    w=tuple(tuple(int(i==j)-m[j]*b[i] for i in range(n)) for j in range(n))
    assert g.terminal(w) and tuple(g.desc(w))==I and g.inv(w)==w
    c=cert['conjugator'];word=tuple(c)+(cert['simple_root'],)+tuple(reversed(c));assert g.elt(word)==w
    out.append(dict(rank=n,beta=b,pairings=tuple(m),I=I,maximum_independent=len(I)==alpha,independence_number=alpha,height=sum(b),length=g.length(w),root_reduction=cert,reduced_word=g.word(w)))
    print(json.dumps(out[-1]),flush=True)
 return dict(rank=n,max_pairing_absolute_value=max_entry,maximum_independent_only=max_only,independence_number=alpha,independent_sets_tested=len(Is),pairings_tested=checked,positive_integer_vectors=positive_integer,norm_two_vectors=norm_two,real_terminal_roots=out,seconds=time.monotonic()-start)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('rank',type=int);p.add_argument('--max-entry',type=int,default=1);p.add_argument('--max-only',action='store_true');a=p.parse_args()
 data=search(a.rank,a.max_entry,a.max_only)
 Path(__file__).with_name(f'cartan_E{a.rank}_m{a.max_entry}_max{int(a.max_only)}.json').write_text(json.dumps(data,indent=2)+'\n')
 print('SUMMARY', {k:v for k,v in data.items() if k!='real_terminal_roots'},'roots',len(data['real_terminal_roots']),flush=True)
