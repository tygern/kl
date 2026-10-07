"""Search reduced factors equal to finite-parabolic longest elements.
All weak prefixes of each requested upper are visited unless a witness is found.
"""
from targeted import en
from collections import deque
from itertools import combinations
from pathlib import Path
import json,time

def multiply(g,w,word):
 for s in word:w=g.right(w,s)
 return w

def longest(g,J):
 w=g.e
 while True:
  ds=g.desc(w);ss=next((s for s in J if s not in ds),None)
  if ss is None:return w
  w=g.right(w,ss)

def search(g,w,threshold=7):
 n=g.n;paras=[]
 for size in range(1,n):
  for J in combinations(range(n),size):
   v=longest(g,J);lv=g.length(v)
   if lv<threshold:continue
   if any(set(K)<=set(J) for K,u,l,word in paras):continue
   paras.append((J,v,lv,g.word(v)))
 todo=deque([(w,())]);seen={w};tested=0
 while todo:
  p,suffix=todo.popleft();lp=g.length(p)
  for J,u,lu,wu in paras:
   if lu>lp:continue
   # A suffix u starts with one of its simple descents; current p must have
   # a right descent in J before trying full multiplication.
   if not set(J)&set(g.desc(p)):continue
   a=multiply(g,p,reversed(wu));tested+=1
   if g.length(a)+lu==lp:
    return dict(found=True,J=J,factor_length=lu,prefix=g.word(a),factor=wu,suffix=suffix,weak_prefixes=len(seen),tests=tested)
  for s in g.desc(p):
   q=g.right(p,s)
   if q not in seen:seen.add(q);todo.append((q,(s,)+suffix))
 return dict(found=False,weak_prefixes=len(seen),tests=tested,minimal_parabolics=len(paras))

if __name__=='__main__':
 g=en(9);data=json.load(open('research/en_families/e9_full_support_eligible.json'));out=[]
 rows=[row for row in data if row['length'] in (32,34)]
 w33=next(row['word'] for row in data if row['length']==33);rows.append(dict(word=w33+[0],length=34))
 for row in rows:
  start=time.monotonic();result=search(g,g.elt(row['word']));result.update(word=row['word'],length=row['length'],seconds=time.monotonic()-start)
  out.append(result);print(json.dumps(result),flush=True)
  g.right.cache_clear();g.desc.cache_clear();g.word.cache_clear()
 Path(__file__).with_name('parabolic_factor_search.json').write_text(json.dumps(out,indent=2)+'\n')
