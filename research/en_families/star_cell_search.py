"""Bounded left/right length-three-star orbit search for an a-value witness."""
from targeted import en
from conjugate_search import left
from find_a_factors import longest
from collections import deque
from pathlib import Path
import json,time,argparse
p=argparse.ArgumentParser();p.add_argument('--cap',type=int,default=30000);p.add_argument('--length',type=int,default=34);a=p.parse_args()
g=en(9);rows=json.load(open('research/en_families/e9_full_support_eligible.json'));row=next(r for r in rows if r['length']==a.length);start=g.elt(row['word'])
parents={start:None};todo=deque([(start,g.inv(start))]);weights={};found=None;t0=time.monotonic()
while todo and len(parents)<a.cap:
 w,wi=todo.popleft()
 for hand,v in [('R',w),('L',wi)]:
  ds=g.desc(v)
  if ds not in weights:weights[ds]=g.length(longest(g,ds))
  if weights[ds]>=7:
   found=(w,hand,ds,weights[ds]);break
 if found:break
 for hand,v,vi in [('R',w,wi),('L',wi,w)]:
  ds=g.desc(v)
  for s,t in g.edges:
   if (s in ds)==(t in ds):continue
   z=g.right(v,s);dz=g.desc(z);letter=s
   if (s in dz)==(t in dz):z=g.right(v,t);letter=t
   dz=g.desc(z);assert (s in dz)!=(t in dz)
   zinv=left(g,vi,letter)
   new,newi=(z,zinv) if hand=='R' else (zinv,z)
   if new not in parents:
    parents[new]=(w,hand,s,t,letter);todo.append((new,newi))
steps=[]
if found:
 w,hand,ds,weight=found;result=dict(found=True,word=g.word(w),descent_hand=hand,descents=ds,factor_length=weight)
 while parents[w] is not None:
  prev,hand,s,t,letter=parents[w];steps.append(dict(hand=hand,pair=[s,t],letter=letter));w=prev
 result['star_steps']=list(reversed(steps))
else:result=dict(found=False)
result.update(start_word=row['word'],states_seen=len(parents),queue_remaining=len(todo),state_cap=a.cap,seconds=time.monotonic()-t0)
Path(__file__).with_name(f'star_cell_E9_length{a.length}.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
