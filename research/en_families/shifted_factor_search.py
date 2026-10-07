from targeted import en
from find_a_factors import search
from pathlib import Path
from collections import deque
import json,time
G=en(9);data=json.load(open('research/en_families/affine_reflection_family.json'));beta=data['beta'];delta=data['delta'];m=data['cartan_pairing']
def reflection(k):
 b=tuple(v+k*d for v,d in zip(beta,delta))
 return tuple(tuple(int(i==j)-m[j]*b[i] for i in range(9)) for j in range(9))
def weakcount(w):
 seen={w};todo=deque([w])
 while todo:
  z=todo.popleft()
  for s in G.desc(z):
   y=G.right(z,s)
   if y not in seen:seen.add(y);todo.append(y)
 return len(seen)
out=[]
for k in (1,2):
 w=reflection(k);print('WEAK',k,weakcount(w),flush=True)
 for extra,threshold in [(None,8),(0,7),(1,7)]:
  top=w if extra is None else G.right(w,extra);t0=time.monotonic();res=search(G,top,threshold)
  res.update(k=k,extra_generator=extra,threshold=threshold,length=G.length(top),seconds=time.monotonic()-t0)
  out.append(res);print(json.dumps(res),flush=True)
  G.right.cache_clear();G.desc.cache_clear();G.word.cache_clear()
Path(__file__).with_name('shifted_factor_search.json').write_text(json.dumps(out,indent=2)+'\n')
