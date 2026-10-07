"""Certify degree-four obstruction to principal finite-D unfolding."""
from affine_d_crowns import AffineD
from collections import Counter
import json
W=AffineD(4,10)
Q=(0,1,3,4,2)*2
w=W.ids[W.word(Q)]
assert W.lengths[w]==10
nodes={w};stack=[w]
while stack:
 z=stack.pop()
 for a in W.down[z]:
  if a not in nodes:nodes.add(a);stack.append(a)
ranks=Counter(W.lengths[z] for z in nodes)
edges=set()
for z in nodes:
 if W.lengths[z]==3:
  support=set(W.words[z])
  if len(support)==2:edges.add(tuple(sorted(support)))
assert edges=={(0,2),(1,2),(2,3),(2,4)}
# Coxeter relations in the faithful affine model.
for s in range(5):
 assert W.word((s,s))==W.e
 for t in range(s):
  m=3 if 2 in (s,t) else 2
  assert W.word((s,t)*m)==W.e
result={'top_word':Q,'length':10,'interval_size':len(nodes),'rank_vector':[ranks[i] for i in range(11)],'detected_rank3_edges':sorted(edges),'maximum_detected_degree':4}
with open('research/broad_affine/unfolding_obstruction.json','w') as f:json.dump(result,f,indent=2)
print(json.dumps(result,indent=2))
