"""Exact affine-D ball and rank-three crown search; standard library only.
Simple 0 flips/swaps first pair; i=1..n-1 swaps i-1,i;
affine n flips/swaps last pair and adds e_(n-1)+e_n.
All elements are faithful affine signed permutation actions.
"""
from collections import defaultdict,Counter
import json,sys

class AffineD:
 def __init__(self,n,L):
  self.n=n;self.e=(tuple(range(1,n+1)),(0,)*n)
  self.els=[self.e];self.ids={self.e:0};self.words=[()];self.lengths=[0]
  for w in self.els:
   wi=self.ids[w]
   if self.lengths[wi]==L:continue
   for s in range(n+1):
    z=self.simple(w,s)
    if z not in self.ids:
     self.ids[z]=len(self.els);self.els.append(z);self.words.append(self.words[wi]+(s,));self.lengths.append(self.lengths[wi]+1)
  self.down=[]
  for wi,Q in enumerate(self.words):
   self.down.append(tuple(set(self.ids[self.word(Q[:j]+Q[j+1:])] for j in range(len(Q)) if self.lengths[self.ids[self.word(Q[:j]+Q[j+1:])]]==len(Q)-1)))
 def simple(self,w,s):
  p,t=map(list,w)
  if s==0:p[0],p[1]=-p[1],-p[0]
  elif s<self.n:p[s-1],p[s]=p[s],p[s-1]
  else:
   for a in p[-2:]:t[abs(a)-1]+=1 if a>0 else -1
   p[-2],p[-1]=-p[-1],-p[-2]
  return tuple(p),tuple(t)
 def word(self,Q):
  w=self.e
  for s in Q:w=self.simple(w,s)
  return w
 def crowns(self):
  maximum=0;examples={};counts=Counter()
  for wi in range(len(self.els)):
   if self.lengths[wi]<3:continue
   coatoms=defaultdict(set)
   for a in self.down[wi]:
    for b in self.down[a]:
     for x in self.down[b]:coatoms[x].add(a)
   for x,cs in coatoms.items():
    counts[len(cs)]+=1
    if len(cs) not in examples:
     examples[len(cs)]={'bottom':self.words[x],'top':self.words[wi],'bottom_length':self.lengths[x],'top_length':self.lengths[wi],'coatoms':[self.words[c] for c in cs]}
    if len(cs)>maximum:
     maximum=len(cs);print('new_maximum',maximum,examples[maximum],flush=True)
  return {'group':f'affineD{self.n}','ball_size':len(self.els),'layer_sizes':dict(Counter(self.lengths)),'rank3_crown_histogram':dict(counts),'examples':examples}
if __name__=='__main__':
 n=int(sys.argv[1]) if len(sys.argv)>1 else 4
 L=int(sys.argv[2]) if len(sys.argv)>2 else 10
 W=AffineD(n,L)
 result=W.crowns()
 with open(f'research/broad_affine/crowns_D{n}_L{L}.json','w') as f:json.dump(result,f,indent=2)
 print(json.dumps(result,indent=2))
