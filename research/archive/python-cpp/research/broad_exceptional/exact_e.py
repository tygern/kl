"""Exact simply laced finite Coxeter groups; integral simple-root images.
E_r diagram: 0-1-2-3-...-(r-2), with node r-1 attached to node 2.
No floating point; elements faithfully act on the root lattice.
"""
from functools import cache
from collections import Counter
from time import perf_counter
import json

class Weyl:
    def __init__(self, rank=6, edges=None):
        self.rank=rank
        self.edges=edges if edges is not None else [(i,i+1) for i in range(rank-2)]+[(2,rank-1)]
        self.adj=[[] for _ in range(rank)]
        for s,t in self.edges:self.adj[s].append(t);self.adj[t].append(s)
        basis=[tuple(int(i==j) for i in range(rank)) for j in range(rank)]
        roots=set(basis);todo=list(basis)
        for a in todo:
            for s in range(rank):
                b=list(a);b[s]=-a[s]+sum(a[t] for t in self.adj[s]);b=tuple(b)
                if b not in roots:roots.add(b);todo.append(b)
        self.roots=sorted(roots);self.ri={a:i for i,a in enumerate(self.roots)}
        self.neg=[self.ri[tuple(-v for v in a)] for a in self.roots]
        self.positive=[all(v>=0 for v in a) for a in self.roots]
        self.add={}
        for i,a in enumerate(self.roots):
            for j,b in enumerate(self.roots):
                c=tuple(x+y for x,y in zip(a,b))
                if c in self.ri:self.add[i,j]=self.ri[c]
        self.sperm=[]
        for s in range(rank):
            row=[]
            for a in self.roots:
                b=list(a);b[s]=-a[s]+sum(a[t] for t in self.adj[s]);row.append(self.ri[tuple(b)])
            self.sperm.append(row)
        e=tuple(self.ri[a] for a in basis)
        self.elements=[e];self.index={e:0};self.words=[()];self.length=[0];self.right=[]
        for i,w in enumerate(self.elements):
            row=[]
            for s in range(rank):
                v=self.mult(w,s)
                if v not in self.index:
                    self.index[v]=len(self.elements);self.elements.append(v)
                    self.length.append(self.length[i]+1);self.words.append(self.words[i]+(s,))
                row.append(self.index[v])
            self.right.append(tuple(row))
        self.desc=[tuple(s for s in range(rank) if not self.positive[w[s]]) for w in self.elements]
        self.inv=[]
        for word in self.words:
            j=0
            for s in reversed(word):j=self.right[j][s]
            self.inv.append(j)
        self.ldesc=[self.desc[self.inv[i]] for i in range(len(self.elements))]
        self.left=[tuple(self.inv[self.right[self.inv[i]][s]] for s in range(rank)) for i in range(len(self.elements))]
        self.fc=[]
        for i in range(len(self.elements)):
            self.fc.append(all(self.fc[self.right[i][s]] for s in self.desc[i]) and all(not (s in self.desc[i] and t in self.desc[i]) for s,t in self.edges))
    def mult(self,w,s):
        v=list(w);v[s]=self.neg[w[s]]
        for t in self.adj[s]:v[t]=self.add[w[s],w[t]]
        return tuple(v)
    def word(self,word):
        w=0
        for s in word:w=self.right[w][s]
        return w
    @cache
    def leq(self,x,w):
        while x!=w:
            if self.length[x]>=self.length[w]:return False
            s=self.desc[w][0]
            if s in self.desc[x]:x=self.right[x][s]
            w=self.right[w][s]
        return True
    @cache
    def lower(self,w):
        if w==0:return frozenset([0])
        s=self.desc[w][0];v=self.right[w][s];base=self.lower(v)
        return base|{self.right[x][s] for x in base}
    @cache
    def down(self,w):
        word=self.words[w];res=set()
        for k in range(len(word)):
            z=self.word(word[:k]+word[k+1:])
            if self.length[z]==self.length[w]-1:res.add(z)
        return tuple(sorted(res))
    def interval(self,x,w):
        return frozenset(z for z in self.lower(w) if self.leq(x,z))
    def weaklybad(self,w):
        for s in self.desc[w]:
            if any(t in self.desc[self.right[w][s]] for t in self.adj[s]):return False
        for s in self.ldesc[w]:
            if any(t in self.ldesc[self.left[w][s]] for t in self.adj[s]):return False
        return True
    def hardpairs(self):
        fcs=[i for i in range(len(self.elements)) if self.fc[i]]
        for w in range(len(self.elements)):
            if self.weaklybad(w) and not self.fc[w]:
                xs=[x for x in fcs if set(self.desc[w])<=set(self.desc[x]) and set(self.ldesc[w])<=set(self.ldesc[x]) and self.leq(x,w)]
                yield w,xs
    def stats(self):
        result=[]
        for w,xs in self.hardpairs():
            result.append(dict(w=w,word=self.words[w],length=self.length[w],R=self.desc[w],L=self.ldesc[w],full_support=len(set(self.words[w]))==self.rank,
                 xs=[dict(x=x,word=self.words[x],length=self.length[x],rank=self.length[w]-self.length[x]) for x in xs]))
        return dict(order=len(self.elements),roots=len(self.roots),fc=sum(self.fc),max_length=max(self.length),length_distribution=dict(Counter(self.length)),bad=result)

if __name__=='__main__':
    start=perf_counter();g=Weyl();print(json.dumps(g.stats(),indent=2));print('seconds',perf_counter()-start)
