"""Exact bounded simply-laced Coxeter computations using integral matrices.
No root-system or group enumeration: only words and their finite lower ideals.
"""
from functools import cache
import json

def add(a,b,shift=0,scale=1):
    c=list(a)+[0]*max(0,len(b)+shift-len(a))
    for i,v in enumerate(b):c[i+shift]+=scale*v
    while c and c[-1]==0:c.pop()
    return tuple(c)

class Coxeter:
    def __init__(self,n,edges):
        self.n=n;self.edges=edges
        self.kl_values={}
        self.adj=[[] for _ in range(n)]
        for s,t in edges:self.adj[s].append(t);self.adj[t].append(s)
        self.e=tuple(tuple(int(i==j) for i in range(n)) for j in range(n))
    @cache
    def right(self,w,s):
        out=list(w);out[s]=tuple(-a for a in w[s])
        for t in self.adj[s]:out[t]=tuple(a+b for a,b in zip(w[s],w[t]))
        return tuple(out)
    @cache
    def desc(self,w):
        return tuple(i for i,a in enumerate(w) if all(x<=0 for x in a))
    @cache
    def word(self,w):
        if w==self.e:return ()
        ds=self.desc(w);assert ds
        s=ds[0]
        return self.word(self.right(w,s))+(s,)
    def length(self,w):return len(self.word(w))
    def elt(self,word):
        w=self.e
        for s in word:w=self.right(w,s)
        return w
    @cache
    def inv(self,w):return self.elt(reversed(self.word(w)))
    @cache
    def fc(self,w):
        ds=self.desc(w)
        return all(t not in ds for s in ds for t in self.adj[s]) and all(self.fc(self.right(w,s)) for s in ds)
    def terminal(self,w):
        for v in (w,self.inv(w)):
            for s in self.desc(v):
                if any(t in self.desc(self.right(v,s)) for t in self.adj[s]):return False
        return True
    @cache
    def lower(self,w):
        if w==self.e:return frozenset([w])
        s=self.desc(w)[0];base=self.lower(self.right(w,s))
        return base|{self.right(z,s) for z in base}
    @cache
    def leq(self,x,w):
        while x!=w:
            if self.length(x)>=self.length(w):return False
            s=self.desc(w)[0]
            if s in self.desc(x):x=self.right(x,s)
            w=self.right(w,s)
        return True
    def eligible(self,w):
        rd=set(self.desc(w));ld=set(self.desc(self.inv(w)))
        return [x for x in self.lower(w) if self.fc(x) and rd<=set(self.desc(x)) and ld<=set(self.desc(self.inv(x)))]
    @cache
    def corrections(self,w,s):
        v=self.right(w,s);out=[]
        for z in self.lower(v):
            d=self.length(v)-self.length(z)
            if d<=0 or d%2==0 or s not in self.desc(z):continue
            p=self.kl(z,v);mu=p[(d-1)//2] if len(p)>(d-1)//2 else 0
            if mu:out.append((z,(d+1)//2,mu))
        return tuple(out)
    @cache
    def kl(self,x,w):
        if x==w:p=(1,)
        elif not self.leq(x,w):p=()
        else:p=self.with_descent(x,w,self.desc(w)[0])
        self.kl_values[x,w]=p
        return p
    def with_descent(self,x,w,s):
        v=self.right(w,s);xs=self.right(x,s);c=int(s in self.desc(x))
        p=add((),self.kl(xs,v),1-c);p=add(p,self.kl(x,v),c)
        for z,d,mu in self.corrections(w,s):
            if self.leq(x,z):p=add(p,self.kl(x,z),d,-mu)
        return p
    def info(self,word,lower=False,kl=False):
        w=self.elt(word);out=dict(word=word,length=self.length(w),R=self.desc(w),L=self.desc(self.inv(w)),terminal=self.terminal(w),fc=self.fc(w))
        if lower:
            out['lower_size']=len(self.lower(w));xs=self.eligible(w)
            out['eligible']=[dict(word=self.word(x),length=self.length(x),rank=self.length(w)-self.length(x),polynomial=self.kl(x,w) if kl else None) for x in xs]
        return out

def star(m):return Coxeter(m+1,[(0,i) for i in range(1,m+1)])
def en(n):return Coxeter(n,[(i,i+1) for i in range(n-2)]+[(2,n-1)])
def affine_d(n):return Coxeter(n+1,[(i,i+1) for i in range(2,n-2)]+[(0,2),(1,2),(n-2,n-1),(n-2,n)])
def bipartition(g):
    color={0:0};todo=[0]
    for s in todo:
        for t in g.adj[s]:
            if t not in color:color[t]=1-color[s];todo.append(t)
            else:assert color[t]!=color[s]
    return tuple(tuple(s for s in range(g.n) if color[s]==c) for c in (0,1))

if __name__=='__main__':
    import argparse,time
    a=argparse.ArgumentParser();a.add_argument('--star',type=int);a.add_argument('--en',type=int);a.add_argument('--affine-d',type=int);a.add_argument('--layers',type=int,default=3);a.add_argument('--color',type=int,default=1);a.add_argument('--lower',action='store_true');a.add_argument('--kl',action='store_true');args=a.parse_args()
    g=star(args.star) if args.star else en(args.en) if args.en else affine_d(args.affine_d)
    colors=bipartition(g);word=sum((colors[(args.color+i)%2] for i in range(args.layers)),())
    t=time.perf_counter();out=g.info(word,args.lower or args.kl,args.kl);out['seconds']=time.perf_counter()-t
    print(json.dumps(out,indent=2))
