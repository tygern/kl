import sys
n=int(sys.argv[1]); word=[int(c) for c in sys.argv[2]]
adj=[set() for _ in range(n)]
def edge(a,b): adj[a].add(b); adj[b].add(a)
for i in range(n-2): edge(i,i+1)
edge(n-1,2)
nb=[tuple(sorted(adj[i])) for i in range(n)]
Id=tuple(tuple(1 if a==b else 0 for a in range(n)) for b in range(n))
def rmul(M,i):
    cols=list(M); ai=cols[i]; cols[i]=tuple(-a for a in ai)
    for j in nb[i]: cols[j]=tuple(a+b for a,b in zip(cols[j],ai))
    return tuple(cols)
S={Id}
for s in word:
    S|={rmul(M,s) for M in S}
# length distribution of ideal: count by number of negative columns? no; just size
print(f"E{n} word length {len(word)}: Bruhat ideal size = {len(S)}")
