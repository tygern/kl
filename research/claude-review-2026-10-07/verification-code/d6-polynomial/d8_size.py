import sys, time
sys.path.insert(0, ".")
from kl_indep import *
D8 = SignedD(8)
w8 = (1, -8, 3, -6, 5, -4, 7, -2)
x8 = word_product(D8, [0,1,3,5,7])
print("x8", x8, "len", D8.length(x8), "w8 len", D8.length(w8), "even signs", sum(a<0 for a in w8)%2==0)
def rw(G, w):
    word=[]
    while G.length(w)>0:
        s=next(s for s in G.gens if G.is_descent(w,s)); word.append(s); w=G.right(w,s)
    return word[::-1]
word = rw(D8, w8)
print("reduced word (Gern labels)", [s+1 for s in word])
t=time.time()
S={D8.identity()}
for s in word:
    S |= {D8.right(x,s) for x in S}
print("ideal size", len(S), "time", round(time.time()-t,1))
from collections import Counter
c=Counter(D8.length(w) for w in S); print(sorted(c.items()))
