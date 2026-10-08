import importlib.util, sys
spec = importlib.util.spec_from_file_location("cu", "check_uniform.py")
# reuse functions without rerunning: copy minimal
exec(open("check_uniform.py").read().split("results = []")[0])
for r in (3,4,5,6):
    n=4*r+1; A=cartan(n); beta=seed(r); g,d=gamma_delta(r); a=2*r-4
    I_r = sorted({0,4*r}|set(range(3,4*r,2)))
    gd=[x+y for x,y in zip(g,d)]
    v=list(beta)
    for k in range(6):
        cols=[refl_root(A,v,simple(n,j)) for j in range(n)]
        R=sorted(j for j in range(n) if any(c<0 for c in cols[j]))
        mixed=[j for j in range(n) if any(c<0 for c in cols[j]) and any(c>0 for c in cols[j])]
        zeros_in_desc=[(j,[i for i in range(n) if cols[j][i]==0]) for j in I_r if 0 in cols[j]]
        print(f"r={r} k={k}: R(b)==I_r: {R==I_r}; mixed-sign columns: {mixed}; zero entries in descent columns: {zeros_in_desc}")
        v=refl_root(A,g,refl_root(A,gd,v))
