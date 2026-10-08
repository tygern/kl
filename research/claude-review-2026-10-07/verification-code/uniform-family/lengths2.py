exec(open("lengths.py").read().split("for r in (3, 4):")[0])
for r,ks in ((3,(4,5)),(4,(2,3)),(5,(0,1,2)),(6,(0,1))):
    n=4*r+1; A=cartan(n); beta=seed(r); g,d=gamma_delta(r)
    gd=[x+y for x,y in zip(g,d)]
    v=list(beta)
    for k in range(max(ks)+1):
        if k in ks:
            L=len(reduced_word(A, matrix_of_reflection(A,v)))
            print(f"r={r} k={k}: length {L}; 56r-93+116k = {56*r-93+116*k}")
        v=refl_root(A,g,refl_root(A,gd,v))
