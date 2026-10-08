import sys
dump, brute = sys.argv[1], sys.argv[2]
B = {}
for line in open(brute):
    x, y, p = line.rstrip('\n').split('\t'); B[(x, y)] = p
D = {}
for line in open(dump):
    x, y, p = line.rstrip('\n').split('\t'); D[(x, y)] = p
missing = [k for k in D if k not in B]
diff = [(k, D[k], B[k]) for k in D if k in B and D[k] != B[k]]
print(f"{dump}: {len(D)} extremal pairs from engine; brute has {len(B)} comparable pairs; missing in brute: {len(missing)}; value mismatches: {len(diff)}")
for d in diff[:10]: print(d)
# also: all brute pairs with x extremal for y should be in D -> check via counting polys with mu=2 etc
from collections import Counter
print("brute polynomial census:", sorted(Counter(B.values()).items(), key=lambda t: -t[1])[:12])
