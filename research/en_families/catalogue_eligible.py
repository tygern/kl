"""Independently check matrix/FC/Bruhat predicates on filtered complete FC lists.
Catalogue completeness is supplied by parent's separate enumerator.
"""
from targeted import en
from pathlib import Path
import json,argparse
p=argparse.ArgumentParser();p.add_argument('rank',type=int);a=p.parse_args();n=a.rank
rows=json.loads(Path(__file__).with_name(f'conjugates_E{n}_50000_934.json').read_text())['terminals']
catalogue=json.load(open(f'research/en_independent/e{n}-fc.json'))['elements'];g=en(n);out=[]
for row in rows:
 if not row['full_support']:continue
 w=g.elt(row['word']);assert g.terminal(w)
 rm=sum(1<<s for s in g.desc(w));lm=sum(1<<s for s in g.desc(g.inv(w)))
 matches=[x for x in catalogue if x['Rmask']&rm==rm and x['Lmask']&lm==lm]
 eligible=[]
 for item in matches:
  x=g.elt(item['word']);assert g.fc(x)
  assert sum(1<<s for s in g.desc(x))==item['Rmask']
  assert sum(1<<s for s in g.desc(g.inv(x)))==item['Lmask']
  if g.leq(x,w):eligible.append(dict(word=item['word'],length=g.length(x),rank=g.length(w)-g.length(x)))
 out.append(dict(**row,reflection_rank=(n-sum(w[i][i] for i in range(n)))//2,mask_candidates=len(matches),eligible=eligible))
 print(out[-1],flush=True)
Path(__file__).with_name(f'e{n}_full_support_eligible.json').write_text(json.dumps(out,indent=2)+'\n')
