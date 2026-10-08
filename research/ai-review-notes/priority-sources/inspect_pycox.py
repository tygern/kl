"""Bounded static data inspection; does not import or execute PyCox."""
import ast, hashlib, json
from pathlib import Path
p=Path(__file__).with_name('chv1r6180.py')
tree=ast.parse(p.read_text())
def literal(n):
    if isinstance(n,ast.BinOp) and isinstance(n.op,ast.Add): return literal(n.left)+literal(n.right)
    if isinstance(n,ast.List): return [literal(x) for x in n.elts]
    if isinstance(n,ast.Dict): return {literal(k):literal(v) for k,v in zip(n.keys,n.values)}
    return ast.literal_eval(n)
out={}
for n in tree.body:
    if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ('E7KLCELLREPS','E8KLCELLREPS') for t in n.targets):
        rows=literal(n.value)
        out[n.targets[0].id]={'rows':len(rows),'keys':sorted(set().union(*(r.keys() for r in rows))), 'total_seed_words':sum(len(r['replstar']) for r in rows),'sum_cell_sizes':sum(r['size'] for r in rows),'initial_elms_empty':all(r['elms']==[] for r in rows),'source_line':n.lineno,'source_end_line':n.end_lineno,'first_record':rows[0]}
out['source_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
meta=json.loads(p.with_name('pycox-commit.json').read_text())
out['commit']=meta['sha'];out['commit_date']=meta['commit']['committer']['date']
print(json.dumps(out,indent=2))
