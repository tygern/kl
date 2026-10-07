from targeted import en
from pathlib import Path
import json,argparse,time
p=argparse.ArgumentParser();p.add_argument('rank',type=int);p.add_argument('length',type=int);p.add_argument('--kl',action='store_true');a=p.parse_args()
rows=json.loads(Path(__file__).with_name(f'conjugates_E{a.rank}_50000_934.json').read_text())['terminals']
rows=[row for row in rows if row['full_support'] and row['length']==a.length]
out=[]
for row in rows:
 g=en(a.rank);start=time.monotonic();data=g.info(row['word'],lower=True,kl=a.kl);data['seconds']=time.monotonic()-start
 data['reflection_rank']=(a.rank-sum(g.elt(row['word'])[i][i] for i in range(a.rank)))//2
 out.append(data);print(json.dumps(data),flush=True)
Path(__file__).with_name(f'E{a.rank}_length{a.length}_eligible.json').write_text(json.dumps(out,indent=2)+'\n')
