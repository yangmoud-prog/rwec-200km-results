from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
OUT=Path(__file__).resolve().parent
A=OUT.parents[1];WORK=A.parents[1]
P4=A/'analysis/04_autonomous_backwater_20261021';P5=A/'analysis/05_routing_validation_slope075_20261022'
B=WORK/'variants/mild_control_slope075_20pilot_v1'
TARGETS=np.array([20,40,60,80,100,120])
GRID=np.arange(-48.,241.)
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,v):
 p=Path(p);p.parent.mkdir(exist_ok=True,parents=True);p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def load(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def frame(name,r):pd.DataFrame(r).to_csv(OUT/'tables'/name,index=False,encoding='utf-8-sig')
def dirs():
 for d in ['data','tables','figures','fits','predictions','provenance','logs']:(OUT/d).mkdir(exist_ok=True)
def metric(y,p):
 y=np.asarray(y);p=np.asarray(p);d=p-y;den=np.sum((y-y.mean())**2)
 corr=np.corrcoef(y,p)[0,1] if y.std()>1e-12 and p.std()>1e-12 else np.nan
 kge=1-np.sqrt((corr-1)**2+(p.std()/y.std()-1)**2+(p.mean()/y.mean()-1)**2) if y.std()>1e-12 and abs(y.mean())>1e-12 else np.nan
 return dict(NSE=float(1-np.sum(d*d)/den) if den>1e-12 else np.nan,RMSE=float(np.sqrt(np.mean(d*d))),KGE=float(kge),bias=float(d.mean()))
def extent(v,x,threshold=.1):
 ids=np.flatnonzero(v>threshold)
 if not len(ids):return 0.,False,0,0.
 j=ids[-1];far=x[j] if j==len(x)-1 else x[j]+(x[j+1]-x[j])*(v[j]-threshold)/(v[j]-v[j+1])
 active=v>threshold;islands=int(np.sum(active&~np.r_[False,active[:-1]]));connected=0.
 if active[0]:
  z=np.flatnonzero(~active)
  if not len(z):connected=x[-1]
  else:
   k=z[0]-1;connected=x[k]+(x[k+1]-x[k])*(v[k]-threshold)/(v[k]-v[k+1])
 return float(far),j==len(x)-1,islands,float(connected)
