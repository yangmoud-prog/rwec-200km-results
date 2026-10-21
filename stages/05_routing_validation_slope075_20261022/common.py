from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd
OUT=Path(__file__).resolve().parent
A=OUT.parents[1]
WORK=A.parents[1]
P4=A/'analysis/04_autonomous_backwater_20261021'
B=WORK/'variants/mild_control_slope075_20pilot_v1'
TARGETS=np.array([20,40,60,80,100,120])
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def now():return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(timespec='seconds')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def loadnp(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def frame(name,rows):pd.DataFrame(rows).to_csv(OUT/'tables'/name,index=False,encoding='utf-8-sig')
def dirs():
    for d in ['tables','data','predictions','fits','figures','provenance','logs']:(OUT/d).mkdir(exist_ok=True)
def extent(dz,x,threshold=.1):
    ids=np.flatnonzero(dz>threshold)
    if not len(ids):return 0.
    j=ids[-1]
    if j==len(x)-1:return float(x[j])
    return float(x[j]+(x[j+1]-x[j])*(dz[j]-threshold)/(dz[j]-dz[j+1]))
def metrics(o,p,t):
    err=p-o;den=np.sum((o-o.mean())**2);r=np.corrcoef(o,p)[0,1] if p.std()>0 else 0.
    return dict(NSE=float(1-np.sum(err*err)/den),RMSE_m3s=float(np.sqrt(np.mean(err*err))),
      KGE=float(1-np.sqrt((r-1)**2+(p.std()/o.std()-1)**2+(p.mean()/o.mean()-1)**2)),
      peak_error_pct=float(100*(p.max()/o.max()-1)),peak_time_error_h=float(t[np.argmax(p)]-t[np.argmax(o)]))
