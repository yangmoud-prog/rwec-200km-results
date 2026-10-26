from pathlib import Path
import json, hashlib, sys
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')

P=Path(__file__).resolve().parent
W=P.parents[3]
A=P.parent
P8=A/'08_forward_C_moving_endpoint_20261026'
B=W/'variants/mild_control_slope075_20pilot_v1'
sys.path.insert(0,str(P8))
GRID=np.arange(289,dtype=float)
TARGETS=np.array([120,100,80,60,40,20,0])
for name in ['data','tables','figures','provenance','predictions']:(P/name).mkdir(exist_ok=True)
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def frame(name,rows):pd.DataFrame(rows).to_csv(P/'tables'/name,index=False,encoding='utf-8-sig')
def load(p):
 with np.load(p) as d:return {k:d[k] for k in d.files}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def extent(h,x):
 mask=h>=.1
 if not mask.any():return 0.,0.,0,False
 idx=np.where(mask)[0];j=idx[-1]
 end=x[j] if j==len(x)-1 else x[j]+(.1-h[j])*(x[j+1]-x[j])/(h[j+1]-h[j])
 connected=0.
 if mask[0]:
  n=np.where(~mask)[0]
  k=n[0]-1 if len(n) else len(x)-1
  connected=x[k] if k==len(x)-1 else x[k]+(.1-h[k])*(x[k+1]-x[k])/(h[k+1]-h[k])
 return end,connected,int(mask[0])+int(np.sum(np.diff(mask.astype(int))==1)),bool(mask[-1])
def metrics(y,p):
 if not np.all(np.isfinite(p)):return dict(NSE=np.nan,RMSE=np.nan,KGE=np.nan,peak_pct=np.nan,peak_time_h=np.nan)
 r=np.corrcoef(y,p)[0,1];sd=np.std(y);mu=np.mean(y)
 return dict(NSE=1-np.sum((p-y)**2)/np.sum((y-y.mean())**2),RMSE=np.sqrt(np.mean((p-y)**2)),KGE=1-np.sqrt((r-1)**2+(np.std(p)/sd-1)**2+(p.mean()/mu-1)**2),peak_pct=100*(p.max()-y.max())/y.max(),peak_time_h=float(np.argmax(p)-np.argmax(y)))
