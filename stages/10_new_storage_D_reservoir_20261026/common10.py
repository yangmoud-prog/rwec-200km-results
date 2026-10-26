from pathlib import Path
import sys,json,hashlib
import numpy as np
import pandas as pd
if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parent
W=P.parents[3]
P8=P.parent/'08_forward_C_moving_endpoint_20261026'
P9=P.parent/'09_joint_state_coupling_20261026'
D=W/'variants/compound_reservoir_D_20261026'
TARGETS=np.array([120,100,80,60,40,20,0])
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def load(p):
 with np.load(p) as d:return {k:d[k] for k in d.files}
def data(domain,e):return load((P9/'data'/f'B_{e}.npz') if domain=='B' else (P8/'data'/f'{e}.npz'))
def frame(name,rows):pd.DataFrame(rows).to_csv(P/'tables'/name,index=False,encoding='utf-8-sig')
def metrics(y,p):
 if not np.all(np.isfinite(p)):return dict(NSE=np.nan,RMSE=np.nan,KGE=np.nan,peak_pct=np.nan,peak_time_h=np.nan)
 r=np.corrcoef(y,p)[0,1];sd=np.std(y);mu=np.mean(y)
 return dict(NSE=1-np.sum((p-y)**2)/np.sum((y-y.mean())**2),RMSE=np.sqrt(np.mean((p-y)**2)),KGE=1-np.sqrt((r-1)**2+(np.std(p)/sd-1)**2+(p.mean()/mu-1)**2),peak_pct=100*(p.max()-y.max())/y.max(),peak_time_h=float(np.argmax(p)-np.argmax(y)))
