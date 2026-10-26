from pathlib import Path
import sys,json,hashlib,os
import numpy as np
import pandas as pd
if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8')
P=Path(__file__).resolve().parent
A=P.parents[1];W=P.parents[3];P10=P.parent/'10_new_storage_D_reservoir_20261026'
TARGETS=np.array([120,100,80,60,40,20,0])
MODELS=['original','new','general','NLM']
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+'.tmp-'+str(os.getpid()));tmp.write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');tmp.replace(p)
def load(p):
 with np.load(p,allow_pickle=False) as d:return {k:d[k] for k in d.files}
def frame(name,rows):pd.DataFrame(rows).to_csv(P/'tables'/name,index=False,encoding='utf-8-sig')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def metrics(y,p,t):
 if len(y)<3 or not np.all(np.isfinite(p)):return dict(NSE=np.nan,RMSE=np.nan,KGE=np.nan,peak_pct=np.nan,peak_time_h=np.nan)
 den=np.sum((y-y.mean())**2);sd=y.std();mu=y.mean();r=np.corrcoef(y,p)[0,1] if p.std()>1e-12 and sd>1e-12 else 0.
 return dict(NSE=float(1-np.sum((p-y)**2)/den) if den>1e-12 else np.nan,RMSE=float(np.sqrt(np.mean((p-y)**2))),KGE=float(1-np.sqrt((r-1)**2+(p.std()/sd-1)**2+(p.mean()/mu-1)**2)) if sd>1e-12 and abs(mu)>1e-12 else np.nan,peak_pct=float(100*(p.max()-y.max())/y.max()),peak_time_h=float(t[np.argmax(p)]-t[np.argmax(y)]))
