from pathlib import Path
import json,sys,hashlib
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
P=Path(__file__).resolve().parent;W=P.parents[3];C=W/'variants/forward_dynamic_C_20261026';A=P.parents[1]
sys.path.insert(0,str(C/'scripts'));import hec200 as hec
GRID=np.arange(553,dtype=float);TARGETS=np.array([120,100,80,60,40,20,0]);NODES=(200-TARGETS)//20
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,obj):hec.dump(Path(p),obj)
def frame(name,rows):pd.DataFrame(rows).to_csv(P/'tables'/name,index=False,encoding='utf-8-sig')
def load(p):
 with np.load(p) as f:return {k:f[k] for k in f.files}
def extent(dh,x,threshold=.1):
 # x ascending from outlet. Farthest upstream positive threshold crossing; no smoothing.
 mask=dh>=threshold;ix=np.flatnonzero(mask)
 if not len(ix):return 0.,0.,0,False
 k=ix[-1];far=x[k]
 if k<len(x)-1:far=x[k]+(x[k+1]-x[k])*(dh[k]-threshold)/(dh[k]-dh[k+1])
 connected=0.
 if mask[0]:
  gap=np.flatnonzero(~mask);j=gap[0]-1 if len(gap) else len(x)-1;connected=x[j]
  if j<len(x)-1:connected=x[j]+(x[j+1]-x[j])*(dh[j]-threshold)/(dh[j]-dh[j+1])
 islands=int(np.sum(mask & np.r_[True,~mask[:-1]]))
 return float(far),float(connected),islands,bool(mask[-1])
def metric(y,p,t=None):
 y=np.asarray(y);p=np.asarray(p);t=np.arange(len(y)) if t is None else np.asarray(t)
 if not np.all(np.isfinite(p)):return dict(NSE=np.nan,RMSE=np.nan,KGE=np.nan,peak_error_pct=np.nan,peak_time_error_h=np.nan)
 den=np.sum((y-y.mean())**2);nse=1-np.sum((p-y)**2)/den if den>0 else np.nan
 corr=np.corrcoef(y,p)[0,1] if np.std(y)>0 and np.std(p)>0 else np.nan
 kge=1-np.sqrt((corr-1)**2+(np.std(p)/np.std(y)-1)**2+(np.mean(p)/np.mean(y)-1)**2) if np.std(y)>0 else np.nan
 return dict(NSE=float(nse),RMSE=float(np.sqrt(np.mean((p-y)**2))),KGE=float(kge),peak_error_pct=float(100*(max(p)-max(y))/max(y)),peak_time_error_h=float(t[np.argmax(p)]-t[np.argmax(y)]))
