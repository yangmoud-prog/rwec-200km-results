from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
OUT=Path(__file__).resolve().parent
A=OUT.parents[1]
WORK=A.parents[1]
B=WORK/'variants/mild_control_2km_60s_pilot_v1'
TARGETS=np.array([20,40,60,80,100,120])
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def frame(name,rows):pd.DataFrame(rows).to_csv(OUT/'tables'/name,index=False,encoding='utf-8-sig')
def loadnp(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dirs():
    for x in ['data','tables','fits','figures','provenance','predictions']:(OUT/x).mkdir(exist_ok=True)
def extent(dz,x,threshold=.1):
    mask=dz>threshold;ids=np.flatnonzero(mask)
    if not len(ids):return 0.,False,0
    j=ids[-1];segments=int(np.sum(np.diff(np.r_[False,mask,False].astype(int))==1))
    if j==len(x)-1:return float(x[j]),True,segments
    f=(dz[j]-threshold)/(dz[j]-dz[j+1])
    return float(x[j]+f*(x[j+1]-x[j])),False,segments
def phase(t,e):
    if t<e['tb_h']:return '上涨前'
    if t<e['Hpeak_h']:return '上涨'
    if t<e['Hfall_start_h']:return '维持'
    if t<e['Hend_h']:return '消落'
    return '边界恢复后'

