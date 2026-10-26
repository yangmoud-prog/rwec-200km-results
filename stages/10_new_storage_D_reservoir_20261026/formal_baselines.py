"""Independent 20 km reach parameters and dt, trained on causal cascade outputs.
No observed future intermediate flow replaces predicted flow in validation or training input.
"""
from common10 import *
from formulas import simulate,hourly
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize
from scipy.stats import qmc
from numba import njit
import time

@njit(cache=True)
def loss(theta,model,dt,inputs,truth,den,v0):
 total=0.
 for ev in range(len(inputs)):
  q,s,v,ok,step,bal=simulate(theta,model,inputs[ev],dt,np.array([v0[ev]]),1)
  if not ok:return 1000.+(len(inputs)-ev)/len(inputs)
  pred=hourly(q,dt,truth.shape[1])[1];total+=np.sum((pred-truth[ev])**2)/den[ev]
 return total/len(inputs)

def prepare(domain,names,cfg):
 geom,bc={'DUF':('uniform','fixed'),'DUR':('uniform','reservoir'),'DCF':('complex','fixed'),'DCR':('complex','reservoir')}[domain];ds=[]
 for e in names:
  f=D/'models'/f"{e}_{geom}_{bc}_{cfg['version']}";quality=read(f/'quality.json')
  if not quality['numerical_gate']:raise RuntimeError('HEC quality gate not passed: '+f.name)
  raw=load(f/'extracted.npz');keep=raw['time']>=0;d={k:(v[keep] if v.ndim and len(v)==len(keep) else v) for k,v in raw.items()};cum=cumulative_trapezoid(d['area'][0],d['x']*1000,initial=0);initial=[np.interp(200-j*20,d['x'],cum)-np.interp(180-j*20,d['x'],cum) for j in range(10)];d['initial_volumes']=np.array(initial)
  expected=np.arange(cfg['end_h']+1,dtype=float)
  if not np.array_equal(d['time'],expected):raise ValueError('Forecast target must include t=0 and every complete hourly interval: '+f.name)
  ds.append(d)
 return ds

def fit(domain):
 cfg=read(D/'config/D_frozen.json');p=read(P/'provenance/frozen_protocol.json');ds=prepare(domain,cfg['train'],cfg);hours=cfg['end_h'];fine=np.arange(0,hours+.125,.25);fits={};log=[]
 protocol=dict(domain=domain,train=cfg['train'],validation=cfg['validation'],method='greedy upstream-to-downstream calibration to training targets, using own causal upstream predictions',parameters_independent_per_20km=True,dt_independent_per_20km=True,dt_candidates_h=p['dt_candidates_h'],objective='training event mean NSE at reach outlet',calls_per_model_per_reach_per_dt=2000,total_budget_per_model=80000,bounds=p['bounds'],validation_read_by_fitter=False,sealed_tests_read=0)
 save(P/f'PROTOCOL_{domain}.json',protocol)
 for model,name in enumerate(['original','new','general','NLM']):
  upstream=np.array([np.interp(fine,d['time'],d['inflow']) for d in ds]);fitlist=[]
  for reach in range(10):
   target=180-20*reach;ys=np.array([d['q'][:,np.argmin(abs(d['x']-target))] for d in ds]);truth=.5*(ys[:,1:]+ys[:,:-1]);den=np.sum((truth-truth.mean(axis=1,keepdims=True))**2,axis=1);v0=np.array([d['initial_volumes'][reach] for d in ds]);best=None
   for dt in p['dt_candidates_h']:
    dst=P/'provenance'/f'fit_{domain}_{name}_reach{reach:02d}_dt{dt:g}.json';grid=np.arange(0,hours+dt/2,dt);inputs=np.array([np.interp(grid,fine,row) for row in upstream]);lo,hi=np.array(p['bounds'][name]).T
    if dst.exists():result=read(dst)
    else:
     history=[];calls=0;fails=0;start=time.perf_counter()
     def fun(z):
      nonlocal calls,fails
      if calls>=2000:return 1e6
      z=np.clip(z,0,1);value=float(loss(lo+z*(hi-lo),model,dt,inputs,truth,den,v0));calls+=1;fails+=int(value>=1000);history.append((value,z.copy()));return value
     starts=list(qmc.LatinHypercube(len(lo),seed=202610260+model+reach*10).random(120));seed={'original':[1,.5,1],'new':[.5,.5,1],'general':[0,1,.5,1],'NLM':[1.3,.05,.7]}[name];starts[0]=(np.array(seed)-lo)/(hi-lo)
     if reach>0:starts.append((np.array(fitlist[-1]['theta'])-lo)/(hi-lo))
     if name=='general':
      # Nested seeds consume ordinary evaluation budget; no extra optimization calls.
      ot=read(P/'provenance'/f'fit_{domain}_original_reach{reach:02d}_dt{dt:g}.json')['theta'];nt=read(P/'provenance'/f'fit_{domain}_new_reach{reach:02d}_dt{dt:g}.json')['theta'];starts += [(np.array([0,ot[0],ot[1],ot[2]])-lo)/(hi-lo),(np.array([nt[0],0,nt[1],nt[2]])-lo)/(hi-lo)]
     for z in starts:fun(z)
     for _,z in sorted(history,key=lambda v:v[0])[:4]:minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=min(400,2000-calls),xatol=1e-7,fatol=1e-10))
     while calls<2000:
      before=calls;z=min(history,key=lambda v:v[0])[1];minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=2000-calls,xatol=1e-8,fatol=1e-11))
      if calls==before:break
     value,z=min(history,key=lambda v:v[0]);result=dict(domain=domain,model=name,model_id=model,reach=reach,outlet_x_km=target,dt_h=dt,theta=(lo+z*(hi-lo)).tolist(),loss=value,training_NSE=1-value,evaluations=calls,failed_evaluations=fails,seconds=time.perf_counter()-start,boundary_parameters=[i for i,v in enumerate(z) if min(v,1-v)<.001]);save(dst,result)
    log.append(result)
    if best is None or result['loss']<best['loss']:best=result
   fitlist.append(best);dt=best['dt_h'];grid=np.arange(0,hours+dt/2,dt);newup=[]
   for ev in range(len(ds)):
    inp=np.interp(grid,fine,upstream[ev]);q,s,v,ok,step,bal=simulate(np.array(best['theta']),model,inp,dt,np.array([v0[ev]]),1);assert ok;newup.append(np.interp(fine,grid,q[1]))
   upstream=np.array(newup);save(P/f'PARAMETERS_{domain}_{name}.json',dict(model=name,domain=domain,reaches=fitlist,complete=len(fitlist)==10));print(domain,name,reach,best['training_NSE'],best['dt_h'],flush=True)
  fits[name]=fitlist
 frame(f'{domain}_calibration.csv',log);save(P/'provenance'/f'{domain}_calibration_frozen.json',dict(complete=True,validation_read=False,train=cfg['train'],budget_per_model=80000))

def predict(d,fit,hours,factor=1.):
 fine=np.arange(0,hours+.125,.25);up=np.interp(fine,d['time'],d['inflow']);qs=[up];native=[];quality=[]
 for reach,f in enumerate(fit):
  dt=f['dt_h']*factor;grid=np.arange(0,hours+dt/2,dt);inp=np.interp(grid,fine,up);q,s,v,ok,step,bal=simulate(np.array(f['theta']),f['model_id'],inp,dt,np.array([d['initial_volumes'][reach]]),1)
  quality.append(dict(reach=reach,ok=bool(ok),failstep=int(step),balance_m3s=float(bal/(dt*3600)),min_Q=float(q[1].min())))
  if not ok:up=np.full_like(fine,np.nan)
  else:up=np.interp(fine,grid,q[1])
  qs.append(up);native.append(dict(dt_h=dt))
 qq=np.array(qs);hour=hourly(qq,.25,hours);return dict(time=fine,q=qq,hourly=hour),quality

def evaluate(domain):
 assert (P/'provenance'/f'{domain}_calibration_frozen.json').exists();cfg=read(D/'config/D_frozen.json');rows=[];qual=[];sens=[];ds=prepare(domain,cfg['train']+cfg['validation'],cfg)
 fits={name:read(P/f'PARAMETERS_{domain}_{name}.json')['reaches'] for name in ['original','new','general','NLM']}
 for e,d in zip(cfg['train']+cfg['validation'],ds):
  subset='train' if e in cfg['train'] else 'validation';truth=np.array([d['q'][:,np.argmin(abs(d['x']-x))] for x in TARGETS]);truth=.5*(truth[:,1:]+truth[:,:-1])
  for name,fitlist in fits.items():
   pred,quality=predict(d,fitlist,cfg['end_h']);np.savez_compressed(P/'predictions'/f'{domain}_{e}_{name}.npz',**pred);complete=all(v['ok'] for v in quality)
   for v in quality:qual.append(dict(domain=domain,event=e,subset=subset,model=name,**v))
   for j,x in enumerate(TARGETS):rows.append(dict(domain=domain,event=e,subset=subset,model=name,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),complete=complete,**metrics(truth[j],pred['hourly'][j+4])))
   if subset=='validation':
    for fac in [.5,2.]:
     pred2,qu2=predict(d,fitlist,cfg['end_h'],fac)
     for j,x in enumerate(TARGETS):sens.append(dict(domain=domain,event=e,model=name,dt_factor=fac,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),complete=all(v['ok'] for v in qu2),**metrics(truth[j],pred2['hourly'][j+4])))
 frame(f'{domain}_event_metrics.csv',rows);frame(f'{domain}_solver_quality.csv',qual);frame(f'{domain}_dt_sensitivity.csv',sens);print(pd.DataFrame(rows).query("subset=='validation'").groupby('model').NSE.mean().to_string())
 save(P/f'{domain}_completed.json',dict(completed=True,validation_after_frozen=True,sealed_tests=0,train=cfg['train'],validation=cfg['validation']))

if __name__=='__main__':
 domain=sys.argv[1]
 if len(sys.argv)>2 and sys.argv[2]=='evaluate':evaluate(domain)
 elif len(sys.argv)>2 and sys.argv[2]=='fit':fit(domain)
 else:fit(domain);evaluate(domain)
