from common10 import *
from formulas import simulate,hourly
from scipy.optimize import minimize
from scipy.stats import qmc
from numba import njit
import time

@njit(cache=True)
def loss(theta,model,dt,inputs,truth,valid,den,v0):
 total=0.;skip=int(round(48/dt));hours=truth.shape[1]+48
 for ev in range(len(inputs)):
  q,s,v,ok,step,bal=simulate(theta,model,inputs[ev],dt,np.array([v0[ev]]),1)
  if not ok:return 1000.+(len(inputs)-ev)/len(inputs)
  pred=hourly(q,dt,hours)[1,48:];err=pred-truth[ev];total+=np.sum(err*err*valid[ev])/den[ev]
 return total/len(inputs)

def prepare(names):return [load(P/'data'/f'{e}.npz') for e in names]

def fit(name):
 cfg=read(P/'provenance/frozen_protocol.json');model=MODELS.index(name);ds=prepare(cfg['train']);maxh=int(max(d['time'][-1] for d in ds));fine=np.arange(-48,maxh+.125,.25);upstream=np.array([np.interp(fine,d['time'],d['inflow']) for d in ds]);fits=[];allrows=[];seedold=read(P.parent/'05_routing_validation_slope075_20261022/fits/A_core_final.json')['core']
 for reach in range(10):
  target=180-20*reach;truth=np.zeros((len(ds),maxh));valid=np.zeros_like(truth)
  for ev,d in enumerate(ds):
   ids=np.flatnonzero(d['time']>=0);ys=d['q'][ids,np.argmin(abs(d['x']-target))];n=len(ys)-1;truth[ev,:n]=.5*(ys[1:]+ys[:-1]);valid[ev,:n]=1
  count=valid.sum(axis=1);mu=(truth*valid).sum(axis=1)/count;den=(((truth-mu[:,None])**2)*valid).sum(axis=1);assert np.all(den>0);v0=np.array([d['initial_volumes'][reach] for d in ds]);best=None
  for dt in cfg['dt_candidates_h']:
   dst=P/'provenance'/f'fit_A_{name}_reach{reach:02d}_dt{dt:g}.json';grid=np.arange(-48,maxh+dt/2,dt);inputs=np.array([np.interp(grid,fine,row) for row in upstream]);lo,hi=np.array(cfg['bounds'][name]).T
   if dst.exists():result=read(dst)
   else:
    history=[];calls=0;fails=0;started=time.perf_counter()
    def fun(z):
     nonlocal calls,fails
     if calls>=2000:return 1e6
     z=np.clip(z,0,1);value=float(loss(lo+z*(hi-lo),model,dt,inputs,truth,valid,den,v0));calls+=1;fails+=int(value>=1000);history.append((value,z.copy()));return value
    starts=list(qmc.LatinHypercube(len(lo),seed=cfg['seed']+model+reach*10).random(120));seed={'original':[np.log(seedold[0]),np.log(seedold[1]),seedold[2]],'new':[.5,.5,1.],'general':[0.,-3.,.5,1.],'NLM':[1.3,.05,.7]}[name];starts[0]=(np.array(seed)-lo)/(hi-lo)
    if reach>0:starts.append((np.array(fits[-1]['theta'])-lo)/(hi-lo))
    if name=='general':
     # Same call budget; nested starts depend only on fresh original/new training fits.
     st=time.monotonic()
     while not all((P/'provenance'/f'fit_A_{m}_reach{reach:02d}_dt{dt:g}.json').exists() for m in ['original','new']):
      if time.monotonic()-st>7200:raise TimeoutError('Nested training fits not available')
      time.sleep(2)
     ot=read(P/'provenance'/f'fit_A_original_reach{reach:02d}_dt{dt:g}.json')['theta'];nt=read(P/'provenance'/f'fit_A_new_reach{reach:02d}_dt{dt:g}.json')['theta'];starts.extend([(np.array([0.,ot[0],ot[1],ot[2]])-lo)/(hi-lo),(np.array([nt[0],0.,nt[1],nt[2]])-lo)/(hi-lo)])
    for z in starts:fun(z)
    for _,z in sorted(history,key=lambda z:z[0])[:4]:minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=min(400,2000-calls),xatol=1e-7,fatol=1e-10))
    while calls<2000:
     before=calls;z=min(history,key=lambda z:z[0])[1];minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=2000-calls,xatol=1e-8,fatol=1e-11))
     if before==calls:raise RuntimeError('Optimizer made no objective calls')
    value,z=min(history,key=lambda z:z[0]);result=dict(model=name,model_id=model,reach=reach,outlet_x_km=target,dt_h=dt,theta=(lo+z*(hi-lo)).tolist(),loss=value,training_NSE=1-value,evaluations=calls,failed_evaluations=fails,seconds=time.perf_counter()-started,boundary_parameters=[i for i,v in enumerate(z) if min(v,1-v)<.001]);save(dst,result)
   allrows.append(result)
   if best is None or result['loss']<best['loss']:best=result
  fits.append(best);dt=best['dt_h'];grid=np.arange(-48,maxh+dt/2,dt);nextup=[]
  for ev in range(len(ds)):
   inp=np.interp(grid,fine,upstream[ev]);q,s,v,ok,step,bal=simulate(np.array(best['theta']),model,inp,dt,np.array([v0[ev]]),1);assert ok;nextup.append(np.interp(fine,grid,q[1]))
  upstream=np.array(nextup);save(P/f'PARAMETERS_A_{name}.json',dict(model=name,reaches=fits,complete=len(fits)==10));print(name,reach,best['training_NSE'],best['dt_h'],flush=True)
 frame(f'calibration_{name}.csv',allrows);save(P/'provenance'/f'{name}_frozen.json',dict(complete=True,validation_read=False,train=cfg['train'],budget=80000))

def predict(d,fit,factor=1.):
 hours=int(d['time'][-1]);fine=np.arange(-48,hours+.125,.25);up=np.interp(fine,d['time'],d['inflow']);qq=[up];vols=[];states=[];interfaces=[];quality=[]
 from scipy.integrate import cumulative_trapezoid
 for reach,f in enumerate(fit):
  dt=f['dt_h']*factor;grid=np.arange(-48,hours+dt/2,dt);inp=np.interp(grid,fine,up);q,s,v,ok,step,bal=simulate(np.array(f['theta']),f['model_id'],inp,dt,np.array([d['initial_volumes'][reach]]),1);quality.append(dict(reach=reach,ok=bool(ok),failed_step=int(step),min_Q=float(np.nanmin(q[1])),balance_m3s=float(bal/(dt*3600))))
  nativein=np.interp(fine,grid,inp);out=np.interp(fine,grid,q[1]) if ok else np.full_like(fine,np.nan);vols.append(d['initial_volumes'][reach]+cumulative_trapezoid(nativein-out,fine*3600,initial=0));states.append(np.interp(fine,grid,s[0]));interfaces.append(cumulative_trapezoid(up-nativein,fine*3600,initial=0));up=out;qq.append(up)
 qq=np.array(qq);return dict(time=fine,q=qq,hourly=hourly(qq,.25,hours+48)[:,48:],reach_bookkeeping_volume=np.array(vols),effective_state=np.array(states),resampling_volume=np.array(interfaces)),quality

if __name__=='__main__':fit(sys.argv[1])
