from common10 import *
from formulas import objective
from scipy.stats import qmc
from scipy.optimize import minimize
import time,argparse

def main(domain):
 p=read(P/'provenance/frozen_protocol.json');names=p['train'];ds=[data(domain,e) for e in names]
 truth=np.array([.5*(d['truth'][:,1:]+d['truth'][:,:-1]) for d in ds]);den=np.sum((truth-truth.mean(axis=2,keepdims=True))**2,axis=2)
 v0=np.array([d['initial_volumes'] for d in ds]);hours=truth.shape[2];records=[];bestall={}
 for model,name in enumerate(['original','new','general','NLM']):
  best=None
  for dt in p['dt_candidates_h']:
   dst=P/'provenance'/f'fit_{domain}_{name}_dt{dt:g}.json'
   if dst.exists():
    result=read(dst)
   else:
    grid=np.arange(0,hours+dt/2,dt);inputs=np.array([np.interp(grid,d['time'],d['inflow']) for d in ds]);lo,hi=np.array(p['bounds'][name]).T
    budget=p['budget_per_model_per_dt'];history=[];calls=0;fails=0;start=time.perf_counter()
    def fun(z):
     nonlocal calls,fails
     if calls>=budget:return 1e6
     z=np.clip(z,0,1);theta=lo+z*(hi-lo);value=float(objective(theta,model,dt,inputs,truth,den,v0));calls+=1;fails+=int(value>=1000);history.append((value,z.copy()));return value
    starts=qmc.LatinHypercube(len(lo),seed=p['seed']+model).random(100)
    seed={'original':[1,.5,1],'new':[.5,.5,1],'general':[0,1,.5,1],'NLM':[1.3,.05,.7]}[name]
    starts[0]=(np.array(seed)-lo)/(hi-lo)
    for z in starts:fun(z)
    initial=sorted(history,key=lambda v:v[0])[:4]
    for _,z in initial:minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=min(200,budget-calls),xatol=1e-6,fatol=1e-9))
    while calls<budget:
     z=min(history,key=lambda v:v[0])[1];before=calls
     minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=budget-calls,xatol=1e-8,fatol=1e-11))
     if calls==before:break
    value,z=min(history,key=lambda v:v[0]);theta=lo+np.clip(z,0,1)*(hi-lo)
    result=dict(domain=domain,model=name,model_id=model,dt_h=dt,theta=theta.tolist(),loss=value,mean_training_NSE=1-value,evaluations=calls,failed_parameter_evaluations=fails,elapsed_s=time.perf_counter()-start,boundary_parameters=[i for i,v in enumerate(z) if min(v,1-v)<.001])
    save(dst,result);print(json.dumps(result),flush=True)
   records.append(result)
   if best is None or result['loss']<best['loss']:best=result
  bestall[name]=best;save(P/f'PARAMETERS_{domain}.json',bestall);frame(f'calibration_{domain}.csv',records)
 save(P/'provenance'/f'parameters_frozen_{domain}.json',dict(complete=True,models=list(bestall),validation_read_by_fitter=False))
if __name__=='__main__':main(sys.argv[1])
