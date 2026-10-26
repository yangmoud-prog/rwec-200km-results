"""Uniform extra optimizer budget, decided from training-only nested-model search failure."""
from common10 import *
from formulas import objective
from scipy.optimize import minimize
from scipy.stats import qmc
import time

def main(domain):
 p=read(P/'provenance/frozen_protocol.json');ds=[data(domain,e) for e in p['train']]
 truth=np.array([.5*(d['truth'][:,1:]+d['truth'][:,:-1]) for d in ds]);den=np.sum((truth-truth.mean(axis=2,keepdims=True))**2,axis=2);v0=np.array([d['initial_volumes'] for d in ds]);hours=truth.shape[2];bestall={};records=[]
 save(P/'provenance'/f'training_optimizer_amendment_{domain}.json',dict(reason='general family nested originals but first search did not recover nested training fits at all dt; add exactly 1000 calls to every model/dt before any Stage10 validation evaluation',extra_calls_each=1000,validation_accessed=False))
 for model,name in enumerate(['original','new','general','NLM']):
  best=None
  for dt in p['dt_candidates_h']:
   dst=P/'provenance'/f'refined_{domain}_{name}_dt{dt:g}.json'
   if dst.exists():result=read(dst)
   else:
    first=read(P/'provenance'/f'fit_{domain}_{name}_dt{dt:g}.json');lo,hi=np.array(p['bounds'][name]).T
    grid=np.arange(0,hours+dt/2,dt);inputs=np.array([np.interp(grid,d['time'],d['inflow']) for d in ds]);history=[];calls=0;fails=0;start=time.perf_counter()
    def fun(z):
     nonlocal calls,fails
     if calls>=1000:return 1e6
     z=np.clip(z,0,1);value=float(objective(lo+z*(hi-lo),model,dt,inputs,truth,den,v0));history.append((value,z.copy()));calls+=1;fails+=int(value>=1000);return value
    seeds=[first['theta']]
    if name=='general':
     old=read(P/'provenance'/f'refined_{domain}_original_dt{dt:g}.json')['theta'];new=read(P/'provenance'/f'refined_{domain}_new_dt{dt:g}.json')['theta'];seeds += [[0,old[0],old[1],old[2]],[new[0],0,new[1],new[2]]]
    starts=[(np.array(s)-lo)/(hi-lo) for s in seeds]
    starts+=list(qmc.LatinHypercube(len(lo),seed=p['seed']+100+model).random(96))
    for z in starts:fun(z)
    initial=sorted(history,key=lambda v:v[0])[:3]
    for _,z in initial:minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=min(260,1000-calls),xatol=1e-7,fatol=1e-10))
    while calls<1000:
     before=calls;z=min(history,key=lambda v:v[0])[1]
     minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=1000-calls,xatol=1e-8,fatol=1e-11))
     if before==calls:break
    value,z=min(history,key=lambda v:v[0]);theta=lo+z*(hi-lo)
    result=dict(domain=domain,model=name,model_id=model,dt_h=dt,theta=theta.tolist(),loss=value,mean_training_NSE=1-value,evaluations=first['evaluations']+calls,failed_parameter_evaluations=first['failed_parameter_evaluations']+fails,elapsed_s=first['elapsed_s']+time.perf_counter()-start,boundary_parameters=[i for i,v in enumerate(z) if min(v,1-v)<.001]);save(dst,result);print(json.dumps(result),flush=True)
   records.append(result)
   if best is None or result['loss']<best['loss']:best=result
  bestall[name]=best
 save(P/f'PARAMETERS_{domain}.json',bestall);frame(f'calibration_{domain}.csv',records);save(P/'provenance'/f'parameters_refined_frozen_{domain}.json',dict(complete=True,objective_calls_per_model=8000,validation_read_by_fitter=False))
if __name__=='__main__':main(sys.argv[1])
