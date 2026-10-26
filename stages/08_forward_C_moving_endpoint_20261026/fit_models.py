from common import *
from models import loss
from scipy.stats import qmc
from scipy.optimize import minimize
import time

NAMES={0:'RWEC0',1:'LM',2:'NLM',3:'MC_oracle',4:'MC_autonomous'}
BOUNDS={0:[(-3,6),(-1.6,2.3),(0,6)],1:[(-3,3),(.0,.5)],2:[(-3,3),(0,.5),(.4,2.0)],3:[(-3,6),(-1.6,2.3),(0,6),(0,1.2),(0,24)],4:[(-3,6),(-1.6,2.3),(0,6),(0,1.2),(0,24)]}
def main():
 plan=read(P/'provenance/C_frozen.json');budget=1000;allfits={};records=[]
 packets=[load(P/'data'/f'{e}.npz') for e in plan['train']];truth=np.array([d['truth'] for d in packets]);den=np.sum((truth-truth.mean(axis=2,keepdims=True))**2,axis=2)
 y0=packets[0]['initial_depth'];v0=packets[0]['initial_volumes'];deadline=read(W/'continuation/two_window_budget_20261026.json')['second_window_deadline_unix']
 save(P/'provenance/calibration_contract.json',dict(max_evaluations_per_model_per_dt=budget,dt_h=[.25,.5,1.,2.],objective='mean training NSE over14events and7fixedsections',bounds=BOUNDS,validation_used_in_selection=False,spatial_parameters='same homogeneous20km parameters applied10times; common selected subreach dt',stages='128LHS starts + bounded Nelder Mead restarts; equal total objective call budget',future_stage_only_in_oracle_length_labels=True))
 for model,name in NAMES.items():
  best=None
  for dt in [.25,.5,1.,2.]:
   if time.time()>deadline-2400:raise RuntimeError('Reserve40min for deliverables before final window closes')
   out=P/'provenance'/f'fit_{name}_dt{dt:g}.json'
   if out.exists():
    fit=read(out);records.append(fit)
    if best is None or fit['loss']<best['loss']:best=fit
    continue
   grid=np.arange(0,552+dt/2,dt);inputs=np.array([hec.boundary(plan['events'][e],grid,'ref')[0] for e in plan['train']])
   lengths=np.array([np.interp(grid,GRID,d['b_true' if model==3 else 'b_pred']) for d in packets])
   bounds=np.array(BOUNDS[model],float);lo,hi=bounds.T;calls=0;failures=0;history=[];start=time.monotonic()
   def objective(z):
    nonlocal calls,failures
    if calls>=budget:return 1.e6
    theta=lo+np.clip(z,0,1)*(hi-lo);val=float(loss(theta,model,dt,inputs,lengths,truth,den,plan['selected_slope'],y0,v0));calls+=1
    if val>=1000:failures+=1
    history.append((val,z.copy()))
    return val
   starts=qmc.LatinHypercube(len(lo),seed=20261026+model).random(128)
   if model==0:starts[0]=(np.array([1.,.5,1.])-lo)/(hi-lo)
   if model in (1,2):starts[0]=(np.array([1.,.05]+([1.] if model==2 else []))-lo)/(hi-lo)
   if model in (3,4):
    base=read(P/'provenance'/f'fit_RWEC0_dt{dt:g}.json')['theta'];starts[0]=(np.r_[base,0.,0.]-lo)/(hi-lo)
   for z in starts:objective(z)
   initial=sorted(history,key=lambda a:a[0])[:4]
   for _,z in initial:
    remaining=budget-calls
    if remaining<=0:break
    allowance=max(1,remaining//(len(initial)-initial.index((_,z)))) if False else min(240,remaining)
    minimize(objective,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=allowance,xatol=1.e-5,fatol=1.e-7))
   while calls<budget:
    seed=min(history,key=lambda a:a[0])[1]
    before=calls;minimize(objective,seed,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=budget-calls,xatol=1.e-6,fatol=1.e-9))
    if calls==before:break
   value,z=min(history,key=lambda a:a[0]);theta=lo+np.clip(z,0,1)*(hi-lo)
   fit=dict(model=name,model_id=model,dt_h=dt,theta=theta.tolist(),loss=value,mean_training_NSE=1-value,evaluations=calls,failed_parameter_evaluations=failures,elapsed_s=time.monotonic()-start,created=hec.now())
   save(out,fit);records.append(fit);print(json.dumps(fit),flush=True)
   if best is None or fit['loss']<best['loss']:best=fit
  allfits[name]=best;save(P/'provenance/frozen_parameters.json',allfits)
  frame('calibration_search.csv',records)
 print('ALL_FIVE_MODELS_CALIBRATED')

if __name__=='__main__':main()
