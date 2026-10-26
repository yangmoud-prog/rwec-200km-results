from common import *
from models import rwec,muskingum
from joint import simulate
from scipy.stats import qmc
from scipy.optimize import minimize
from numba import njit
import time

NAMES={0:'RWEC0',1:'LM',2:'NLM',4:'Serial8',5:'Joint_oracle',6:'Joint'}
BOUNDS={0:[(-3,6),(-1.6,2.3),(0,6)],1:[(-3,3),(0,.5)],2:[(-3,3),(0,.5),(.4,2.)],4:[(-3,6),(-1.6,2.3),(0,6),(0,1.2),(0,24)],5:[(-3,6),(-1.6,2.5),(-1.38629436112,1.38629436112)],6:[(-3,6),(-1.6,2.5),(-1.38629436112,1.38629436112)]}

@njit(cache=True)
def hourly(q,dt,joint):
 n=int(round(1/dt));out=np.zeros((q.shape[0],288))
 for h in range(288):
  for i in range(n):
   t=h*n+i
   if joint:out[:,h]+=q[:,t+1]/n
   else:out[:,h]+=(q[:,t]+q[:,t+1])/(2*n)
 return out

@njit(cache=True)
def objective_num(theta,model,dt,inputs,bs,truth,den,yinit,y0,v0):
 total=0.
 for ev in range(len(inputs)):
  if model==0 or model==4:
   alpha=theta[3] if model==4 else 0.;relax=theta[4] if model==4 else 0.
   q,ss,tt,bb,v,ok,fail,bal,local=rwec(inputs[ev],bs[ev],dt,np.exp(theta[0]),np.exp(theta[1]),theta[2],alpha,relax,.000075,y0[ev],v0[ev])
  elif model in (1,2):
   q,v,ok,fail,bal=muskingum(inputs[ev],dt,np.exp(theta[0]),theta[1],theta[2] if model==2 else 1.)
  else:
   q,b,h,qb,sb,su,w,it,er,fail,bal=simulate(inputs[ev],dt,np.exp(theta[0]),np.exp(theta[1]),np.exp(theta[2]),yinit[ev],bs[ev],1 if model==5 else 0);ok=fail==0
  if not ok:return 1000.+(len(inputs)-ev)/len(inputs)
  pred=hourly(q,dt,model>=5)[4:11]
  total+=np.sum(np.sum((pred-truth[ev])**2,axis=1)/den[ev])
 return total/(len(inputs)*7)

def main():
 assert (P/'provenance/closure_completed.json').exists()
 protocol=read(P/'provenance/frozen_protocol.json');names=protocol['train'];data=[load(P/'data'/f'B_{e}.npz') for e in names]
 truth=np.array([.5*(d['truth'][:,1:]+d['truth'][:,:-1]) for d in data]);den=np.sum((truth-truth.mean(axis=2,keepdims=True))**2,axis=2)
 yinit=np.array([d['yinit'] for d in data]);y0=np.array([d['initial_depth'] for d in data]);v0=np.array([d['initial_volumes'] for d in data]);records=[];bestall={}
 protocol.update(observation_operator='common hourly interval-average discharge: HEC and traditional outputs trapezoidal, joint native conservative means; endpoint positions evaluated at hourly endpoints',bounds=BOUNDS)
 save(P/'provenance/frozen_protocol.json',protocol)
 for model,name in NAMES.items():
  best=None
  for dt in protocol['dt_candidates_h']:
   dst=P/'provenance'/f'fit_{name}_dt{dt:g}.json'
   if dst.exists():
    result=read(dst);records.append(result)
    if best is None or result['loss']<best['loss']:best=result
    continue
   grid=np.arange(0,288+dt/2,dt);inputs=np.array([np.interp(grid,GRID,d['inflow']) for d in data]);bs=np.array([np.interp(grid,GRID,d['b_true' if model==5 else 'b_pred']) if model in (4,5) else np.zeros(len(grid)) for d in data])
   lo,hi=np.array(BOUNDS[model]).T;budget=protocol['max_objective_calls_per_model_per_dt'];history=[];calls=0;fails=0;start=time.perf_counter()
   def fun(z):
    nonlocal calls,fails
    if calls>=budget:return 1e6
    theta=lo+np.clip(z,0,1)*(hi-lo);value=float(objective_num(theta,model,dt,inputs,bs,truth,den,yinit,y0,v0));calls+=1;fails+=int(value>=1000);history.append((value,z.copy()));return value
   starts=qmc.LatinHypercube(len(lo),seed=20261026+model).random(48)
   seed={0:[1.,.5,1.],1:[1.,.05],2:[1.3,.05,.7],4:[1.,.5,1.,.15,6.],5:[np.log(3),np.log(2),0.],6:[np.log(3),np.log(2),0.]}[model]
   if model==4:
    base=read(P/'provenance'/f'fit_RWEC0_dt{dt:g}.json')['theta'];seed=base+[0.,6.]
   starts[0]=(np.array(seed)-lo)/(hi-lo)
   for z in starts:fun(z)
   initial=sorted(history,key=lambda a:a[0])[:3]
   for _,z in initial:
    minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=min(85,budget-calls),xatol=1e-5,fatol=1e-8))
   while calls<budget:
    z=min(history,key=lambda a:a[0])[1];before=calls
    minimize(fun,z,method='Nelder-Mead',bounds=[(0,1)]*len(lo),options=dict(maxfev=budget-calls,xatol=1e-7,fatol=1e-10))
    if calls==before:break
   value,z=min(history,key=lambda a:a[0]);theta=lo+np.clip(z,0,1)*(hi-lo)
   result=dict(model=name,model_id=model,dt_h=dt,theta=theta.tolist(),loss=value,mean_training_NSE=1-value,evaluations=calls,failed_parameter_evaluations=fails,elapsed_s=time.perf_counter()-start,boundary_parameters=[i for i,v in enumerate(z) if min(v,1-v)<.001])
   save(dst,result);records.append(result);print(json.dumps(result),flush=True)
   if best is None or result['loss']<best['loss']:best=result
  bestall[name]=best;save(P/'PARAMETERS.json',bestall);frame('calibration_search.csv',records)
 print('FROZEN_PARAMETERS_READY',flush=True)
if __name__=='__main__':main()
