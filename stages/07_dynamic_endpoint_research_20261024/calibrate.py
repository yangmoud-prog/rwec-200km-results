from common import *
from dynamic_core import batch_loss
from scipy.optimize import differential_evolution,minimize
import argparse,time

MODELS={'RWEC-0':0,'RWEC-L-predicted':1,'RWEC-L-oracle':2,'Muskingum-linear':3}
SEEDS=[20261024,20261025,20261026]

def dataset(group,kind):
 meta=[m for m in read(OUT/'data/metadata.json') if m['group']==group and m['subset']=='train']
 inputs=[];lengths=[];truth=[];masks=[];denoms=[]
 for m in meta:
  d=load(OUT/'data'/m['file']);p=load(OUT/'predictions'/f"endpoint_{group}_{m['event']}.npz")
  inputs.append(d['inflow']);lengths.append(d['true_lb'] if kind==2 else p['central']);truth.append(d['truth'])
  mask=(d['truth']>=0)&(GRID[:,None]>0);masks.append(mask)
  denoms.append([max(1.,np.sum((d['truth'][mask[:,j],j]-d['truth'][mask[:,j],j].mean())**2)) for j in range(6)])
 return tuple(np.ascontiguousarray(x) for x in [inputs,lengths,truth,masks,denoms])

def fit(group,model,iterations=65):
 kind=MODELS[model];data=dataset(group,kind);bounds=[(np.log(1.e-4),np.log(20.)),(np.log(.15),np.log(150.)),(0.,8.)]
 if kind in [1,2]:bounds.append((0.,4.))
 if kind==3:bounds=[(np.log(.5),np.log(12.)),(0.,1.)]
 t0=time.perf_counter();fun=lambda p:float(batch_loss(np.asarray(p),*data,kind))
 warm=(np.array([np.log(1.),np.log(4.),1.,0.])[:len(bounds)] if kind!=3 else np.array([np.log(2.),.2]))
 fun(warm);runs=[]
 for seed in SEEDS:
  start=time.perf_counter()
  de=differential_evolution(fun,bounds,popsize=7,maxiter=iterations,tol=2.e-6,atol=1.e-8,polish=False,seed=seed,updating='immediate')
  local=minimize(fun,de.x,method='Nelder-Mead',bounds=bounds,options={'maxfev':1000,'xatol':1.e-7,'fatol':1.e-9})
  best=local if local.fun<de.fun else de
  runs.append(dict(seed=seed,theta=best.x.tolist(),loss=float(best.fun),de_success=bool(de.success),de_message=str(de.message),de_nfev=int(de.nfev),local_success=bool(local.success),local_message=str(local.message),local_nfev=int(local.nfev),elapsed_s=time.perf_counter()-start))
  save(OUT/'fits'/f'{group}_{model}.partial.json',dict(group=group,model=model,runs=runs))
  print(group,model,seed,'loss',best.fun,'time_s',runs[-1]['elapsed_s'],flush=True)
 # For dynamic candidates, also refine from the fully calibrated nested baseline.
 nested=None
 if kind in [1,2]:
  base=read(OUT/'fits'/f'{group}_RWEC-0.json');xx=np.r_[base['theta'],0.]
  nr=minimize(fun,xx,method='Nelder-Mead',bounds=bounds,options={'maxfev':1000,'xatol':1.e-7,'fatol':1.e-9})
  nested=dict(theta=nr.x.tolist(),loss=float(nr.fun),success=bool(nr.success),nfev=int(nr.nfev))
  runs.append(dict(seed='nested-baseline',**nested))
 best=min(runs,key=lambda r:r['loss']);v=best['theta']
 parameters=dict(K20_h=float(np.exp(v[0])),X=float(v[1]*.5/np.exp(v[0]))) if kind==3 else dict(k=float(np.exp(v[0])),n=float(np.exp(v[1])),tau20_h=float(v[2]),beta=float(v[3]) if kind else 0.)
 result=dict(group=group,model=model,kind=kind,theta=v,parameters=parameters,loss=best['loss'],runs=runs,bounds=bounds,elapsed_s=time.perf_counter()-t0,training_events=len(data[0]),joint_fit=True,test_events=0,criterion='equal event and section normalized squared error, true-positive common mask',bound_hits=[i for i,(lo,hi) in enumerate(bounds) if min(abs(v[i]-lo),abs(v[i]-hi))<1.e-4])
 save(OUT/'fits'/f'{group}_{model}.json',result)
 return result

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--group',required=True);parser.add_argument('--iterations',type=int,default=65);args=parser.parse_args()
 for model in MODELS:
  path=OUT/'fits'/f'{args.group}_{model}.json'
  if path.exists():print('Reuse completed fit',path.name,flush=True);continue
  fit(args.group,model,args.iterations)
