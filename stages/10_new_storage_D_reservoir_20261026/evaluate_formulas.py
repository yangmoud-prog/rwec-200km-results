from common10 import *
from formulas import simulate,hourly,equilibrium
import time

def run(d,fit,dt=None):
 dt=fit['dt_h'] if dt is None else dt;hours=int(round(d['time'][-1]));grid=np.arange(0,hours+dt/2,dt);inp=np.interp(grid,d['time'],d['inflow']);start=time.perf_counter()
 q,s,v,ok,step,bal=simulate(np.array(fit['theta']),fit['model_id'],inp,dt,d['initial_volumes']);elapsed=time.perf_counter()-start
 out=hourly(q,dt,hours)
 if not ok:out[:]=np.nan
 return dict(time=grid,q=q,hourly=out,storage=v,effective_storage=s),dict(complete=bool(ok),failed_step=int(step),runtime_s=elapsed,balance_flow_max_m3s=float(bal/(dt*3600)),min_Q_m3s=float(q.min()),min_effective_storage_m3=float(s.min()),dt_h=float(dt))

def main(domains):
 p=read(P/'provenance/frozen_protocol.json');rows=[];quality=[];sens=[];pars=[]
 for domain in domains:
  assert (P/'provenance'/f'parameters_refined_frozen_{domain}.json').exists()
  fits=read(P/f'PARAMETERS_{domain}.json')
  for name,f in fits.items():
   th=f['theta'];row=dict(domain=domain,model=name,dt_h=f['dt_h'],theta=th)
   if name!='NLM':
    a,b,n,tau=(1,np.exp(th[0]),np.exp(th[1]),th[2]) if name=='original' else ((np.exp(th[0]),1,np.exp(th[1]),th[2]) if name=='new' else (np.exp(th[0]),np.exp(th[1]),np.exp(th[2]),th[3]))
    r=equilibrium(a,b,n);rho=(r/(r+1))**(n+1)+a*(r/(a*r+b))**(n+1)
    row.update(a=a,b=b,n=n,tau_h=tau,equilibrium_ratio=r,effective_steady_storage_per_20km_m3=r*1000*f['dt_h']*3600,linearized_memory_factor=rho,local_recession_time_h=-f['dt_h']/np.log(rho))
   else:row.update(tref_h=np.exp(th[0]),X=th[1],nonlinear_exponent=th[2])
   pars.append(row)
  for e in p['train']+p['validation']:
   d=data(domain,e);subset='train' if e in p['train'] else 'validation';truth=.5*(d['truth'][:,1:]+d['truth'][:,:-1])
   for name,fit in fits.items():
    pred,qual=run(d,fit);np.savez_compressed(P/'predictions'/f'{domain}_{e}_{name}.npz',**pred);quality.append(dict(domain=domain,event=e,subset=subset,model=name,**qual))
    for j,x in enumerate(TARGETS):rows.append(dict(domain=domain,event=e,subset=subset,model=name,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),complete=qual['complete'],**metrics(truth[j],pred['hourly'][j+4])))
    if subset=='validation':
     for factor in [.5,2]:
      pp,qq=run(d,fit,fit['dt_h']*factor)
      for j,x in enumerate(TARGETS):sens.append(dict(domain=domain,event=e,model=name,dt_factor=factor,dt_h=fit['dt_h']*factor,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),complete=qq['complete'],**metrics(truth[j],pp['hourly'][j+4])))
 df=pd.DataFrame(rows);frame('formula_event_metrics.csv',df);frame('formula_solver_quality.csv',quality);frame('formula_parameters_interpreted.csv',pars);frame('formula_fixed_parameter_dt_sensitivity.csv',sens)
 cols=['NSE','RMSE','KGE','peak_pct','peak_time_h'];summary=df.groupby(['domain','subset','model','x_from_outlet_km','distance_from_inlet_km'])[cols].agg(['mean','median','min','max']);summary.columns=['_'.join(c) for c in summary.columns];frame('formula_section_statistics.csv',summary.reset_index())
 paired=df.pivot(index=['domain','subset','event','x_from_outlet_km','distance_from_inlet_km'],columns='model',values='NSE').reset_index()
 for name in ['new','general','NLM']:paired[name+'_minus_original']=paired[name]-paired['original']
 frame('formula_event_improvements.csv',paired)
 save(P/'provenance/formula_evaluation_completed.json',dict(domains=domains,validation_evaluated_after_parameters_frozen=True,validation_used_for_recalibration=False,sealed_tests_read=0,rows=len(rows)))
 print(df.query("subset=='validation'").groupby(['domain','model'])[['NSE','RMSE','KGE']].mean().to_string())
if __name__=='__main__':main(sys.argv[1:] or ['B','C'])
