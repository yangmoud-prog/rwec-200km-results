from common import *
from dynamic_core import route,muskingum
from calibrate import MODELS
import time,sys

NODES=[9,8,7,6,5,4]

def simulate(d,model,fit,length,dt=1.):
 t=np.arange(-48.,240.+.1*dt,dt);I=np.interp(t,GRID,d['inflow']);lb=np.interp(t,GRID,length);p=fit['parameters']
 start=time.perf_counter()
 if model=='Muskingum-linear':
  q=muskingum(I,dt,p['K20_h'],p['X']);res=0.;limits=0;mint=0.
 else:q,s,tr,res,limits,mint=route(I,lb,dt,p['k'],p['n'],p['tau20_h'],p['beta'])
 elapsed=time.perf_counter()-start
 return q[::int(round(1/dt))],dict(compute_s=elapsed,water_balance_max_m3=float(res),volume_limits=int(limits),minimum_transit_storage_m3=float(mint))

def main():
 meta=read(OUT/'data/metadata.json');rows=[];checks=[];sensitivity=[];bounds=[];tmask=GRID>0
 for group in ['A','B']:
  fits={model:read(OUT/'fits'/f'{group}_{model}.json') for model in MODELS}
  for m in [a for a in meta if a['group']==group]:
   d=load(OUT/'data'/m['file']);pred=load(OUT/'predictions'/f"endpoint_{group}_{m['event']}.npz")
   for model,fit in fits.items():
    length=d['true_lb'] if model=='RWEC-L-oracle' else pred['central']
    q,quality=simulate(d,model,fit,length);checks.append(dict(group=group,event=m['event'],subset=m['subset'],model=model,**quality))
    np.savez_compressed(OUT/'predictions'/f"flow_{group}_{m['event']}_{model}.npz",time=GRID,q=q[:,NODES])
    for j,xc in enumerate(TARGETS):
     y=d['truth'][:,j];p=q[:,NODES[j]];mask=tmask&(y>=0)
     for window,mask2 in [('positive_common',mask),('all_hours',tmask)]:
      out=metric(y[mask2],p[mask2]);ids=np.flatnonzero(tmask);py=ids[np.argmax(y[tmask])];pp=ids[np.argmax(p[tmask])]
      rows.append(dict(group=group,event=m['event'],subset=m['subset'],model=model,xc_km=int(xc),window=window,n=int(mask2.sum()),whole_forecast_positive=bool(np.all(y[tmask]>=0)),
       **out,peak_error_pct=float(100*(p[pp]-y[py])/y[py]),peak_time_error_h=float(GRID[pp]-GRID[py])))
    if m['subset']=='validation':
     for dt in [.5,.25]:
      qq,qual=simulate(d,model,fit,length,dt)
      for j,xc in enumerate(TARGETS):
       mask=tmask&(d['truth'][:,j]>=0)
       sensitivity.append(dict(group=group,event=m['event'],model=model,dt_h=dt,xc_km=int(xc),**metric(d['truth'][mask,j],qq[mask,NODES[j]]),RMSE_vs_hourly=float(np.sqrt(np.mean((qq[tmask,NODES[j]]-q[tmask,NODES[j]])**2))),**qual))
   if group=='A':
    # Truly unforced before origin: same I and numerically same past water surface.
    if np.max(abs(d['delta_stage'][GRID<=0]))<=.001:
     for j,xc in enumerate(TARGETS):
      diff=d['truth'][tmask,j]-d['reference'][tmask,j]
      bounds.append(dict(event=m['event'],subset=m['subset'],xc_km=int(xc),past_profile_difference_max_m=float(abs(d['delta_stage'][GRID<=0]).max()),
       paired_RMSE_gap_m3s=float(np.sqrt(np.mean(diff**2))),common_prediction_RMSE_lower_bound_m3s=float(np.sqrt(np.mean(diff**2))/2),max_flow_gap_m3s=float(abs(diff).max())))
 frame('flow_event_metrics.csv',rows);frame('routing_quality.csv',checks);frame('fixed_parameter_timestep.csv',sensitivity);frame('indistinguishable_pair_bounds.csv',bounds)
 df=pd.DataFrame(rows);v=df[(df.subset=='validation')&(df.window=='positive_common')]
 agg=v.groupby(['group','model','xc_km'])[['NSE','RMSE','KGE','peak_error_pct','peak_time_error_h']].mean().reset_index();frame('flow_validation_summary.csv',agg)
 comparison=[]
 for (group,xc),g in v.groupby(['group','xc_km']):
  piv=g.pivot(index='event',columns='model',values='NSE');base=piv['RWEC-0']
  for model in MODELS:
   comparison.append(dict(group=group,xc_km=int(xc),model=model,mean_delta_NSE=float((piv[model]-base).mean()),median_delta_NSE=float((piv[model]-base).median()),improved_events=int(((piv[model]-base)>1.e-6).sum()),events=len(piv)))
 frame('paired_improvement.csv',comparison)
 # Regression against previous original RWEC implementation with identical parameters.
 sys.path.insert(0,str(P5));from routing_core import engine
 d=load(OUT/'data'/meta[0]['file']);p=read(OUT/'fits/A_RWEC-0.json')['parameters']
 q=route(d['inflow'],np.zeros(len(GRID)),1.,p['k'],p['n'],p['tau20_h'],0.)[0]
 old=engine(d['inflow'],1.,p['k'],p['n'],p['tau20_h'],0.,0,np.zeros((2,101)),np.zeros((2,101)),np.zeros(101),np.zeros(len(GRID)))[0]
 err=float(abs(q-old).max());assert err<1.e-6,err
 # Training neighbors must exclude the event itself and all validation events.
 neighbors=pd.read_csv(OUT/'tables/endpoint_neighbors.csv');leak=[]
 for group in ['A','B']:
  tr={m['event'] for m in meta if m['group']==group and m['subset']=='train'}
  for _,r in neighbors[neighbors.group==group].iterrows():
   if r.neighbor not in tr or r.neighbor==r.event:leak.append(r.to_dict())
 assert not leak
 save(OUT/'provenance/numerical_and_causality_checks.json',dict(beta_zero_regression_max_m3s=err,neighbor_leakage_count=len(leak),prediction_query_keys=['inflow_forecast','past_H_through_origin'],future_H_prediction_input=False,sealed_tests_opened=0,
  note='Oracle true length is explicitly isolated. Geometric true volumes are used only for diagnostics. Moving-volume balance uses measured profiles at hourly resolution.'))
 print(v.groupby(['group','model'])[['NSE','RMSE','KGE']].mean().to_string(),flush=True)
 save(OUT/'checkpoint.json',dict(stage='first_pass_computations_complete',A_train=49,A_validation=21,B_train=14,B_validation=6,sealed_test_results_opened=0))

if __name__=='__main__':main()
