from common10 import *
from fit_A import predict,prepare
from formulas import simulate,equilibrium
import time

def masks(nodes):
 h=len(nodes)-1;tm=np.arange(h)+.5;negative=np.flatnonzero(nodes<0);cut=float(negative[-1]+12) if len(negative) else 0.
 return dict(full=np.ones(h,dtype=bool),nonnegative_hour_endpoints=(nodes[:-1]>=0)&(nodes[1:]>=0),after_last_negative_plus12h=tm>=cut+.5)

def main():
 cfg=read(P/'provenance/frozen_protocol.json');assert read(P/'provenance/all_parameters_frozen.json')['complete'];fits={m:read(P/f'PARAMETERS_A_{m}.json')['reaches'] for m in MODELS};rows=[];quality=[];sens=[];storage=[];balance=[];runtime=[];negative=[]
 for e in cfg['train']+cfg['validation']:
  d=load(P/'data'/f'{e}.npz');subset='train' if e in cfg['train'] else 'validation';keep=d['time']>=0;tm=np.arange(int(d['time'][-1]))+.5;var=d['q'][keep];ref=d['q_ref'][keep];true=.5*(var[:-1]+var[1:]);true_ref=.5*(ref[:-1]+ref[1:]);preds={}
  for j,x in enumerate(TARGETS):
   col=int(np.argmin(abs(d['x']-x)));mm=masks(var[:,col]);negative.append(dict(event=e,subset=subset,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),minimum_Q_m3s=float(var[:,col].min()),negative_nodes=int((var[:,col]<0).sum()),full_hours=len(tm),excluded_hours=int((~mm['nonnegative_hour_endpoints']).sum()),whole_event_nonnegative=bool(np.all(var[:,col]>=0))))
  for m,fit in fits.items():
   started=time.perf_counter();pred,qual=predict(d,fit);elapsed=time.perf_counter()-started;preds[m]=pred;np.savez_compressed(P/'predictions'/f'{e}_{m}.npz',**pred);complete=all(q['ok'] for q in qual);runtime.append(dict(event=e,subset=subset,model=m,seconds=elapsed))
   for q in qual:quality.append(dict(event=e,subset=subset,model=m,**q))
   for j,x in enumerate(TARGETS):
    col=int(np.argmin(abs(d['x']-x)));mm=masks(var[:,col]);pq=pred['hourly'][j+4]
    for scene,truth in [('dynamic',true[:,col]),('fixed5_same_forecast',true_ref[:,col])]:
     for window,mask in mm.items():
      rows.append(dict(event=e,subset=subset,model=m,scenario=scene,window=window,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),complete=complete,scored_hours=int(mask.sum()),**metrics(truth[mask],pq[mask],tm[mask])))
   # Water balance at native solver and across independently selected time grids.
   after=pred['time']>=0;ti=pred['time'][after];qq=pred['q'][:,after];vv=np.cumsum(pred['reach_bookkeeping_volume'],axis=0)[:,after];pp=np.cumsum(pred['resampling_volume'],axis=0)[:,after];from scipy.integrate import cumulative_trapezoid
   net=cumulative_trapezoid(qq[0]-qq[-1],ti*3600,initial=0);err=(vv[-1]-vv[-1,0])-net;fullvol=float(np.trapz(np.maximum(qq[0]-1000,0),ti*3600));balance.append(dict(event=e,subset=subset,model=m,max_unadjusted_cascade_residual_m3=float(abs(err).max()),fraction_of_flood_volume=float(abs(err).max()/max(fullvol,1)),projection_corrected_max_m3=float(abs(err+pp[-1]-pp[-1,0]).max())))
   cum=cumulative_trapezoid(d['area'][keep],d['x']*1000,axis=1,initial=0)
   for j,x in enumerate(TARGETS):
    count=(200-int(x))//20;geom=cum[:,-1]-cum[:,np.argmin(abs(d['x']-x))];modelvol=vv[count-1,::4];vg=geom-geom[0];vm=modelvol-modelvol[0];storage.append(dict(event=e,subset=subset,model=m,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),effective_volume_anomaly_RMSE_m3=float(np.sqrt(np.mean((vm-vg)**2))),effective_volume_rate_RMSE_m3s=float(np.sqrt(np.mean((np.diff(vm)/3600-np.diff(vg)/3600)**2)))))
   if subset=='validation':
    for fac in [.5,2.]:
     pred2,qual2=predict(d,fit,fac)
     for j,x in enumerate(TARGETS):
      col=int(np.argmin(abs(d['x']-x)))
      for window,mask in masks(var[:,col]).items():sens.append(dict(event=e,model=m,dt_factor=fac,window=window,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),complete=all(q['ok'] for q in qual2),**metrics(true[mask,col],pred2['hourly'][j+4,mask],tm[mask])))
  if subset=='validation':
   for m in MODELS:
    pq=preds[m]['hourly']
    for j,x in enumerate(TARGETS):
     col=int(np.argmin(abs(d['x']-x)));delta=true[:,col]-true_ref[:,col];mask=masks(var[:,col])['nonnegative_hour_endpoints'];yp=np.zeros_like(delta)
     # Same inflow-only prediction for both unknown boundary scenarios; true difference is an evaluation target.
     storage.append(dict(event=e,subset='paired_response',model=m,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),predicted_pair_difference_is_zero=True,paired_response_RMSE_m3s=float(np.sqrt(np.mean(delta[mask]**2))),paired_response_NSE=metrics(delta[mask],yp[mask],tm[mask])['NSE']))
 df=pd.DataFrame(rows);frame('event_metrics.csv',df)
 for w in df.window.unique():frame('event_metrics_'+w+'.csv',df[df.window==w])
 frame('validation_dynamic_event_metrics.csv',df.query("subset=='validation' and scenario=='dynamic'"));frame('solver_quality.csv',quality);frame('fixed_parameter_dt_sensitivity.csv',sens);frame('negative_flow_windows.csv',negative);frame('storage_diagnostics.csv',storage);frame('cascade_mass_balance.csv',balance);frame('runtime.csv',runtime)
 cols=['NSE','RMSE','KGE','peak_pct','peak_time_h'];sm=df.groupby(['subset','scenario','window','model','x_from_outlet_km','distance_from_inlet_km'])[cols].agg(['mean','median','min','max']);sm.columns=['_'.join(v) for v in sm.columns];frame('section_statistics.csv',sm.reset_index());paired=df.pivot(index=['event','subset','scenario','window','x_from_outlet_km','distance_from_inlet_km'],columns='model',values='NSE').reset_index()
 for m in ['new','general','NLM']:paired[m+'_minus_original']=paired[m]-paired.original
 frame('paired_improvements.csv',paired);eg=paired.groupby(['event','subset','scenario','window'])[[m+'_minus_original' for m in ['new','general','NLM']]].mean().reset_index();frame('event_mean_gains.csv',eg)
 frame('score_completeness.csv',df.groupby(['subset','scenario','window','model']).agg(target_rows=('NSE','size'),finite_NSE=('NSE','count'),complete_predictions=('complete','sum')).reset_index())
 # All calibrated parameter sets must retain their own positive equilibrium.
 checks=[];parameters=[]
 for m,fit in fits.items():
  for f in fit:
   t=f['theta'];dt=f['dt_h'];inp=np.full(int(96/dt)+1,1000.);q,s,v,ok,step,bal=simulate(np.array(t),f['model_id'],inp,dt,np.array([1e8]),1);error=float(abs(q-1000).max());checks.append(dict(model=m,reach=f['reach'],passed=bool(ok and error<1e-5 and s.min()>=0 and q.min()>=0),constant_Q_error_m3s=error,min_S_m3=float(s.min()),balance_m3=float(bal)))
   row=dict(model=m,reach=f['reach'],x_from_outlet_km=f['outlet_x_km'],distance_from_inlet_km=200-f['outlet_x_km'],dt_h=dt,training_NSE=f['training_NSE'],boundary_parameter_count=len(f['boundary_parameters']))
   if m=='NLM':row.update(tref_h=float(np.exp(t[0])),X=t[1],m=t[2])
   else:
    a,b,n,tau=(1.,np.exp(t[0]),np.exp(t[1]),t[2]) if m=='original' else (np.exp(t[0]),1.,np.exp(t[1]),t[2]) if m=='new' else (np.exp(t[0]),np.exp(t[1]),np.exp(t[2]),t[3]);r=equilibrium(a,b,n);rho=(r/(r+1))**(n+1)+a*(r/(a*r+b))**(n+1);row.update(a=float(a),b=float(b),n=float(n),tau_h=tau,steady_storage_m3=float(r*1000*dt*3600),local_memory_factor=float(rho),local_decay_h=float(-dt/np.log(rho)))
   parameters.append(row)
 frame('final_parameter_steady_checks.csv',checks);frame('parameters_interpreted.csv',parameters)
 d=load(P/'data'/f"{cfg['validation'][0]}.npz");causal=[]
 for m,fit in fits.items():
  p,_=predict(d,fit);changed={k:v.copy() for k,v in d.items()}
  for k in ['q','q_ref','z','z_ref','area','area_ref']:changed[k][:]=np.nan
  other,_=predict(changed,fit);assert np.array_equal(p['q'],other['q']);changed={k:v.copy() for k,v in d.items()};changed['inflow'][changed['time']>96]*=1.5;other,_=predict(changed,fit);assert np.array_equal(p['q'][:,p['time']<=96],other['q'][:,other['time']<=96]);causal.append(dict(model=m,future_HEC_fields_ignored=True,upstream_suffix_preserves_earlier_predictions=True))
 save(P/'provenance/causality_check.json',dict(checks=causal,sealed_tests_read=0,downstream_future_H_input=False,negative_mask_is_evaluation_only=True))
 finite_balance=[v['fraction_of_flood_volume'] for v in balance if np.isfinite(v['fraction_of_flood_volume'])]
 save(P/'provenance/evaluation_complete.json',dict(complete=True,events=70,validation=21,train=49,metric_rows=len(df),steady_all_passed=all(c['passed'] for c in checks),negative_forecast_reach_runs=sum(v['min_Q']<0 for v in quality),failed_reach_runs=sum(not v['ok'] for v in quality),sensitivity_failed_target_rows=sum(not v['complete'] for v in sens),cascade_residual_max_fraction=max(finite_balance) if finite_balance else None,parameter_revision_after_validation=False))
 print(df.query("subset=='validation' and scenario=='dynamic'").groupby(['window','model'])[['NSE','RMSE','KGE']].mean().to_string(),flush=True)
if __name__=='__main__':main()
