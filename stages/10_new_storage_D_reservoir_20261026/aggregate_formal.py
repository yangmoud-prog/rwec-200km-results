from common10 import *
from formal_baselines import prepare,predict

def main():
 domains=['DUF','DUR','DCF','DCR'];cfg=read(D/'config/D_frozen.json');df=pd.concat([pd.read_csv(P/'tables'/f'{d}_event_metrics.csv') for d in domains],ignore_index=True)
 frame('D_event_metrics.csv',df);cols=['NSE','RMSE','KGE','peak_pct','peak_time_h'];summary=df.groupby(['domain','subset','model','x_from_outlet_km','distance_from_inlet_km'])[cols].agg(['mean','median','min','max']);summary.columns=['_'.join(c) for c in summary.columns];frame('D_section_statistics.csv',summary.reset_index())
 paired=df.pivot(index=['domain','subset','event','x_from_outlet_km','distance_from_inlet_km'],columns='model',values='NSE').reset_index()
 for name in ['new','general','NLM']:paired[name+'_minus_original']=paired[name]-paired['original']
 frame('D_paired_improvements.csv',paired)
 factorial=df.pivot(index=['subset','event','model','x_from_outlet_km','distance_from_inlet_km'],columns='domain',values='RMSE').reset_index();factorial['geometry_at_fixed']=factorial.DCF-factorial.DUF;factorial['backwater_uniform']=factorial.DUR-factorial.DUF;factorial['backwater_complex']=factorial.DCR-factorial.DCF;factorial['interaction']=factorial.DCR-factorial.DCF-factorial.DUR+factorial.DUF;frame('D_factorial_optimized_error.csv',factorial)
 quality=pd.concat([pd.read_csv(P/'tables'/f'{d}_solver_quality.csv') for d in domains]);frame('D_solver_quality.csv',quality);sens=pd.concat([pd.read_csv(P/'tables'/f'{d}_dt_sensitivity.csv') for d in domains]);frame('D_fixed_parameter_dt_sensitivity.csv',sens)
 # Storage/release phase indicators are reference-output diagnostics only; never selection or prediction features.
 phase=[];dyn=[]
 for domain,geom in [('DUR','uniform'),('DCR','complex')]:
  for e in cfg['events']:
   eid=e['event_id'];folder=D/'models'/f"{eid}_{geom}_reservoir_{cfg['version']}";d=load(folder/'extracted.npz');pair=load(folder/'paired.npz');keep=d['time']>=0;time=d['time'][keep];b=pair['b_true'][keep];dv=pair['extra_volume'][keep];rate=np.diff(dv)/3600;bm=.5*(b[:-1]+b[1:]);hecmean=.5*(d['q'][keep][:-1]+d['q'][keep][1:]);trueb=np.array([np.interp(bb,d['x'],q) for bb,q in zip(bm,hecmean)])
   for name in ['original','new','general','NLM']:
    pred=load(P/'predictions'/f'{domain}_{eid}_{name}.npz');hm=pred['hourly'];pb=np.array([np.interp(bb,np.arange(0,201,20),hm[::-1,i]) for i,bb in enumerate(bm)]);dyn.append(dict(domain=domain,event=eid,subset=e['development_split'],model=name,position='same HEC true endpoint; evaluation location only',**metrics(trueb,pb)))
    for j,x in enumerate(TARGETS):
     iy=int(np.argmin(abs(d['x']-x)));trueq=hecmean[:,iy];pq=hm[j+4]
     for label,mask in [('storage',rate>10),('release',rate<-10),('weak_or_recovered',abs(rate)<=10)]:
      if mask.sum()<3:continue
      phase.append(dict(domain=domain,event=eid,subset=e['development_split'],model=name,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),phase=label,hours=int(mask.sum()),RMSE_m3s=float(np.sqrt(np.mean((pq[mask]-trueq[mask])**2))),mean_error_m3s=float(np.mean(pq[mask]-trueq[mask]))))
 frame('D_storage_release_phase_errors.csv',phase);frame('D_same_endpoint_baseline_metrics.csv',dyn)
 # Prediction-input audit: future HEC fields can be replaced without affecting the baseline.
 d=prepare('DCR',[cfg['validation'][0]],cfg)[0];fit=read(P/'PARAMETERS_DCR_original.json')['reaches'];a,_=predict(d,fit,cfg['end_h']);changed={k:v.copy() if isinstance(v,np.ndarray) else v for k,v in d.items()}
 for k in ['q','z','area','reservoir_H','reservoir_V','release','velocity','total_storage']:
  if k in changed:changed[k][:]=np.nan
 z,_=predict(changed,fit,cfg['end_h']);assert np.array_equal(a['q'],z['q'])
 changed={k:v.copy() if isinstance(v,np.ndarray) else v for k,v in d.items()};cut=360;changed['inflow'][changed['time']>cut]*=1.5;z,_=predict(changed,fit,cfg['end_h']);assert np.max(abs(a['q'][:,a['time']<=cut]-z['q'][:,z['time']<=cut]))==0
 save(P/'provenance/causality_check.json',dict(future_HEC_state_scramble_identical=True,upstream_future_suffix_perturbation_preserves_earlier_predictions=True,forecast_inputs=['upstream forecast','known initial geometric volumes','frozen parameters'],coupled_candidates_run=False,reason='oracle geometric and flux closure gates failed',sealed_tests_read=0))
 print(df.query("subset=='validation'").groupby(['domain','model'])[['NSE','RMSE','KGE']].mean().to_string())
if __name__=='__main__':main()
