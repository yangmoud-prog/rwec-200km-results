from common import *
M=pd.read_csv(OUT/'tables/mechanism_events.csv');P=pd.read_csv(OUT/'tables/mechanism_process.csv')
averages=P.groupby(['group','event','xc_km']).leff_km.mean().rename('mean_leff_km').reset_index();M=M.merge(averages,on=['group','event','xc_km'])
M['peak_lag_h']=M.peak_time_h-M.inflow_peak_time_h
params=[];scores=[]
for group in ['A','B']:
 d=M[M.group==group];tr=d[d.subset=='train']
 for model in ['reach_only','reach_plus_backwater']:
  def features(sub):
   a=np.column_stack([np.ones(len(sub)),200-sub.xc_km.to_numpy()])
   return a if model=='reach_only' else np.column_stack([a,sub.mean_leff_km])
  coef=np.linalg.lstsq(features(tr),tr.peak_lag_h.to_numpy(),rcond=None)[0]
  params.append(dict(group=group,model=model,intercept_h=coef[0],reach_h_per_km=coef[1],backwater_h_per_km=coef[2] if len(coef)>2 else 0.))
  for subset in ['train','validation']:
   sub=d[d.subset==subset];scores.append(dict(group=group,model=model,subset=subset,**metric(sub.peak_lag_h.to_numpy(),features(sub)@coef)))
frame('peak_lag_relation_parameters.csv',params);frame('peak_lag_relation_metrics.csv',scores)
print(pd.DataFrame(scores).to_string(index=False))
