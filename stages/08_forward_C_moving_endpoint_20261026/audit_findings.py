from common import *
p=P;d=pd.read_csv(p/'tables/fixed_section_event_metrics.csv');s=pd.read_csv(p/'tables/model_dt_sensitivity.csv')
phase=d[(d.subset=='validation')&(d.x_from_outlet_km<=20)].groupby(['model','window'])[['NSE','RMSE']].mean().reset_index()
dt=s[s.x_from_outlet_km<=20].groupby(['model','dt_factor'])[['NSE','RMSE','valid']].mean().reset_index()
frame('near_outlet_phase_means.csv',phase.to_dict('records'));frame('near_outlet_dt_means.csv',dt.to_dict('records'))
print(phase.to_string(index=False));print(dt.to_string(index=False))
print(pd.read_csv(p/'tables/calibrated_parameters.csv').to_string(index=False))
plan=read(P/'provenance/C_frozen.json');fail=[];length=[]
for eid in plan['validation']:
 a=load(P/'data'/f'{eid}.npz');pr=load(P/'predictions'/f'{eid}_hydrographs.npz');active=GRID<=360
 length.append(dict(event=eid,MAE_full_km=float(np.mean(abs(a['b_true']-a['b_pred']))),MAE_0_360h_km=float(np.mean(abs(a['b_true'][active]-a['b_pred'][active]))),max_error_km=float(max(abs(a['b_true']-a['b_pred'])))))
 late=GRID>=plan['events'][eid]['D_h']
 for name in ['MC_oracle','MC_autonomous']:
  for x in [20,0]:
   j=np.argmin(abs(a['x']-x));q=pr[name][int((200-x)/20)];err=q-a['q'][:,j];idx=np.where(late)[0][np.argmax(abs(err[late]))]
   fail.append(dict(event=eid,model=name,x_from_outlet_km=x,distance_from_inlet_km=200-x,largest_late_error_m3s=float(err[idx]),time_h=float(GRID[idx]),truth_m3s=float(a['q'][idx,j]),prediction_m3s=float(q[idx]),true_b_km=float(a['b_true'][idx]),b_change_previous_h_km=float(a['b_true'][idx]-a['b_true'][max(idx-1,0)])))
frame('late_release_failure_cases.csv',fail);frame('endpoint_error_window_check.csv',length)
print(pd.DataFrame(fail).sort_values('largest_late_error_m3s',ascending=False).head(8).to_string(index=False))
