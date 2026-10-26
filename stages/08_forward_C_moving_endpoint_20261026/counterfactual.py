from common import *
import sys
plan=read(P/'provenance/C_frozen.json')
if sys.argv[1]=='prepare':
 original=read(C/'config/pilot_events.json')[0]
 e=hec.event(original['Qp_m3s'],original['D_h'],original['r'],1.,original['tb_ratio'],original['Tr_ratio'],original['Tf_ratio'],'CP01_ALT1',original['Th_ratio'])
 case='counterfactual_CP01_1m_var';hec.build(case,e,'var',end_h=552,slope=plan['selected_slope'])
 save(C/'config/counterfactual_queue.json',[dict(case=case,event='CP01_ALT1',kind='var',phase='identifiability')])
else:
 a=load(C/'results'/f"pilot_s{round(plan['selected_slope']*1e6):03d}_CP01_var/full_fields.npz");b=load(C/'results/counterfactual_CP01_1m_var/full_fields.npz')
 assert np.max(abs(a['q'][:,0]-b['q'][:,0]))<.02
 past=a['time']<=0;assert np.max(abs(a['z'][past]-b['z'][past]))<1.e-5
 use=a['time']>=-1.e-6;rows=[]
 for x in TARGETS:
  j=np.argmin(abs(a['x_m']-x*1000));delta=a['q'][use,j]-b['q'][use,j]
  rows.append(dict(x_from_outlet_km=int(x),distance_from_inlet_km=int(200-x),pair_RMSE_m3s=float(np.sqrt(np.mean(delta**2))),deterministic_equal_probability_min_RMSE_m3s=float(.5*np.sqrt(np.mean(delta**2))),max_difference_m3s=float(max(abs(delta))),max_unavoidable_absolute_error_m3s=float(.5*max(abs(delta)))))
 frame('identifiability_counterfactual.csv',rows)
 save(P/'provenance/identifiability_audit.json',dict(identical_upstream_forecast=True,identical_preorigin_full_state=True,stage_rises_m=[4,1],different_future_stage_only=True,minimum_equal_probability_MSE='||Q4-Q1||^2/(4N)',model_training_used=False))
 print(pd.DataFrame(rows).to_string(index=False))
