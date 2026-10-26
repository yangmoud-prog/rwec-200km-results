from common import *

def main():
 c=W/'variants/forward_dynamic_C_20261026/results'
 p1=c/'counterfactual_CP01_1m_var/full_fields.npz';p4=c/'pilot_s075_CP01_var/full_fields.npz';pr=c/'pilot_s075_CP01_ref/full_fields.npz'
 one=load(p1);four=load(p4);ref=load(pr);t=one['time'];use=t>=-1e-6;x=one['x_m'][::-1]/1000;q1=one['q'][use,::-1];q4=four['q'][use,::-1];z1=one['z'][use,::-1];z4=four['z'][use,::-1];zr=ref['z'][use,::-1];b1=np.array([extent(v,x)[0] for v in z1-zr]);b4=np.array([extent(v,x)[0] for v in z4-zr]);rows=[]
 for xx in TARGETS:
  j=int(xx/2);difference=q4[:,j]-q1[:,j]
  rows.append(dict(x_from_outlet_km=int(xx),distance_from_inlet_km=200-int(xx),equal_probability_RMSE_lower_bound_m3s=float(np.sqrt(np.mean(difference**2))/2),maximum_unavoidable_absolute_error_m3s=float(np.max(abs(difference))/2),mean_two_scenario_envelope_width_m3s=float(np.mean(abs(difference))),max_envelope_width_m3s=float(np.max(abs(difference)))))
 frame('C_identifiability_two_scenario_bounds.csv',rows)
 np.savez_compressed(P/'predictions/C_counterfactual.npz',time=t[use],x=x,q1=q1,q4=q4,b1=b1,b4=b4,h1=z1[:,0],h4=z4[:,0])
 cfg=read(P8/'provenance/C_frozen.json');unc=[]
 for e in cfg['validation']:
  d=load(P8/'data'/f'{e}.npz');active=d['b_true']>.1;inside=(d['b_true']>=d['b_low'])&(d['b_true']<=d['b_high'])
  unc.append(dict(event=e,horizon_h=552,method='existing external8NN empirical10-90percentile, not joint model',full_coverage=float(inside.mean()),active_coverage=float(inside[active].mean()),mean_width_km=float(np.mean(d['b_high']-d['b_low'])),endpoint_MAE_km=float(np.mean(abs(d['b_pred']-d['b_true'])))))
 frame('C_external_predictor_uncertainty_reference.csv',unc)
 zero=int(np.argmin(abs(t)))
 save(P/'provenance/C_identifiability.json',dict(source_sha256=[sha(p) for p in [p1,p4,pr]],same_past_stage_max_difference_m=float(np.max(abs(one['z'][t<=1e-6,-1]-four['z'][t<=1e-6,-1]))),same_initial_state_stage_max_difference_m=float(np.max(abs(one['z'][zero]-four['z'][zero]))),same_upstream_HEC_flow_max_difference_m3s=float(np.max(abs(one['q'][:,0]-four['q'][:,0]))),endpoint_equal_probability_RMSE_lower_bound_km=float(np.sqrt(np.mean((b4-b1)**2))/2),meaning='two physically admissible unknown future-stage scenarios, no assumed real probability or calibrated confidence interval',C_joint_deterministic_prediction_claimed=False,C_data_used_to_fit_B=False))
 print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(unc).mean(numeric_only=True).to_string())
if __name__=='__main__':main()
