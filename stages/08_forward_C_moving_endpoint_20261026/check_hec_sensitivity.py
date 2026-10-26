from common import *
plan=read(P/'provenance/C_frozen.json');rows=[];ok=True
for label in ['dt30','dx1000']:
 for kind in ['ref','var']:
  original=f"pilot_s{round(plan['selected_slope']*1e6):03d}_CP01_{kind}"
  case=f'sensitivity_CP01_{label}_{kind}';a=load(C/'results'/original/'full_fields.npz');b=load(C/'results'/case/'full_fields.npz');qual=read(C/'results'/case/'quality.json')
  for x in TARGETS:
   ia=np.argmin(abs(a['x_m']-x*1000));ib=np.argmin(abs(b['x_m']-x*1000));qa=a['q'][:,ia];qb=b['q'][:,ib];za=a['z'][:,ia];zb=b['z'][:,ib]
   rmse=float(np.sqrt(np.mean((qb-qa)**2)));peak=float(abs(qb.max()-qa.max())/qa.max());stage=float(np.max(abs(zb-za)))
   passed=rmse/max(qa.max(),1)<.02 and peak<.02 and stage<.05 and all(qual[k] for k in ['basic_quality_pass','initialization_pass','tail_recovered']) and qb.min()>0
   rows.append(dict(case=case,kind=kind,x_from_outlet_km=int(x),distance_from_inlet_km=int(200-x),flow_RMSE_m3s=rmse,peak_relative_difference=peak,max_stage_difference_m=stage,passed=bool(passed)))
   ok &= passed
frame('HEC_sensitivity.csv',rows);save(P/'provenance/HEC_sensitivity_gate.json',dict(passed=bool(ok),criteria='flowRMSE/peak<2%; peak change<2%; max water-level change<0.05m; completequality+initialization+recovery+positiveflow',cases=4,created=hec.now()))
assert ok,'Numerical sensitivity gate failed; preserve and diagnose before batch'
print(pd.DataFrame(rows).groupby('case')[['flow_RMSE_m3s','peak_relative_difference','max_stage_difference_m']].max().to_string())
