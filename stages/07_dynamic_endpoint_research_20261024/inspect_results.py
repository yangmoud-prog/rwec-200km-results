from common import *

tables=['endpoint_events','mechanism_events','moving_boundary_balance','volume_proxy_metrics','routing_quality','fixed_parameter_timestep','indistinguishable_pair_bounds','paired_improvement','flow_validation_summary']
for name in tables:
 d=pd.read_csv(OUT/'tables'/f'{name}.csv');print('\n',name,len(d))
 if name=='endpoint_events':print(d[d.threshold_m==.1].groupby('group')[['maximum_length_km','mean_length_km','island_hours','censored_hours','zero_outlet_but_nonzero_length_hours']].mean().to_string())
 if name=='mechanism_events':print(d.groupby(['group','xc_km'])[['geometry_derivative_correlation','geometry_derivative_vs_deltaQ_RMSE']].mean().to_string())
 if name=='moving_boundary_balance':print(d.groupby('group').mean(numeric_only=True).to_string())
 if name=='volume_proxy_metrics':print(d[d.subset=='validation'].groupby(['group','diagnostic'])[['NSE','RMSE']].mean().to_string())
 if name=='routing_quality':print(d.groupby(['group','model'])[['water_balance_max_m3','volume_limits','minimum_transit_storage_m3']].max().to_string())
 if name=='fixed_parameter_timestep':print(d.groupby(['group','model','dt_h'])[['NSE','RMSE','RMSE_vs_hourly']].mean().to_string())
 if name=='indistinguishable_pair_bounds':print(d.groupby('xc_km')[['common_prediction_RMSE_lower_bound_m3s','max_flow_gap_m3s']].mean().to_string())
 if name=='paired_improvement':print(d[d.model=='RWEC-L-predicted'].to_string(index=False))
 if name=='flow_validation_summary':print(d.to_string(index=False))
for g in ['A','B']:
 for m in ['RWEC-0','RWEC-L-predicted','RWEC-L-oracle','Muskingum-linear']:
  f=read(OUT/'fits'/f'{g}_{m}.json');print(g,m,f['parameters'],'seedloss',[r['loss'] for r in f['runs']],'bounds',f['bound_hits'])
print('Bvalidation',[m['event'] for m in read(OUT/'data/metadata.json') if m['group']=='B' and m['subset']=='validation'])
