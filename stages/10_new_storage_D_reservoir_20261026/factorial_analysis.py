from common10 import *
from scipy.integrate import cumulative_trapezoid

def main():
 cfg=read(D/'config/D_frozen.json');rows=[];storage=[];ranges=[]
 for e in cfg['events']:
  eid=e['event_id'];ds={}
  for domain,geom,bc in [('DUF','uniform','fixed'),('DUR','uniform','reservoir'),('DCF','complex','fixed'),('DCR','complex','reservoir')]:
   f=D/'models'/f"{eid}_{geom}_{bc}_{cfg['version']}";d=load(f/'extracted.npz');keep=d['time']>=0
   ds[domain]={k:(v[keep] if v.ndim and len(v)==len(keep) else v) for k,v in d.items()}
   if bc=='reservoir':
    pp=load(f/'paired.npz');ranges.append(dict(event=eid,subset=e['development_split'],geometry=geom,min_native_Q=float(d['minflow_native'].min()),max_reservoir_rise=float(d['reservoir_H'].max()-5),max_endpoint_km=float(pp['b_true'].max()),max_extra_volume_m3=float(pp['extra_volume'].max()),max_detached_components=int(pp['islands'].max())))
  for x in TARGETS:
   idx=int(np.argmin(abs(ds['DUF']['x']-x)));q={k:v['q'][:,idx] for k,v in ds.items()};components={'geometry_at_fixed':q['DCF']-q['DUF'],'backwater_uniform':q['DUR']-q['DUF'],'backwater_complex':q['DCR']-q['DCF'],'interaction':q['DCR']-q['DCF']-q['DUR']+q['DUF']}
   for key,delta in components.items():rows.append(dict(event=eid,subset=e['development_split'],x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),component=key,RMS_m3s=float(np.sqrt(np.mean(delta**2))),min_difference_m3s=float(delta.min()),max_difference_m3s=float(delta.max()),peak_difference_m3s=float((q['DCR'].max()-q['DCF'].max()) if key=='backwater_complex' else np.nan)))
  vs={k:v['total_storage']-v['total_storage'][0] for k,v in ds.items()};components={'geometry_at_fixed':vs['DCF']-vs['DUF'],'backwater_uniform':vs['DUR']-vs['DUF'],'backwater_complex':vs['DCR']-vs['DCF'],'interaction':vs['DCR']-vs['DCF']-vs['DUR']+vs['DUF']}
  for key,dv in components.items():storage.append(dict(event=eid,subset=e['development_split'],component=key,volume_RMS_m3=float(np.sqrt(np.mean(dv**2))),volume_max_m3=float(dv.max()),volume_min_m3=float(dv.min()),storage_rate_RMS_m3s=float(np.sqrt(np.mean((np.diff(dv)/3600)**2)))))
 frame('D_factorial_physical_discharge.csv',rows);frame('D_factorial_physical_storage.csv',storage);frame('D_dynamic_ranges.csv',ranges)
 print(pd.DataFrame(rows).query('x_from_outlet_km==20').groupby('component')[['RMS_m3s','min_difference_m3s','max_difference_m3s']].mean().to_string())
if __name__=='__main__':main()
