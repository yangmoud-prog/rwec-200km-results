from build_D import *
import re

def extract(folder):
 cfg=json.loads((folder/'settings.json').read_text(encoding='utf-8'));design=json.loads((folder/'D_design.json').read_text(encoding='utf-8'))
 with h5py.File(folder/f'{NAME}.p01.hdf') as f:
  ts=f[TS];rawt=ts['Time'][:];elapsed=(rawt.astype(float)-float(rawt[0]))*24;time_jitter=float(np.max(np.abs(elapsed-np.rint(elapsed))))
  if time_jitter>0.001:raise ValueError('HEC output is not the frozen hourly grid')
  t=np.rint(elapsed)-ORIGIN_H
  if not np.array_equal(np.diff(t),np.ones(len(t)-1)):raise ValueError('Missing or duplicate hourly HEC output')
  attrs=f['Geometry/Cross Sections/Attributes'][:];xx=np.array([float(v.decode()) for v in attrs['RS']]);order=np.argsort(xx);x=xx[order]/1000
  q=ts['Cross Sections/Flow'][:].astype(float)[:,order];z=ts['Cross Sections/Water Surface'][:].astype(float)[:,order];vel=ts['Cross Sections/Velocity Total'][:].astype(float)[:,order];cumq=ts['Cross Sections/Flow Volume Cumulative'][:].astype(float)[:,order]
  pt=f['Geometry/Cross Sections/Property Tables'];info=pt['XSEC Info'][:];vals=pt['XSEC Value'][:];area=np.zeros_like(q);coverage=True
  for col,j in enumerate(order):
   st,n,_=info[j];v=vals[st:st+n];aa=v[:,1:4].sum(axis=1);area[:,col]=np.interp(z[:,col],v[:,0],aa);coverage=coverage and z[:,col].min()>=v[0,0] and z[:,col].max()<=v[-1,0]
  summary=f['Results/Unsteady/Output/Output Blocks/Base Output/Summary Output/Cross Sections'];minflow=summary['Minimum Flow'][:][0,order]
  if design['boundary']=='reservoir':
   raw=ts['Storage Areas/Reservoir/Storage Area Variables'][:].astype(float);tail=ts['Storage Areas/TailReceiver/Storage Area Variables'][:].astype(float)
   hr=raw[:,0];release=raw[:,3];vr=design['reservoir_area_m2']*hr
  else:raw=np.empty((len(t),0));tail=raw.copy();hr=np.full(len(t),5.);release=q[:,0];vr=np.zeros(len(t))
 s=np.trapz(area,x*1000,axis=1);boundary=np.genfromtxt(folder/'boundary.csv',delimiter=',',names=True,encoding='utf-8-sig');inp=np.interp(t,boundary['t_since_flood_h'],boundary['inflow_m3s'])
 d=dict(time=t,x=x,q=q,z=z,area=area,velocity=vel,cumulative_q=cumq,total_storage=s,reservoir_H=hr,reservoir_V=vr,release=release,reservoir_raw=raw,tail_raw=tail,inflow=inp,minflow_native=minflow)
 np.savez_compressed(folder/'extracted.npz',**d)
 pre=t<=0;tailmask=t>=t.max()-24;flowbalance=np.diff(s)/3600-.5*(q[1:,-1]+q[:-1,-1]-q[1:,0]-q[:-1,0]);delta_s=s-s[0];net=cumulative_trapezoid(q[:,-1]-q[:,0],t*3600,initial=0)
 textlog=(folder/f'{NAME}.p01.computeMsgs.txt').read_text(errors='replace');m=re.search(r'Overall Volume Accounting Error as percentage:\s*([-+\deE.]+)',textlog);engine_error=float(m.group(1)) if m else None
 qual=dict(case=folder.name,geometry=design['geometry'],boundary=design['boundary'],complete=bool(abs(t[-1]-cfg['end_since_flood_h'])<.01),finite=bool(np.all(np.isfinite(q)) and np.all(np.isfinite(z))),minimum_native_Q_m3s=float(minflow.min()),initial_Q_max_error_m3s=float(abs(q[pre]-1000).max()),initial_H_error_m=float(abs(hr[pre]-5).max()),tail_Q_max_error_m3s=float(abs(q[tailmask]-1000).max()),H_rise_m=float(hr.max()-5),HTab_covered=bool(coverage),HEC_volume_error_pct=engine_error,independent_river_storage_balance_RMSE_m3s=float(np.sqrt(np.mean(flowbalance**2))),independent_river_volume_residual_fraction=float(abs(delta_s-net).max()/max(1,np.ptp(s))),native_summary_covers_computation_step=True)
 if design['boundary']=='reservoir':
  residual=np.diff(vr)/3600-.5*(raw[1:,1]+raw[:-1,1]);qual.update(reservoir_balance_RMSE_m3s=float(np.sqrt(np.mean(residual**2))),reservoir_area_from_net_flux_m2=float(np.dot(hr-hr[0],cumulative_trapezoid(raw[:,1],dx=3600,initial=0))/np.dot(hr-hr[0],hr-hr[0])),tail_max_H_m=float(tail[:,0].max()),tail_unsubmerged=bool(tail[:,0].max()<design['weir_crest_m']),raw_SA_output_unit_factor=float(np.median(vr/raw[:,5])),release_rule_RMSE_m3s=float(np.sqrt(np.mean((release-design['release_C']*design['weir_length_m']*np.maximum(hr-design['weir_crest_m'],0)**1.5)**2))))
 qual['maximum_raw_hourly_time_jitter_h']=time_jitter;qual['time_axis_method']='round verified hourly elapsed HEC output before applying 48 h origin'
 qual['numerical_gate']=bool(qual['complete'] and qual['finite'] and qual['HTab_covered'] and minflow.min()>0 and qual['initial_Q_max_error_m3s']<1 and qual['tail_Q_max_error_m3s']<1 and engine_error is not None and abs(engine_error)<.1 and qual['independent_river_volume_residual_fraction']<.03 and qual.get('tail_unsubmerged',True) and qual.get('reservoir_balance_RMSE_m3s',0)<2)
 dump(folder/'quality.json',qual);return d,qual

def pair(dynamic,reference):
 d=np.load(dynamic/'extracted.npz');r=np.load(reference/'extracted.npz')
 if not np.array_equal(d['time'],r['time']) or not np.array_equal(d['x'],r['x']):raise ValueError('Paired HEC grids differ')
 dh=d['z']-r['z'];x=d['x'];ends=[];conn=[];islands=[]
 for row in dh:
  mask=row>=.1;ids=np.flatnonzero(mask);b=0.;bc=0.
  if len(ids):
   j=ids[-1];b=x[j] if j==len(x)-1 else x[j]+(.1-row[j])*(x[j+1]-x[j])/(row[j+1]-row[j])
  if mask[0]:
   js=np.flatnonzero(~mask);j=js[0]-1 if len(js) else len(x)-1;bc=x[j] if j==len(x)-1 else x[j]+(.1-row[j])*(x[j+1]-x[j])/(row[j+1]-row[j])
  ends.append(b);conn.append(bc);islands.append(int(mask[0])+int(np.sum(np.diff(mask.astype(int))==1)))
 ends=np.array(ends);dq=d['q']-r['q'];dv=np.trapz(d['area']-r['area'],x*1000,axis=1);indices=[int(np.argmin(abs(x-k))) for k in [20,40,60,80]]
 result=dict(event=dynamic.name.split('_')[0],geometry=json.loads((dynamic/'D_design.json').read_text())['geometry'],endpoint_max_km=float(ends.max()),endpoint_range_km=float(np.ptp(ends)),detached_times=int(np.sum(np.array(islands)>1)),minimum_Q_m3s=float(d['minflow_native'].min()),max_suppression_20_80_m3s=float(-dq[:,indices].min()),max_recovery_20_80_m3s=float(dq[:,indices].max()),extra_storage_max_m3=float(dv.max()),extra_storage_min_m3=float(dv.min()),storage_increasing_and_releasing=bool(np.diff(dv).max()>1e5 and np.diff(dv).min()<-1e5))
 result['physical_gate']=bool(result['minimum_Q_m3s']>0 and result['endpoint_range_km']>20 and result['max_suppression_20_80_m3s']>100 and result['max_recovery_20_80_m3s']>100 and result['storage_increasing_and_releasing'])
 np.savez_compressed(dynamic/'paired.npz',time=d['time'],x=x,deltaH=dh,deltaQ=dq,extra_volume=dv,b_true=ends,b_connected=conn,islands=islands)
 dump(dynamic/'pair_quality.json',result);return result

if __name__=='__main__':
 for case in sys.argv[1:]:
  _,qual=extract(ROOT/'models'/case);print(json.dumps(qual),flush=True)
