from common import *
import h5py
from scipy.integrate import cumulative_trapezoid
from scipy.signal import correlate,correlation_lags

def setup():
 dirs();split=read(P4/'provenance/event_split.json');frozen=read(B/'config/pilot_frozen.json')
 trainA=set(split['train']);devA=trainA|set(split['validation'])
 plan=read(A/'config/formal_plan.json');pairs={}
 for c in plan['formal_cases']:
  if c['role']=='main' and c['event_id'] in devA:pairs.setdefault(c['event_id'],{})[c['kind']]=c
 manifest=[];meta=[];endpoints=[];mechanism=[];profiles=[];moving=[]
 for label,ids,train in [('A',sorted(devA),trainA),('B',frozen['order'],set(frozen['train']))]:
  for eid in ids:
   if label=='A':
    path=P4/'data'/f'A_{eid}.npz';d=load(path);t=d['time_h'];x=d['full_x_km'];dz=d['delta_stage_m'];q=d['flow_var_m3s'];qr=d['flow_ref_m3s'];qgrid=np.arange(0,201,10.)
    v=load(A/'results'/pairs[eid]['var']['case']/'reach_timeseries.npz');h=v['stage_m'][:,0]
    rpath=A/'models'/pairs[eid]['ref']['case']/'Backwater200.p01.hdf'
    with h5py.File(rpath) as f:zref=f['Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Cross Sections/Water Surface'][:].astype(float)[:,::-1]
    z=zref+dz;slope=.00005
    inf=load(P5/'predictions'/f'A_{eid}_RWEC-0.npz')['q'][:len(t),0]
   else:
    path=P5/'data'/f'B_{eid}.npz';d=load(path);t=d['time'];x=d['x'];dz=d['delta_stage'];q=d['q_control'];qr=d['q_ref'];qgrid=x;h=d['z_control'][:,0];z=d['z_control'];zref=d['z_ref'];slope=.000075;inf=d['inflow']
   manifest.append(dict(group=label,event=eid,path=str(path),sha256=sha(path)))
   assert t[0]<=-48 and t[-1]>=240 and np.allclose(np.diff(t),1,atol=1e-5)
   y=z-x[None,:]*1000*slope;yr=zref-x[None,:]*1000*slope
   area=500*y+3*y*y;ar=500*yr+3*yr*yr;dv=cumulative_trapezoid(area-ar,x*1000,axis=1,initial=0)
   vols=cumulative_trapezoid(area,x*1000,axis=1,initial=0);lbsets={}
   subset='train' if eid in train else 'validation'
   for threshold in [.05,.1,.2]:
    vals=[extent(v,x,threshold) for v in dz];ll=np.array([v[0] for v in vals]);lbsets[threshold]=ll
    f=(t>=0)&(t<=240)
    endpoints.append(dict(group=label,event=eid,subset=subset,threshold_m=threshold,maximum_length_km=float(ll[f].max()),mean_length_km=float(ll[f].mean()),censored_hours=int(sum(v[1] for v in np.array(vals,dtype=object)[f])),island_hours=int(sum(v[2]>1 for v in np.array(vals,dtype=object)[f])),zero_outlet_but_nonzero_length_hours=int(np.sum((dz[f,0]<=threshold)&(ll[f]>0)))))
   lb=lbsets[.1];connected=np.array([extent(v,x,.1)[3] for v in dz])
   for xc in TARGETS:
    ix=int(np.argmin(abs(x-xc)));jq=int(np.argmin(abs(qgrid-xc)));f=(t>=0)&(t<=240)
    vv=dv[:,-1]-dv[:,ix];dq=q[:,jq]-qr[:,jq];dvdt=np.gradient(vv,t*3600)
    leff=np.maximum(0,lb-xc)*1000;rect=500*np.maximum(h-5,0)*leff
    corr=float(np.corrcoef(dq[f],-dvdt[f])[0,1]) if np.std(dvdt[f])>0 else np.nan
    mechanism.append(dict(group=label,event=eid,subset=subset,xc_km=int(xc),positive_event=bool(np.all(q[f,jq]>=0)),
     maximum_effective_length_km=float(leff[f].max()/1000),peak_delta_volume_m3=float(vv[f].max()),minimum_delta_volume_m3=float(vv[f].min()),
     suppression_m3s=float(dq[f].min()),recovery_m3s=float(dq[f].max()),geometry_derivative_vs_deltaQ_RMSE=float(np.sqrt(np.mean((dq[f]+dvdt[f])**2))),geometry_derivative_correlation=corr,
     peak_time_h=float(t[np.flatnonzero(f)[np.argmax(q[f,jq])]]),inflow_peak_time_h=float(t[np.argmax(inf)]),
     length_peak_time_h=float(t[np.flatnonzero(f)[np.argmax(lb[f])]]),volume_peak_time_h=float(t[np.flatnonzero(f)[np.argmax(vv[f])]])))
    for k in np.flatnonzero(f):profiles.append(dict(group=label,event=eid,subset=subset,time_h=float(t[k]),xc_km=int(xc),lb_km=float(lb[k]),leff_km=float(leff[k]/1000),delta_volume_m3=float(vv[k]),rectangle_proxy_m3=float(rect[k]),deltaQ_m3s=float(dq[k])))
   # Reynolds transport diagnostic on the free upstream portion [lb(t),200km].
   sm=np.array([vols[k,-1]-np.interp(lb[k],x,vols[k]) for k in range(len(t))])
   moving_term=np.zeros(len(t));qphys=np.zeros(len(t))
   for k in range(1,len(t)):
    avgarea=.5*(area[k]+area[k-1]);cum=cumulative_trapezoid(avgarea,x*1000,initial=0)
    sweep=np.interp(lb[k],x,cum)-np.interp(lb[k-1],x,cum)
    bm=.5*(lb[k]+lb[k-1]);qb=.5*(np.interp(bm,qgrid,q[k])+np.interp(bm,qgrid,q[k-1]))
    moving_term[k]=sweep/3600;qphys[k]=qb
   f=(t>0)&(t<=240);idsf=np.flatnonzero(f);storage_rate=np.diff(sm)/3600
   influx=.5*(inf[1:]+inf[:-1]);res=storage_rate-influx+qphys[1:]+moving_term[1:]
   m=f[1:];moving.append(dict(group=label,event=eid,subset=subset,moving_mass_residual_RMSE_m3s=float(np.sqrt(np.mean(res[m]**2))),
    omitted_sweep_residual_RMSE_m3s=float(np.sqrt(np.mean((res[m]-moving_term[1:][m])**2))),maximum_sweep_flux_m3s=float(abs(moving_term[f]).max()),
    negative_relative_flux_hours=int(np.sum((qphys[f]+moving_term[f])<0))))
   sel=np.array([np.argmin(abs(t-v)) for v in GRID]);assert np.allclose(t[sel],GRID,atol=1e-5)
   future=GRID>=0;past=GRID<=0
   truth=np.column_stack([q[sel,int(np.argmin(abs(qgrid-xc)))] for xc in TARGETS]);ref=np.column_stack([qr[sel,int(np.argmin(abs(qgrid-xc)))] for xc in TARGETS])
   np.savez_compressed(OUT/'data'/f'{label}_{eid}.npz',time=GRID,inflow=inf[sel],past_H=h[sel][past],true_lb=lb[sel],true_connected=connected[sel],
    truth=truth,reference=ref,delta_stage=dz[sel],x=x,z=z[sel],zref=zref[sel],dv=dv[sel],moving_sweep=moving_term[sel],moving_Q=qphys[sel],moving_S=sm[sel])
   meta.append(dict(group=label,event=eid,subset=subset,file=f'{label}_{eid}.npz'))
  print('Prepared',label,len(ids),flush=True)
 save(OUT/'data/metadata.json',meta);save(OUT/'provenance/input_manifest.json',dict(inputs=manifest,sealed_test_results_opened=0))
 frame('endpoint_events.csv',endpoints);frame('mechanism_events.csv',mechanism);frame('mechanism_process.csv',profiles);frame('moving_boundary_balance.csv',moving)
 # Train-only geometry proxy fits; validation diagnostics use no refitting.
 P=pd.DataFrame(profiles);rr=[];fits=[]
 for (label,xc),g in P.groupby(['group','xc_km']):
  tr=g[g.subset=='train'];weights=1/tr.groupby('event').time_h.transform('size').to_numpy();ytr=tr.delta_volume_m3.to_numpy()
  for model,col in [('length_only','leff_km'),('height_length','rectangle_proxy_m3')]:
   xx=tr[col].to_numpy();gain=float(max(0,np.sum(weights*xx*ytr)/max(1e-20,np.sum(weights*xx*xx))))
   fits.append(dict(group=label,xc_km=int(xc),diagnostic=model,gain=gain))
   for eid,ev in g.groupby('event'):rr.append(dict(group=label,xc_km=int(xc),event=eid,subset=ev.subset.iloc[0],diagnostic=model,**metric(ev.delta_volume_m3.to_numpy(),gain*ev[col].to_numpy())))
 frame('volume_proxy_parameters.csv',fits);frame('volume_proxy_metrics.csv',rr)
 save(OUT/'checkpoint.json',dict(stage='endpoint_and_mechanism_calculated',A_events=70,B_events=20,sealed_test_results_opened=0))
 print('Endpoint/geometry/moving-boundary diagnostics completed',flush=True)

if __name__=='__main__':setup()
