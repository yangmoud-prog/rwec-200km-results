"""Falsification of two zero-fit profile candidates with true endpoint and true reference.
This is an information-rich diagnostic, never autonomous forecast input.
"""
from common10 import *
from scipy.integrate import cumulative_trapezoid
import h5py

def area_exact(H,sections):
 out=np.zeros_like(H)
 for j,s in enumerate(sections):
  pts=s['points'];h=H[:,j]
  for (x0,z0),(x1,z1) in zip(pts[:-1],pts[1:]):
   d0=np.maximum(h-z0,0);d1=np.maximum(h-z1,0);full=h>=max(z0,z1);out[:,j]+=np.where(full,(x1-x0)*(d0+d1)/2,(x1-x0)*np.maximum(d0,d1)**2/(2*max(abs(z1-z0),1e-20)))
 return out

def atb(a,x,b):return np.array([np.interp(v,x,row) for v,row in zip(b,a)])
def moving_terms(a,x,b):
 cum=cumulative_trapezoid(a,x*1000,axis=1,initial=0);sb=atb(cum,x,b);sweep=[]
 for i in range(1,len(b)):
  l,h=sorted([b[i-1],b[i]]);xx=np.r_[l,x[(x>l)&(x<h)],h]
  sweep.append(np.sign(b[i]-b[i-1])*np.trapz(np.interp(xx,x,.5*(a[i]+a[i-1])),xx*1000))
 return sb,np.array(sweep)

def evaluate(folder,ref):
 d=load(folder/'extracted.npz');r=load(ref/'extracted.npz');paired=load(folder/'paired.npz');design=read(folder/'D_design.json');sections=sorted(design['sections'],key=lambda v:v['x']);x=d['x'];keep=d['time']>=0
 for a in [d,r,paired]:
  for k,v in a.items():
   if isinstance(v,np.ndarray) and v.ndim and len(v)==len(keep):a[k]=v[keep]
 b=paired['b_true'];dh=d['z']-r['z'];h0=dh[:,0];truea=area_exact(d['z'],sections);referencea=area_exact(r['z'],sections)
 true_sb,true_sweep=moving_terms(truea,x,b);qb=atb(d['q'],x,b);q0=d['q'][:,0];net=.5*(qb[1:]+qb[:-1]-q0[1:]-q0[:-1]);hec_balance=np.diff(true_sb)/3600-true_sweep/3600-net
 # The reference has precisely the same time-dependent fixed-boundary flood.
 # Second candidate uses friction-gradient geometry; oracle fields give the most favourable physical diagnostic.
 with h5py.File(ref/'Backwater200.p01.hdf') as f:
  tab=f['Geometry/Cross Sections/Property Tables'];info=tab['XSEC Info'][:];values=tab['XSEC Value'][:];decay=np.zeros_like(dh)
  for j in range(len(x)):
   idx=len(x)-1-j;st,n,_=info[idx];v=values[st:st+n];z=v[:,0];K=v[:,7:10].sum(axis=1);DK=np.gradient(K,z);aa=v[:,1:4].sum(axis=1);T=v[:,16]
   kk=np.maximum(np.interp(r['z'][:,j],z,K),1e-9);dkk=np.maximum(np.interp(r['z'][:,j],z,DK),1e-9);ar=np.maximum(np.interp(r['z'][:,j],z,aa),1e-9);tw=np.maximum(np.interp(r['z'][:,j],z,T),1e-9);fr2=r['q'][:,j]**2*tw/(9.81*ar**3)
   decay[:,j]=2*r['q'][:,j]**2*dkk/(kk**3*np.maximum(1-fr2,.01))
 integrated=cumulative_trapezoid(decay,x*1000,axis=1,initial=0)
 rows=[]
 for candidate in ['exponential','hydraulic_coordinate']:
  delta=np.zeros_like(dh)
  for i in range(len(b)):
   if b[i]<=0 or h0[i]<=.1:continue
   coordinate=x/b[i] if candidate=='exponential' else integrated[i]/max(np.interp(b[i],x,integrated[i]),1e-12)
   delta[i]=h0[i]*np.exp(-np.log(h0[i]/.1)*coordinate)
  shapez=r['z']+delta;shapea=area_exact(shapez,sections);sb,sweep=moving_terms(shapea,x,b)
  # Independent physical flow prediction from shape storage, real outlet flow and exact same moving boundary.
  qb_mean=q0[1:]/2+q0[:-1]/2+np.diff(sb)/3600-sweep/3600;qb_true_mean=.5*(qb[1:]+qb[:-1]);err=qb_mean-qb_true_mean
  active=(b>0)[:,None]&(x[None,:]<=b[:,None]);extra=np.trapz(shapea-referencea,x*1000,axis=1);extra_true=np.trapz(truea-referencea,x*1000,axis=1)
  stage_rmse=float(np.sqrt(np.mean((shapez[active]-d['z'][active])**2))) if active.any() else 0.;vol_error=float(np.sqrt(np.mean((extra-extra_true)**2))/max(np.ptp(extra_true),1));rate_error=float(np.sqrt(np.mean(err**2)));flow_nse=float(metrics(qb_true_mean,qb_mean)['NSE']);hecerr=float(np.sqrt(np.mean(hec_balance**2)))
  row=dict(event=folder.name.split('_')[0],geometry=design['geometry'],candidate=candidate,reference='true fixed-boundary unsteady surface; diagnostic only',endpoint='true HEC threshold location; diagnostic only',stage_active_RMSE_m=stage_rmse,extra_volume_NRMSE_range=vol_error,storage_rate_flux_RMSE_m3s=rate_error,true_HEC_balance_RMSE_m3s=hecerr,endpoint_flow_NSE=flow_nse,endpoint_flow_RMSE_m3s=rate_error,exact_area_vs_HEC_HTab_RMSE_m2=float(np.sqrt(np.mean((truea-d['area'])**2))),physical_gate=bool(stage_rmse<=.1 and vol_error<=.05 and rate_error<=max(10,5*hecerr) and flow_nse>=.95))
  regular=(abs(np.diff(b))<=2)&(b[1:]>0)&(b[:-1]>0)
  row.update(endpoint_jump_over_2km_per_hour_count=int(np.sum(abs(np.diff(b))>2)),regular_active_hours=int(regular.sum()),regular_active_flux_RMSE_m3s=float(np.sqrt(np.mean(err[regular]**2))) if regular.any() else None,regular_HEC_flux_RMSE_m3s=float(np.sqrt(np.mean(hec_balance[regular]**2))) if regular.any() else None)
  rows.append(row);np.savez_compressed(P/'predictions'/f"oracle_{folder.name}_{candidate}.npz",time=d['time'],x=x,b=b,z_true=d['z'],z_shape=shapez,delta_true=dh,delta_shape=delta,SB_true=true_sb,SB_shape=sb,extra_volume_true=extra_true,extra_volume_shape=extra,qb_true=qb_true_mean,qb_shape=qb_mean,sweep_true=true_sweep,sweep_shape=sweep,HEC_balance=hec_balance)
 return rows

def main():
 protocol=dict(created_before_oracle_results=True,fit_parameters=0,candidates=['exponential','hydraulic_coordinate'],thresholds=dict(stage_active_RMSE_m=.1,extra_volume_NRMSE_range=.05,flux_RMSE_m3s='max(10,5*HEC hourly residual)',endpoint_flow_NSE=.95),interpretation='A pass only removes this shape falsification; it does not establish autonomous forecast skill. A fail blocks using this shape as a joint-state geometric closure.',future_truth='only diagnostic upper-bound oracle; no fitted model is driven by these inputs')
 save(P/'provenance/oracle_closure_protocol.json',protocol);rows=[]
 for geom in ['uniform','complex']:
  for event in ['DP01','DP02','DP03','DP04']:
   folder=D/'models'/f'{event}_{geom}_reservoir_v3';ref=D/'models'/f'{event}_{geom}_fixed_v2'
   if folder.with_name(folder.name+'_full').exists():folder=folder.with_name(folder.name+'_full');ref=ref.with_name(ref.name+'_full')
   if not (folder/'paired.npz').exists():continue
   rows+=evaluate(folder,ref)
 frame('D_oracle_physical_closure.csv',rows)
 if rows:print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__':main()
