from common import *
from models import rwec,muskingum
from moving_balance import audit,sample_flow
import time

def simulate(d,e,fit,b,slope,dt=None,alpha_override=None):
 dt=fit['dt_h'] if dt is None else dt;g=np.arange(0,552+dt/2,dt);u=hec.boundary(e,g,'ref')[0];bg=np.interp(g,GRID,b);th=fit['theta'];mid=fit['model_id']
 tic=time.perf_counter()
 if mid in (0,3,4):
  alpha=0. if mid==0 else th[3]
  if alpha_override is not None:alpha=alpha_override
  r=rwec(u,bg,dt,np.exp(th[0]),np.exp(th[1]),th[2],alpha,0. if mid==0 else th[4],slope,d['initial_depth'],d['initial_volumes'])
  q,ss,tt,bb,vol,ok,step,maxbal,maxlocal=r
 else:
  q,st,ok,step,maxbal=muskingum(u,dt,np.exp(th[0]),th[1],1. if mid==1 else th[2]);vol=d['initial_volumes'][:,None]+st-st[:,[0]];ss=st;tt=np.zeros_like(st);bb=np.zeros_like(st);maxlocal=maxbal
 runtime=time.perf_counter()-tic
 return dict(grid=g,dt=dt,q=q,vol=vol,ss=ss,tt=tt,bb=bb,ok=bool(ok),failstep=int(step),fixed_balance=float(maxbal),maxlocal=float(maxlocal),runtime=runtime,b=bg)

def main():
 plan=read(P/'provenance/C_frozen.json');fits=read(P/'provenance/frozen_parameters.json');rows=[];moving=[];balances=[];sensitivity=[];failures=[];ablation=[]
 assert len(fits)==5
 for eid in plan['order']:
  d=load(P/'data'/f'{eid}.npz');subset='train' if eid in plan['train'] else 'validation';casepred={};basepred={}
  for name,fit in fits.items():
   b=d['b_true'] if name=='MC_oracle' else d['b_pred'];run=simulate(d,plan['events'][eid],fit,b,plan['selected_slope'])
   if not run['ok']:
    failures.append(dict(event=eid,subset=subset,model=name,step=run['failstep'],time_h=run['grid'][run['failstep']],reason='Invalid storage/negative discharge/nonconvergence; no clipping or fallback'));qh=np.full((11,553),np.nan)
   else:qh=np.array([np.interp(GRID,run['grid'],a) for a in run['q']])
   casepred[name]=qh
   for idx,x in enumerate(TARGETS):
    truth=d['truth'][idx];pred=qh[NODES[idx]];dq=truth-np.interp(x,d['x'],d['qref'][0])*0 # overwritten with matched fixed-position difference
    j=np.argmin(abs(d['x']-x));dq=d['q'][:,j]-d['qref'][:,j]
    for win,mask in [('full',np.ones(553,dtype=bool)),('active_0_360h',GRID<=360),('suppression',dq<=-50),('release',dq>=50)]:
     met=metric(truth[mask],pred[mask],GRID[mask]) if np.sum(mask)>=3 else dict(NSE=np.nan,RMSE=np.nan,KGE=np.nan,peak_error_pct=np.nan,peak_time_error_h=np.nan)
     rows.append(dict(event=eid,subset=subset,model=name,x_from_outlet_km=int(x),distance_from_inlet_km=int(200-x),window=win,n=int(mask.sum()),valid=run['ok'],**met))
   # All models evaluated at identical true and predicted endpoint locations; truth never feeds routing.
   for pos,position in [('true_common',d['b_true']),('predicted_common',d['b_pred'])]:
    py=sample_flow(qh,position);hy=np.array([np.interp(bb,d['x'],qq) for bb,qq in zip(position,d['q'])])
    moving.append(dict(event=eid,subset=subset,model=name,position=pos,**metric(hy,py,GRID)))
   moving.append(dict(event=eid,subset=subset,model=name,position='end_to_end_predicted_vs_true',**metric(d['qb_true'],sample_flow(qh,d['b_pred']),GRID)))
   if run['ok']:
    bal=audit(run['q'],run['vol'],run['b'],run['dt']);cum=cumulative_trapezoid(d['area'],d['x']*1000,axis=1,initial=0)
    hecsm=np.array([v[-1]-np.interp(bb,d['x'],v) for bb,v in zip(b,cum)]);sm=np.interp(GRID,run['grid'],bal['sm'])
    balances.append(dict(event=eid,subset=subset,model=name,dt_h=fit['dt_h'],fixed_max_residual_m3=run['fixed_balance'],moving_max_residual_m3=bal['max_balance_m3'],moving_RMS_m3s=bal['rms_balance_m3s'],no_sweep_RMS_m3s=bal['rms_no_sweep_m3s'],sampled_flux_RMS_m3s=bal['rms_sampled_balance_m3s'],storage_RMSE_million_m3=float(np.sqrt(np.mean((sm-hecsm)**2))/1e6),runtime_s=run['runtime'],jump_steps=bal['jump_steps_over20km']))
    np.savez_compressed(P/'predictions'/f'{eid}_{name}_states.npz',time=run['grid'],q=run['q'],vol=run['vol'],storage_RWEC=run['ss'],storage_transit=run['tt'],storage_backwater=run['bb'],b=run['b'],hec_storage_at_same_position=hecsm,**bal)
   if name.startswith('MC_'):
    off=simulate(d,plan['events'][eid],fit,b,plan['selected_slope'],alpha_override=0.)
    qo=np.array([np.interp(GRID,off['grid'],a) for a in off['q']]) if off['ok'] else np.full_like(qh,np.nan)
    for j,x in enumerate(TARGETS):
     ablation.append(dict(event=eid,subset=subset,model=name,x_from_outlet_km=int(x),distance_from_inlet_km=int(200-x),NSE_with=metric(d['truth'][j],qh[NODES[j]])['NSE'],NSE_without=metric(d['truth'][j],qo[NODES[j]])['NSE']))
   if subset=='validation':
    for factor in [.5,2.]:
     other=simulate(d,plan['events'][eid],fit,b,plan['selected_slope'],dt=fit['dt_h']*factor)
     qp=np.array([np.interp(GRID,other['grid'],a) for a in other['q']]) if other['ok'] else np.full_like(qh,np.nan)
     for j,x in enumerate(TARGETS):sensitivity.append(dict(event=eid,model=name,dt_factor=factor,dt_h=fit['dt_h']*factor,x_from_outlet_km=int(x),distance_from_inlet_km=int(200-x),valid=other['ok'],**metric(d['truth'][j],qp[NODES[j]])))
  np.savez_compressed(P/'predictions'/f'{eid}_hydrographs.npz',time=GRID,**casepred)
  print(eid,subset,'evaluated',flush=True)
 frame('fixed_section_event_metrics.csv',rows);frame('moving_endpoint_event_metrics.csv',moving);frame('water_balance_event_metrics.csv',balances);frame('model_dt_sensitivity.csv',sensitivity);frame('prediction_failures.csv',failures);frame('fixed_parameter_storage_ablation.csv',ablation)
 df=pd.DataFrame(rows);full=df[df.window=='full'];agg=full.groupby(['subset','model','x_from_outlet_km','distance_from_inlet_km']).agg(events=('event','size'),valid_events=('valid','sum'),NSE=('NSE','mean'),RMSE=('RMSE','mean'),KGE=('KGE','mean'),peak_error_pct=('peak_error_pct','mean'),peak_time_error_h=('peak_time_error_h','mean')).reset_index()
 # Do not silently turn failure omission into a complete mean.
 agg.loc[agg.valid_events<agg.events,['NSE','RMSE','KGE','peak_error_pct','peak_time_error_h']]=np.nan
 frame('fixed_section_mean_metrics.csv',agg.to_dict('records'))
 pivot=full.pivot_table(index=['event','subset','x_from_outlet_km','distance_from_inlet_km'],columns='model',values='NSE',dropna=False).dropna(how='all').reset_index()
 for name in ['MC_oracle','MC_autonomous','NLM','LM']:pivot[f'{name}_minus_RWEC0']=pivot[name]-pivot['RWEC0']
 frame('event_raw_improvements.csv',pivot.to_dict('records'))
 print(agg[agg.subset=='validation'].to_string(index=False))

if __name__=='__main__':main()
