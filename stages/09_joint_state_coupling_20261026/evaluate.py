from common import *
from models import rwec,muskingum
from joint import simulate,shape,atx,length_scale
import time

def run(d,name,fit,dt=None):
 dt=fit['dt_h'] if dt is None else dt;theta=fit['theta'];grid=np.arange(0,288+dt/2,dt);inflow=np.interp(grid,GRID,d['inflow']);st=time.perf_counter();native=name.startswith('Joint') or name=='Feedback_off'
 if native:
  mode=1 if name=='Joint_oracle' else (2 if name=='Feedback_off' else 0);bs=np.interp(grid,GRID,d['b_true']) if mode==1 else np.zeros(len(grid))
  q,b,h,qb,sb,su,w,it,er,fail,bal=simulate(inflow,dt,*np.exp(theta),d['yinit'],bs,mode);extra=dict(b=b,h=h,qb=qb,sb=sb,su=su,sweep_rate=w,iterations=it,iteration_residual=er)
 elif name in ('RWEC0','Serial8'):
  bs=np.interp(grid,GRID,d['b_pred']) if name=='Serial8' else np.zeros(len(grid));q,r,t,b,v,ok,fail,bal,local=rwec(inflow,bs,dt,np.exp(theta[0]),np.exp(theta[1]),theta[2],theta[3] if name=='Serial8' else 0,theta[4] if name=='Serial8' else 0,.000075,d['initial_depth'],d['initial_volumes']);extra=dict(total_storage=v.sum(axis=0))
 else:
  q,v,ok,fail,bal=muskingum(inflow,dt,np.exp(theta[0]),theta[1],theta[2] if name=='NLM' else 1.);extra=dict(total_storage=v.sum(axis=0))
 elapsed=time.perf_counter()-st
 if not native:bal/=dt*3600
 if native:volume=np.cumsum(np.c_[np.zeros(11),q[:,1:]*dt],axis=1)
 else:volume=cumulative_trapezoid(q,grid,axis=1,initial=0)
 hp=np.array([np.diff(np.interp(GRID,grid,v)) for v in volume]);instant=np.array([np.interp(GRID,grid,row) for row in q])
 result=dict(time=grid,q=q,hourly_mean=hp,hourly_display=instant,**extra)
 return result,dict(failed_step=int(fail),failed_time_h=float(fail*dt) if fail else None,max_balance_flow_m3s=float(bal),runtime_s=elapsed,min_fixed_Q=float(np.nanmin(q)),negative_count=int(np.sum(q<0)),dt_h=dt)

def dynamic_q(qhour,b):
 return np.array([np.interp(v,np.arange(0,201,20),qhour[::-1,i]) for i,v in enumerate(b)])

def main():
 cfg=read(P/'provenance/frozen_protocol.json');fits=read(P/'PARAMETERS.json');fits['Feedback_off']=dict(fits['Joint']);rows=[];dynamic=[];state=[];quality=[];sens=[];physics=[]
 for e in cfg['order']:
  d=load(P/'data'/f'B_{e}.npz');subset='train' if e in cfg['train'] else 'validation';truth=.5*(d['truth'][:,1:]+d['truth'][:,:-1]);bmid=.5*(d['b_true'][1:]+d['b_true'][:-1]);hecmean=.5*(d['q'][1:]+d['q'][:-1]);dynamictruth=np.array([np.interp(v,d['x'],qq) for v,qq in zip(bmid,hecmean)])
  for name,fit in fits.items():
   pred,qual=run(d,name,fit);quality.append(dict(event=e,subset=subset,model=name,**qual));np.savez_compressed(P/'predictions'/f'B_{e}_{name}.npz',**pred)
   for j,x in enumerate(TARGETS):rows.append(dict(event=e,subset=subset,model=name,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),complete=qual['failed_step']==0,**metrics(truth[j],pred['hourly_mean'][j+4])))
   pq=dynamic_q(pred['hourly_mean'],bmid);dynamic.append(dict(event=e,subset=subset,model=name,location='same true endpoint midpoint; evaluation only',complete=qual['failed_step']==0,**metrics(dynamictruth,pq)))
   if 'b' in pred:
    bp=np.interp(GRID,pred['time'],pred['b']);sb=np.interp(GRID,pred['time'],pred['sb']);su=np.interp(GRID,pred['time'],pred['su']);hp=np.interp(GRID,pred['time'],pred['h']);qb=np.interp(GRID,pred['time'],pred['qb']);cum=cumulative_trapezoid(d['area'],d['x']*1000,axis=1,initial=0);hecsb=np.array([np.interp(b,d['x'],v) for b,v in zip(bp,cum)]);hecsu=cum[:,-1]-hecsb
    active=d['b_true']>.1
    boundary_error=np.array([(hh-5)*np.exp(-bb*1000/length_scale(qq,np.exp(fit['theta'][2])))-.1 for bb,hh,qq in zip(bp,hp,qb) if bb>0 and np.isfinite(qq)])
    state.append(dict(event=e,subset=subset,model=name,complete=qual['failed_step']==0,endpoint_MAE_km=float(np.mean(abs(bp-d['b_true']))),endpoint_RMSE_km=float(np.sqrt(np.mean((bp-d['b_true'])**2))),endpoint_active_MAE_km=float(np.mean(abs(bp[active]-d['b_true'][active]))),outlet_H_RMSE_m=float(np.sqrt(np.mean((hp-d['z'][:,0])**2))),SB_same_predicted_location_RMSE_m3=float(np.sqrt(np.mean((sb-hecsb)**2))),SU_same_predicted_location_RMSE_m3=float(np.sqrt(np.mean((su-hecsu)**2))),total_storage_anomaly_RMSE_m3=float(np.sqrt(np.mean(((sb+su-sb[0]-su[0])-(cum[:,-1]-cum[0,-1]))**2))),SB_peak_time_error_h=float(np.nanargmax(sb)-np.argmax(d['sb_true'])),max_iterations=float(np.max(pred['iterations'])),mean_iterations=float(np.mean(pred['iterations'][1:])),mean_absolute_sweep_rate=float(np.nanmean(abs(pred['sweep_rate']))),oracle_boundary_threshold_RMSE_m=float(np.sqrt(np.mean(boundary_error**2))) if name=='Joint_oracle' else np.nan))
  if subset=='validation':
   for name in ['RWEC0','NLM','Serial8','Joint_oracle','Joint','Feedback_off']:
    fit=fits[name]
    for factor in [.5,2.]:
     pred,qual=run(d,name,fit,fit['dt_h']*factor)
     for j,x in enumerate(TARGETS):sens.append(dict(event=e,model=name,dt_factor=factor,dt_h=fit['dt_h']*factor,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),complete=qual['failed_step']==0,failed_time_h=qual['failed_time_h'],**metrics(truth[j],pred['hourly_mean'][j+4])))
 frame('fixed_section_event_metrics.csv',rows);frame('dynamic_location_event_metrics.csv',dynamic);frame('joint_state_event_metrics.csv',state);frame('solver_quality.csv',quality);frame('fixed_parameter_dt_sensitivity.csv',sens)
 df=pd.DataFrame(rows);df.groupby(['subset','model','x_from_outlet_km','distance_from_inlet_km'],sort=False).agg(events=('event','count'),completed=('complete','sum'),NSE=('NSE','mean'),RMSE=('RMSE','mean'),KGE=('KGE','mean'),peak_pct=('peak_pct','mean'),peak_abs_pct=('peak_pct',lambda v:v.abs().mean()),peak_time_h=('peak_time_h','mean'),peak_abs_time_h=('peak_time_h',lambda v:v.abs().mean())).reset_index().to_csv(P/'tables/fixed_section_means.csv',index=False,encoding='utf-8-sig')
 pivot=df.pivot(index=['event','subset','x_from_outlet_km','distance_from_inlet_km'],columns='model',values='NSE').reset_index()
 for model in ['RWEC0','NLM','Serial8','Feedback_off']:pivot['Joint_minus_'+model]=pivot['Joint']-pivot[model]
 frame('paired_NSE_changes.csv',pivot)
 print(df[df.subset=='validation'].groupby('model').agg(NSE=('NSE','mean'),RMSE=('RMSE','mean'),KGE=('KGE','mean'),completed=('complete','sum')).to_string());print(pd.DataFrame(state).query("subset=='validation'").groupby('model').mean(numeric_only=True).to_string())
if __name__=='__main__':main()
