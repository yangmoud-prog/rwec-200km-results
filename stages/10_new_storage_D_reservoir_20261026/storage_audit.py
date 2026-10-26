"""Record already-frozen model storage; no new fitting or changed discharge.
Effective/bookkeeping storage remains distinct from HEC geometric volume.
"""
from common10 import *
from formal_baselines import prepare
from formulas import simulate
from scipy.integrate import cumulative_trapezoid
def main():
 cfg=read(D/'config/D_frozen.json');rows=[];balances=[]
 for domain in ['DUF','DUR','DCF','DCR']:
  ds=prepare(domain,cfg['train']+cfg['validation'],cfg)
  for event,d in zip(cfg['train']+cfg['validation'],ds):
   subset='train' if event in cfg['train'] else 'validation';x=d['x'];truecum=cumulative_trapezoid(d['area'],x*1000,axis=1,initial=0);true_storage=np.array([truecum[:,-1]-truecum[:,np.argmin(abs(x-target))] for target in TARGETS]);true_anomaly=true_storage-true_storage[:,[0]]
   for name in ['original','new','general','NLM']:
    fit=read(P/f'PARAMETERS_{domain}_{name}.json')['reaches'];old=load(P/'predictions'/f'{domain}_{event}_{name}.npz');fine=old['time'];up=np.interp(fine,d['time'],d['inflow']);volumes=[];effective=[];projection=[];native_error=[];qrows=[up]
    for reach,f in enumerate(fit):
     dt=f['dt_h'];grid=np.arange(0,cfg['end_h']+dt/2,dt);inp=np.interp(grid,fine,up);q,s,vol,ok,step,bal=simulate(np.array(f['theta']),f['model_id'],inp,dt,np.array([d['initial_volumes'][reach]]),1);assert ok
     inf=np.interp(fine,grid,inp);out=np.interp(fine,grid,q[1]);vv=d['initial_volumes'][reach]+cumulative_trapezoid(inf-out,fine*3600,initial=0)
     # Compare native timestamps with the actual solver state, allowing only its initial-volume offset.
     nv=vol[0]-vol[0,0]+d['initial_volumes'][reach];native_error.append(float(abs(np.interp(grid,fine,vv)-nv).max()))
     projection.append(cumulative_trapezoid(up-inf,fine*3600,initial=0));volumes.append(vv);effective.append(np.interp(fine,grid,s[0]));up=out;qrows.append(out)
    qq=np.array(qrows);assert np.array_equal(qq,old['q'])
    volumes=np.array(volumes);effective=np.array(effective);projection=np.array(projection);native_endpoints=np.cumsum(volumes,axis=0);projected_endpoints=np.cumsum(volumes+projection,axis=0);net=cumulative_trapezoid(qq[0]-qq[-1],fine*3600,initial=0);anomaly=native_endpoints[-1]-native_endpoints[-1,0];rawerr=anomaly-net;projected=projected_endpoints[-1]-projected_endpoints[-1,0]-net;floodvol=float(np.trapz(np.maximum(qq[0]-1000,0),fine*3600))
    balances.append(dict(domain=domain,event=event,subset=subset,model=name,native_grid_state_match_max_m3=max(native_error),cascade_unadjusted_volume_residual_max_m3=float(abs(rawerr).max()),cascade_residual_fraction_flood_volume=float(abs(rawerr).max()/max(floodvol,1)),total_resampling_projection_max_m3=float(abs(projection.sum(axis=0)).max()),with_explicit_projection_residual_max_m3=float(abs(projected).max()),projection_is_numerical_not_geometric_storage=True))
    # Primary diagnostic does NOT hide resampling error behind the projection state.
    hourly_vol=native_endpoints[3:10,::4];pred_anomaly=hourly_vol-hourly_vol[:,[0]];pred_rate=np.diff(hourly_vol,axis=1)/3600;true_rate=np.diff(true_storage,axis=1)/3600
    for j,target in enumerate(TARGETS):
     err=pred_anomaly[j]-true_anomaly[j];de=pred_rate[j]-true_rate[j]
     rows.append(dict(domain=domain,event=event,subset=subset,model=name,x_from_outlet_km=int(target),distance_from_inlet_km=200-int(target),storage_anomaly_RMSE_m3=float(np.sqrt(np.mean(err**2))),storage_anomaly_NRMSE_range=float(np.sqrt(np.mean(err**2))/max(np.ptp(true_anomaly[j]),1)),storage_rate_RMSE_m3s=float(np.sqrt(np.mean(de**2))),storage_rate_mean_error_m3s=float(de.mean()),interpretation='effective bookkeeping anomaly diagnostic; not independently validated geometric state'))
    np.savez_compressed(P/'predictions'/f'storage_{domain}_{event}_{name}.npz',time=fine,reach_volume_bookkeeping=volumes,effective_state=effective,inter_reach_resampling_volume=projection,HEC_time=d['time'],HEC_geometric_storage=true_storage,target_x_km=TARGETS)
 frame('D_storage_state_event_metrics.csv',rows);frame('D_cascade_mass_balance.csv',balances)
 df=pd.DataFrame(rows);sm=df.groupby(['domain','subset','model','x_from_outlet_km','distance_from_inlet_km'])[['storage_anomaly_RMSE_m3','storage_anomaly_NRMSE_range','storage_rate_RMSE_m3s','storage_rate_mean_error_m3s']].agg(['mean','median','min','max']);sm.columns=['_'.join(c) for c in sm.columns];frame('D_storage_state_statistics.csv',sm.reset_index())
 save(P/'STORAGE_AUDIT.json',dict(predictions_reproduced_identically=len(balances),model_changed=False,parameters_changed=False,maximum_cascade_residual_fraction=max(r['cascade_residual_fraction_flood_volume'] for r in balances),maximum_unadjusted_volume_residual_m3=max(r['cascade_unadjusted_volume_residual_max_m3'] for r in balances),maximum_projection_corrected_residual_m3=max(r['with_explicit_projection_residual_max_m3'] for r in balances),interpretation='Native reach mass balance does not remove inter-reach time-resampling discrepancies. Their numerical projection volumes are recorded explicitly; they are not physical backwater storage and do not prove a geometric closure.'))
 print(read(P/'STORAGE_AUDIT.json'))
if __name__=='__main__':main()
