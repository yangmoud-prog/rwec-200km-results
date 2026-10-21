from common import *
from validation import data_A,data_B,weights,fast_base,NODES,COLS
from routing_core import forecast
from scipy.integrate import cumulative_trapezoid
import sys,time

def run(label):
    groups=data_A() if label=='A' else data_B();lookup=loadnp(OUT/'data'/f'{label}_lookup.npz');fits=read(OUT/'fits'/f'{label}_all_final.json')
    core=fits[0]['core'];grid=np.array([0.,.001,.005,.01,.025,.05,.1,.3,1.])
    gains=[];stress=[];profiles=[];source=[];bounds=[]
    losses={}
    for eid,g in groups.items():
        weights(g);train=g['subset']=='train';x=np.arange(0.,200.1,2.)
        if label=='A':
            d=loadnp(P4/'data'/f'A_{eid}.npz');tt=d['time_h'];actual_dz=d['delta_stage_m'];dv=d['delta_volume_m3'];actual_v=np.column_stack([dv[:,-1]-dv[:,int(xc//10)] for xc in TARGETS]);truth_h=d['boundary_increment_m']+5
            n=min(len(c['time']) for c in g['cases']);t=g['time'][:n];mask=t>=0
            yy=np.array([c['truth'][:n][mask] for c in g['cases']]);mu=yy.mean(0)
            for j,xc in enumerate(NODES):
                spread=yy[:,:,j].max(0)-yy[:,:,j].min(0)
                bounds.append(dict(event=eid,subset=g['subset'],xc_km=int(xc),scenarios=len(g['cases']),common_start_h=0,common_end_h=float(t[-1]),
                    conditional_RMSE_floor_m3s=float(np.sqrt(np.mean((yy[:,:,j]-mu[:,j])**2))),maximum_error_floor_m3s=float(spread.max()/2)))
        else:
            d=loadnp(OUT/'data'/f'B_{eid}.npz');tt=d['time'];actual_dz=d.get('delta_stage');dv=d.get('dv');actual_v=None if dv is None else np.column_stack([dv[:,-1]-dv[:,int(xc//2)] for xc in TARGETS]);truth_h=d['z_control'][:,0]
        for mode in [1,2]:
            for alpha in grid:
                r=forecast(g['inflow'],1.,core,float(alpha),mode,lookup);p=r[0][:,COLS]
                score=np.sum(g['W']*(p-g['mu'])**2,axis=0)
                # Constant irreducibility is alpha-independent, so omit only from selection.
                gains.append(dict(event=eid,subset=g['subset'],mode=mode,alpha=float(alpha),loss_mean_without_constant=float(score.mean()),nonconverged_steps=int(r[6].sum()),out_of_lookup_steps=int(r[9])))
                if train:losses.setdefault((mode,float(alpha)),[]).append(score)
                if label=='A' and alpha in [.01,.05]:
                    for c in g['cases']:
                        n=len(c['time']);mask=c['time']>=0
                        for j,xc in enumerate(TARGETS):
                            stress.append(dict(event=eid,subset=g['subset'],case=c['case'],condition=c['condition'],model=f'RWEC-{"L" if mode==1 else "V"}_alpha{alpha:g}',xc_km=int(xc),alpha=float(alpha),nonconverged_steps=int(r[6].sum()),**metrics(c['truth'][mask,j],r[0][:n][mask,(200-xc)//20],c['time'][mask])))
                    if alpha==.01:
                        np.savez_compressed(OUT/'predictions'/f'A_{eid}_nonzero_mode{mode}.npz',time=g['time'],q=r[0],h=r[1],lb_km=r[2]/1000,delta_stage=r[3])
            # Geometric response errors of the actually calibrated autonomous model.
            selected=next(f for f in fits if f['mode']==mode and f['core_fixed'])
            r=forecast(g['inflow'],1.,selected['core'],selected['alpha'],mode,lookup)
            n=min(len(tt),len(r[0]));mask=tt[:n]>=0
            if actual_dz is not None:
                residual=r[3][:n][mask]-actual_dz[:n][mask]
                vh=cumulative_trapezoid(r[4],x*1000,axis=1,initial=0)
                for j,xc in enumerate(TARGETS):
                    predv=vh[:n,-1]-vh[:n,int(xc//2)];volerr=predv[mask]-actual_v[:n,j][mask]
                    profiles.append(dict(event=eid,subset=g['subset'],model=selected['model'],xc_km=int(xc),
                        delta_profile_RMSE_m=float(np.sqrt(np.mean(residual[:,int(xc//2):]**2))),geometry_volume_RMSE_m3=float(np.sqrt(np.mean(volerr**2))),
                        geometry_volume_bias_m3=float(volerr.mean()),H_RMSE_m=float(np.sqrt(np.mean((r[1][:n][mask]-truth_h[:n][mask])**2)))))
    stability=[]
    for mode in [1,2]:
        array=np.array([np.mean(losses[mode,float(a)],axis=0) for a in grid])
        for j,xc in enumerate(NODES):
            k=int(np.argmin(array[:,j]));stability.append(dict(mode=mode,xc_km=int(xc),best_training_grid_alpha=float(grid[k]),gain_vs_alpha0_loss=float(array[0,j]-array[k,j]),grid_min=float(grid[0]),grid_max=float(grid[-1]),lower_bound_selected=bool(k==0),upper_bound_selected=bool(k==len(grid)-1)))
    frame(f'{label}_alpha_grid_by_event.csv',gains);frame(f'{label}_position_alpha_stability.csv',stability);frame(f'{label}_profile_volume_errors.csv',profiles)
    if stress:frame('A_nonzero_gain_stress_metrics.csv',stress)
    if bounds:frame('A_information_bounds_flood_window.csv',bounds)
    print(label,'NONZERO AND PARAMETER DIAGNOSTICS COMPLETE',flush=True)
if __name__=='__main__':run(sys.argv[1])
