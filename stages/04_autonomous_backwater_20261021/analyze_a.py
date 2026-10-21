"""Existing paired data: mechanism and information limits, no HEC execution."""
import math
import h5py,numpy as np,pandas as pd
from scipy.integrate import cumulative_trapezoid
from common import *

def run():
    dirs();plan=read(A/'config/formal_plan.json');split=read(A/'analysis/03_internal_discovery_20261021/provenance/event_split.json')
    train=set(split['train']);validation=set(split['validation']);assert len(train)==49 and len(validation)==21
    main={}
    for c in plan['formal_cases']:
        if c['split']=='calibration' and c['role']=='main':main.setdefault(c['event_id'],{})[c['kind']]=c
    assert len(main)==70
    bounds=[];process=[];events=[];source=[];derived=[];cache={}
    for eid,pair in sorted(main.items()):
        raw={}
        for kind,c in pair.items():
            p=A/'results'/c['case']/'reach_timeseries.npz';raw[kind]=loadnp(p)
            source.append(dict(path=str(p.relative_to(A)),sha256=sha(p),mtime_ns=p.stat().st_mtime_ns))
            cache[c['case']]=raw[kind]
        assert raw['var']['flow_m3s'].shape==raw['ref']['flow_m3s'].shape
        c=pair['var'];p=A/'analysis/01_backwater_storage_20261020/derived_cases'/f'{c["case"]}.npz';z=loadnp(p)
        source.append(dict(path=str(p.relative_to(A)),sha256=sha(p),mtime_ns=p.stat().st_mtime_ns))
        t=z['time_h'];x=z['full_x_km'];dz=z['delta_stage_full_m'].astype(float);vol=z['delta_volume_m3']
        # Reference geometry is a diagnostic observation, never a forecast input.
        path=A/'models'/pair['ref']['case']/'Backwater200.p01.hdf'
        with h5py.File(path,'r') as f:
            refz=f['Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Cross Sections/Water Surface'][:].astype(float)[:,::-1]
        yr=refz-x[None,:]*1000*.00005;da=(500+6*yr)*dz+3*dz**2
        vf=cumulative_trapezoid(da,x*1000,axis=1,initial=0)
        assert np.max(abs(vf[:,::5]-vol))<.1
        ext=[extent(d,x) for d in dz];lb=np.array([v[0] for v in ext]);censored=np.array([v[1] for v in ext]);segments=np.array([v[2] for v in ext])
        h=z['stage_increment_boundary_m'];hpos=np.maximum(h,0)
        subset='train' if eid in train else 'validation'
        packet=dict(event_id=eid,subset=subset,time=t,lb=lb,censored=censored,vol=vol,h=h,phase=np.array([phase(tt,c['event']) for tt in t]))
        derived.append(packet)
        for xc in TARGETS:
            j=int(xc//10);leff=np.maximum(0,np.minimum(lb,200)-xc)*1000
            dv=vol[:,-1]-vol[:,j];rect=500*hpos*leff
            outside=np.array([vf[k,-1]-np.interp(max(xc,lb[k]),x,vf[k]) for k in range(len(t))])
            for k,tt in enumerate(t):
                process.append(dict(event_id=eid,subset=subset,xc_km=int(xc),Lr_km=int(200-xc),time_h=float(tt),phase=packet['phase'][k],
                    lb_reference_km=float(lb[k]),right_censored=bool(censored[k]),positive_islands=int(segments[k]),
                    effective_length_m=float(leff[k]),geometric_delta_volume_m3=float(dv[k]),
                    boundary_stage_increment_m=float(h[k]),rectangle_proxy_m3=float(rect[k]),
                    signed_volume_beyond_threshold_m3=float(outside[k])))
            events.append(dict(event_id=eid,subset=subset,xc_km=int(xc),max_lb_km=float(lb.max()),max_leff_km=float(leff.max()/1000),
                peak_delta_volume_m3=float(dv.max()),min_delta_volume_m3=float(dv.min()),
                times_zero_leff_nonzero_volume=int(((leff==0)&(abs(dv)>1e6)).sum()),
                right_censored_samples=int(censored.sum()),negative_target_saved_samples=int((raw['var']['flow_m3s'][:,j]<0).sum()),
                max_nonconnected_islands=int(segments.max())))
        information_group(eid,subset,'main_pair',[pair['ref']['case'],c['case']],cache,bounds)
        np.savez_compressed(OUT/'data'/f'A_{eid}.npz',time_h=t,lb_reference_km=lb,delta_stage_m=dz,full_x_km=x,
            delta_volume_m3=vol,boundary_increment_m=h,flow_var_m3s=raw['var']['flow_m3s'],flow_ref_m3s=raw['ref']['flow_m3s'])
    for eid,g in sorted(plan['groups'].items()):
        ids=[main[eid]['ref']['case'],main[eid]['var']['case']]+g['addon']
        assert len(ids)==10
        for case in g['addon']:
            p=A/'results'/case/'reach_timeseries.npz';cache[case]=loadnp(p)
            source.append(dict(path=str(p.relative_to(A)),sha256=sha(p),mtime_ns=p.stat().st_mtime_ns))
        information_group(eid,'train' if eid in train else 'validation','ten_boundary_scenarios',ids,cache,bounds)
    frame('A_information_bounds.csv',bounds);frame('A_length_volume_process.csv',process);frame('A_event_geometry.csv',events)
    P=pd.DataFrame(process);fits=[];metrics=[]
    # Diagnostic geometric fits: do not equate these gains with RWEC potential storage.
    for xc in TARGETS:
        sub=P[P.xc_km==xc];tr=sub[sub.subset=='train'];weights=1/tr.groupby('event_id').time_h.transform('size').to_numpy()
        for name,predictor in [('length_only','effective_length_m'),('boundary_height_times_length','rectangle_proxy_m3')]:
            xx=tr[predictor].to_numpy();yy=tr.geometric_delta_volume_m3.to_numpy()
            gain=float(np.sum(weights*xx*yy)/max(1e-30,np.sum(weights*xx*xx)))
            fits.append(dict(xc_km=int(xc),diagnostic=name,gain=gain,training_events=49,
                gain_units='m2' if name=='length_only' else 'dimensionless',runtime_forecast_claimed=False))
            for eid,g in sub.groupby('event_id'):
                y=g.geometric_delta_volume_m3.to_numpy();q=g[predictor].to_numpy()*gain;de=np.sum((y-y.mean())**2)
                metrics.append(dict(event_id=eid,subset=g.subset.iloc[0],xc_km=int(xc),diagnostic=name,
                    RMSE_m3=float(np.sqrt(np.mean((q-y)**2))),NSE=float(1-np.sum((q-y)**2)/de),
                    bias_m3=float(np.mean(q-y))))
    frame('A_geometric_diagnostic_parameters.csv',fits);frame('A_geometric_diagnostic_metrics.csv',metrics)
    for name,df,keys,cols in [('A_information_summary.csv',pd.DataFrame(bounds),['group','subset','xc_km'],['optimal_pooled_RMSE_floor_m3s','mean_NSE_ceiling','minimax_saved_point_error_floor_m3s']),
        ('A_geometry_summary.csv',pd.DataFrame(events),['subset','xc_km'],['max_lb_km','peak_delta_volume_m3','times_zero_leff_nonzero_volume']),
        ('A_geometric_diagnostic_summary.csv',pd.DataFrame(metrics),['diagnostic','subset','xc_km'],['RMSE_m3','NSE'])]:
        s=df.groupby(keys)[cols].agg(['mean','median','min','max']);s.columns=['_'.join(c) for c in s.columns]
        s.reset_index().to_csv(OUT/'tables'/name,index=False,encoding='utf-8-sig')
    save(OUT/'provenance/A_input_manifest.json',dict(inputs=source,sealed_test_results_opened=0,
        reference_native_HDF_read_only=70,source_selection='70 development pairs plus32 existing controls only',
        same_prescribed_inflow_and_initial_state_grouped=True,
        raw_inflow_solver_roundoff_not_used_as_a_scenario_identifier=True))
    save(OUT/'provenance/event_split.json',split)
    save(OUT/'checkpoint.json',dict(stage='A_mechanism_and_identifiability_complete',B_pilot_status='preparing',test_read=0))
    print('A完成：70配对、4组十边界场景的预测信息下界；6个内部位置',flush=True)

def information_group(eid,subset,label,ids,cache,rows):
    n=min(len(cache[k]['time_since_flood_h']) for k in ids)
    t=cache[ids[0]]['time_since_flood_h'][:n]
    for key in ids:assert np.allclose(cache[key]['time_since_flood_h'][:n],t,atol=.02/3600,rtol=0)
    for xc in np.r_[0,TARGETS]:
        Y=np.array([cache[k]['flow_m3s'][:n,int(xc//10)] for k in ids])
        mu=Y.mean(axis=0);mse=float(np.mean((Y-mu)**2));den=np.sum((Y-Y.mean(axis=1,keepdims=True))**2,axis=1)
        w=1/den;opt=np.sum(Y*w[:,None],axis=0)/w.sum();loss=np.mean(np.sum((Y-opt)**2,axis=1)/den)
        spread=Y.max(axis=0)-Y.min(axis=0)
        rows.append(dict(event_id=eid,subset=subset,group=label,xc_km=int(xc),scenarios=len(ids),
            common_time_start_h=float(t[0]),common_time_end_h=float(t[-1]),
            irreducible_MSE_m6s2=mse,optimal_pooled_RMSE_floor_m3s=math.sqrt(mse),
            minimax_saved_point_error_floor_m3s=float(spread.max()/2),mean_NSE_ceiling=float(1-loss),
            max_scenario_spread_m3s=float(spread.max()),
            negative_truth_samples=int((Y<0).sum()),
            condition='identical prescribed upstream hydrograph, geometry and physical initial state; future forcing hidden',
            ensemble_weights='equal empirical scenario weights; not asserted real-world probabilities'))

if __name__=='__main__':run()
