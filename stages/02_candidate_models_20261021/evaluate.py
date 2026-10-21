"""Development metrics, fixed-parameter controls and conditional stability."""
from pathlib import Path
import math,time
import numpy as np
import models
from campaign import OUT,ROOT,LENGTHS,SEED,read,save,load_np,csvwrite,optimize,stamp

def metrics(q,obs,t):
    q=np.asarray(q);obs=np.asarray(obs)
    den=float(np.sum((obs-obs.mean())**2));v=float(np.std(obs));mean=float(obs.mean())
    nse=1-float(np.sum((q-obs)**2))/den if den>1e-8*len(obs) else float('nan')
    r=float(np.corrcoef(q,obs)[0,1]) if np.std(q)>1e-10 and v>1e-10 else float('nan')
    kge=1-math.sqrt((r-1)**2+(float(np.std(q))/v-1)**2+(float(q.mean())/mean-1)**2) if v>1e-10 and abs(mean)>1e-10 and np.isfinite(r) else float('nan')
    peak=float(obs.max());peakerr=(float(q.max())-peak)/peak*100 if abs(peak)>1e-12 else float('nan')
    return dict(NSE=nse,KGE=kge,RMSE=float(np.sqrt(np.mean((q-obs)**2))),
        peak_error_pct=peakerr,peak_time_error_h=float(t[q.argmax()]-t[obs.argmax()]),
        bias_m3s=float(np.mean(q-obs)),MAE=float(np.mean(abs(q-obs))))

def aggregate(rows,keys,metric_keys):
    groups={}
    for row in rows:groups.setdefault(tuple(row[k] for k in keys),[]).append(row)
    result=[]
    for key,rs in groups.items():
        for met in metric_keys:
            val=np.array([r[met] for r in rs],dtype=float);val=val[np.isfinite(val)]
            result.append(dict(**dict(zip(keys,key)),metric=met,n_records=len(rs),n_finite=len(val),
                mean=float(val.mean()) if len(val) else float('nan'),median=float(np.median(val)) if len(val) else float('nan'),
                min=float(val.min()) if len(val) else float('nan'),max=float(val.max()) if len(val) else float('nan')))
    return result

def phase_masks(e,t):
    tb=e['tb_h'];top=tb+e['Tr_h'];fall=top+e['Th_h'];end=e['Hend_h']
    return {'before_stage_rise':t<tb,'stage_rise':(t>=tb)&(t<top),
        'stage_hold':(t>=top)&(t<fall),'stage_drawdown':(t>=fall)&(t<end),
        'after_stage_return':t>=end,'flood_recession':(t>=e['Tp_h'])&(t<=e['D_h'])}

def main():
    assert read(OUT/'checkpoint.json')['status']=='primary_calibration_complete'
    fits=[read(OUT/'fits'/f'model{m}_final.json') for m in range(3)]
    pars=np.array([[models.decode(x) for x in f['x']] for f in fits])
    shared=[(0.,0.)]+[models.decode_shared(fits[m]['shared_x'],m) for m in [1,2]]
    metadata=read(OUT/'data/metadata.json')
    rows=[];phases=[];qadefects=[];param_rows=[]
    for m in range(3):
        for li,L in enumerate(LENGTHS):
            k,n,tau=pars[m,li];alpha,theta=shared[m]
            x=fits[m]['x'][li]
            param_rows.append(dict(model=f'RWEC-{m}',L_km=int(L/1000),k=k,n=n,tau_h=tau,
                alpha=alpha,theta=theta,ell_m=float('inf') if theta==0 else 200000/theta,
                near_k_bound=bool(min(x[0],1-x[0])<.01),near_n_bound=bool(min(x[1],1-x[1])<.01),
                near_tau_bound=bool(min(x[2],1-x[2])<.01),
                k_at_bound=bool(min(x[0],1-x[0])<1e-5),n_at_bound=bool(min(x[1],1-x[1])<1e-5),
                tau_at_bound=bool(min(x[2],1-x[2])<1e-5)))
    csvwrite(OUT/'tables/parameters.csv',param_rows)
    for tag in ['development','controls']:
        d=load_np(OUT/'data'/f'{tag}.npz');n=len(d['ns']);maxN=d['qin'].shape[-1]
        pred=np.zeros((3,20,n,maxN));states=np.zeros_like(pred)
        for m in range(3):
            for li,L in enumerate(LENGTHS):
                k,nc,tau=pars[m,li];alpha,theta=shared[m]
                for j,N in enumerate(d['ns']):
                    t=np.arange(N)-48.;qin=d['qin'][li,j,:N];h=d['h'][j,:N];obs=d['obs'][j,:N]
                    q,s,ib,res=models.route(qin,h,N,k,nc,tau,alpha,theta,float(L))
                    assert np.all(np.isfinite(q)) and np.max(abs(res))<.01
                    pred[m,li,j,:N]=q;pred[m,li,j,N:]=q[-1]
                    states[m,li,j,:N]=s;states[m,li,j,N:]=s[-1]
                    raw=.5*(np.interp(t-tau,t,qin)+np.interp(t-tau-1,t,qin))
                    rec=metadata[tag][j];base=dict(dataset=tag,model=f'RWEC-{m}',L_km=int(L/1000),
                        event_id=rec['event_id'],case=rec['case'])
                    met=metrics(q,obs,t)
                    rows.append(dict(**base,**met,n_times=int(N),negative_target_samples=int((obs<0).sum()),
                        min_target_flow_m3s=float(obs.min()),negative_input_samples=int((qin<0).sum()),
                        clipped_lagged_input_samples=int((raw[1:]<0).sum()),min_simulated_flow_m3s=float(q.min()),
                        nonnegative_model_NSE_upper_bound=1-float(np.sum(np.minimum(obs,0)**2))/float(d['denom'][j]),
                        max_internal_balance_residual_m3=float(abs(res).max())))
                    for name,mask in phase_masks(rec['event'],t).items():
                        if mask.sum()>=2:phases.append(dict(**base,phase=name,n_times=int(mask.sum()),**metrics(q[mask],obs[mask],t[mask])))
        np.savez_compressed(OUT/'predictions'/f'{tag}.npz',flow=pred,storage_million_m3=states,ns=d['ns'])
        print('Evaluated',tag,'all three models.',flush=True)
    metric_keys=['NSE','KGE','RMSE','peak_error_pct','peak_time_error_h','bias_m3s','MAE']
    csvwrite(OUT/'tables/event_length_metrics.csv',rows)
    csvwrite(OUT/'tables/phase_event_metrics.csv',phases)
    csvwrite(OUT/'tables/metric_summary.csv',aggregate(rows,['dataset','model'],metric_keys))
    csvwrite(OUT/'tables/metric_by_length.csv',aggregate(rows,['dataset','model','L_km'],metric_keys))
    csvwrite(OUT/'tables/phase_summary.csv',aggregate(phases,['dataset','model','phase'],metric_keys))
    csvwrite(OUT/'tables/phase_by_length.csv',aggregate(phases,['dataset','model','L_km','phase'],metric_keys))
    base={(r['dataset'],r['case'],r['L_km']):r for r in rows if r['model']=='RWEC-0'}
    improvement=[]
    for r in rows:
        if r['model']=='RWEC-0':continue
        b=base[r['dataset'],r['case'],r['L_km']]
        improvement.append(dict(dataset=r['dataset'],model=r['model'],case=r['case'],event_id=r['event_id'],L_km=r['L_km'],
            delta_NSE=r['NSE']-b['NSE'],delta_KGE=r['KGE']-b['KGE'],
            RMSE_reduction_m3s=b['RMSE']-r['RMSE'],RMSE_reduction_pct=(b['RMSE']-r['RMSE'])/b['RMSE']*100,
            improved_NSE=bool(r['NSE']-b['NSE']>1e-6),worse_NSE=bool(r['NSE']-b['NSE']< -1e-6)))
    csvwrite(OUT/'tables/improvements.csv',improvement)
    csvwrite(OUT/'tables/improvements_by_length.csv',aggregate(improvement,['dataset','model','L_km'],['delta_NSE','delta_KGE','RMSE_reduction_m3s','RMSE_reduction_pct']))
    dev=metadata['development'];preset=['E023','E006','E007']
    candidates=sorted({r['event_id'] for r in improvement if r['dataset']=='development'}-set(preset))
    worst=min(candidates,key=lambda e:np.mean([r['delta_NSE'] for r in improvement if r['dataset']=='development' and r['model']=='RWEC-2' and r['event_id']==e]))
    save(OUT/'provenance/representative_events.json',dict(preset=preset,worst_additional=worst,
        rule='Worst event mean delta NSE of RWEC-2 versus0, among remaining development events; all20lengths equal weight.'))
    controls_analysis(metadata,rows)
    stability(fits,pars,shared)
    step_sensitivity(pars,shared,metadata)
    product_profile(fits,pars,shared)
    summary=[]
    for m in range(3):
        rr=[r for r in rows if r['dataset']=='development' and r['model']==f'RWEC-{m}']
        cc=[r for r in rows if r['dataset']=='controls' and r['model']==f'RWEC-{m}']
        summary.append(dict(model=f'RWEC-{m}',free_parameters=60+m,
            development_mean_NSE=float(np.mean([r['NSE'] for r in rr])),development_median_NSE=float(np.median([r['NSE'] for r in rr])),
            development_mean_KGE=float(np.mean([r['KGE'] for r in rr])),development_mean_RMSE=float(np.mean([r['RMSE'] for r in rr])),
            controls_mean_NSE=float(np.mean([r['NSE'] for r in cc])),controls_mean_RMSE=float(np.mean([r['RMSE'] for r in cc])),
            alpha=shared[m][0],theta=shared[m][1],ell_m=None if shared[m][1]==0 else 200000/shared[m][1],
            charged_reach_nfev=fits[m]['charged_reach_nfev'],near_bounds={v:sum(r[v] for r in param_rows if r['model']==f'RWEC-{m}') for v in ['near_k_bound','near_n_bound','near_tau_bound']},
            improved_of1400=0 if m==0 else sum(r['improved_NSE'] for r in improvement if r['dataset']=='development' and r['model']==f'RWEC-{m}'),
            worse_of1400=0 if m==0 else sum(r['worse_NSE'] for r in improvement if r['dataset']=='development' and r['model']==f'RWEC-{m}')))
    save(OUT/'summary.json',dict(created=stamp(),models=summary,representative_events=preset+[worst],
        test_events_used=0,development_events=70,controlled_runs=32,main_event_length_cases=1400,
        accuracy_scope='In-sample development fits and unfitted boundary controls on four development flood backgrounds; not independent flood validation.',
        all_predictions_finite=True,all_failures_and_reversals_retained=True))
    save(OUT/'checkpoint.json',dict(status='analysis_complete_plotting_pending',updated=stamp(),test_results_read=False))
    print('SUMMARY',summary,flush=True)

def controls_analysis(metadata,rows):
    d=load_np(OUT/'data/controls.npz');p=load_np(OUT/'predictions/controls.npz')['flow']
    result=[]
    pairs={'A':[(0,4),(1,5),(2,6),(3,7)],'B':[(0,2),(1,3),(4,6),(5,7)],'C':[(0,1),(2,3),(4,5),(6,7)]}
    for eid in sorted({r['event_id'] for r in metadata['controls']}):
        idx=[j for j,r in enumerate(metadata['controls']) if r['event_id']==eid]
        assert len(idx)==8
        N=int(d['ns'][idx[0]]);assert all(d['ns'][j]==N for j in idx)
        diffs=np.zeros((3,4,4,20,N))
        for fi,(factor,combinations) in enumerate(pairs.items()):
            for ci,(lo,hi) in enumerate(combinations):
                a,b=idx[lo],idx[hi];truth=d['obs'][b,:N]-d['obs'][a,:N]
                diffs[fi,ci,0]=truth[None,:]
                for m in range(3):
                    dq=p[m,:,b,:N]-p[m,:,a,:N];diffs[fi,ci,m+1]=dq
                    for li,L in enumerate(LENGTHS):
                        met=metrics(dq[li],truth,np.arange(N)-48)
                        result.append(dict(event_id=eid,factor=factor,low=f'C{lo+1:02d}',high=f'C{hi+1:02d}',
                            model=f'RWEC-{m}',L_km=int(L/1000),RMSE=met['RMSE'],NSE=met['NSE'],bias_m3s=met['bias_m3s'],
                            max_abs_response_error_m3s=float(abs(dq[li]-truth).max()),reference_response_range_m3s=float(np.ptp(truth))))
        np.savez_compressed(OUT/'predictions'/f'control_response_{eid}.npz',differences=diffs,time_h=np.arange(N)-48)
    csvwrite(OUT/'tables/control_response_metrics.csv',result)
    csvwrite(OUT/'tables/control_response_summary.csv',aggregate(result,['model','factor'],['RMSE','NSE','bias_m3s','max_abs_response_error_m3s']))
    csvwrite(OUT/'tables/control_response_by_length.csv',aggregate(result,['model','factor','L_km'],['RMSE','NSE','bias_m3s','max_abs_response_error_m3s']))

def stability(fits,pars,shared):
    path=OUT/'fits/conditional_stability.json'
    if path.exists():return
    d=load_np(OUT/'data/development.npz');allidx=np.arange(70,dtype=np.int64)
    groups=[('event_half_even',allidx[::2],range(20)),('event_half_odd',allidx[1::2],range(20)),
        ('length_10_60km',allidx,range(6)),('length_70_130km',allidx,range(6,13)),('length_140_200km',allidx,range(13,20))]
    group_rows=[];bylength=[];log=[]
    for m in [1,2]:
        for gi,(name,ids,lis) in enumerate(groups):
            def fn(x):
                a,theta=models.decode_shared(x,m)
                return np.mean([models.loss(d['qin'][i],d['h'],d['obs'],d['ns'],d['denom'],ids,*pars[m,i],a,theta,float(LENGTHS[i])) for i in lis])
            fit=optimize(fn,fits[m]['shared_x'],160,SEED+50000+m*100+gi)
            a,theta=models.decode_shared(fit['x'],m)
            group_rows.append(dict(model=f'RWEC-{m}',group=name,n_events=len(ids),n_lengths=len(list(lis)),alpha=a,
                theta=theta,ell_m=float('inf') if theta==0 else 200000/theta,mean_NSE=1-fit['loss'],nfev=fit['nfev']))
            log.append(dict(model=m,group=name,fit=fit))
        for i,L in enumerate(LENGTHS):
            def one(x):return models.loss(d['qin'][i],d['h'],d['obs'],d['ns'],d['denom'],allidx,*pars[m,i],float(x[0])*10,shared[m][1],float(L))
            fit=optimize(one,[shared[m][0]/10],80,SEED+51000+m*100+i)
            bylength.append(dict(model=f'RWEC-{m}',L_km=int(L/1000),conditional_alpha=fit['x'][0]*10,
                primary_shared_alpha=shared[m][0],theta_fixed=shared[m][1],mean_NSE=1-fit['loss'],nfev=fit['nfev']))
    csvwrite(OUT/'tables/conditional_shared_parameters.csv',group_rows)
    csvwrite(OUT/'tables/conditional_alpha_by_length.csv',bylength)
    save(path,dict(conditional_only=True,base_parameters_fixed_from_all70=True,not_independent_validation=True,groups=log))
    print('Conditional stability diagnostics complete.',flush=True)

def step_sensitivity(pars,shared,metadata):
    d=load_np(OUT/'data/development.npz');rows=[];curves={}
    for eid in ['E023','E006','E007','E069']:
        j=next(i for i,r in enumerate(metadata['development']) if r['event_id']==eid);N=int(d['ns'][j]);t=np.arange(N)
        for li in [0,9,19]:
            L=float(LENGTHS[li])
            for m in range(3):
                coarse=models.route(d['qin'][li,j,:N],d['h'][j,:N],N,*pars[m,li],*shared[m],L)[0]
                for factor in [1,2,4]:
                    ft=np.arange((N-1)*factor+1)/factor
                    fq=np.interp(ft,t,d['qin'][li,j,:N]);fh=np.interp(ft,t,d['h'][j,:N])
                    fine=models.route(fq,fh,len(ft),*pars[m,li],*shared[m],L,1/factor)[0]
                    hourly=np.r_[fine[0],fine[1:].reshape(N-1,factor).mean(axis=1)]
                    rows.append(dict(event_id=eid,L_km=int(L/1000),model=f'RWEC-{m}',dt_h=1/factor,
                        RMSE_vs_1h_m3s=float(np.sqrt(np.mean((hourly-coarse)**2))),max_abs_vs_1h_m3s=float(abs(hourly-coarse).max()),
                        rel_RMSE_vs_observed_peak=float(np.sqrt(np.mean((hourly-coarse)**2))/d['obs'][j,:N].max())))
                    if eid=='E007' and li==19:
                        curves[f'model{m}_factor{factor}']=hourly
                        curves['time_h']=t-48
    csvwrite(OUT/'tables/timestep_sensitivity.csv',rows)
    np.savez_compressed(OUT/'predictions/timestep_example.npz',**curves)

def product_profile(fits,pars,shared):
    d=load_np(OUT/'data/development.npz');ids=np.arange(70,dtype=np.int64)
    alpha,theta=shared[2];product=alpha*200000/theta;rows=[]
    for th in [20,40,60,80,100,theta,140,160,180,200]:
        ell=200000/th;a=product/ell
        value=float(np.mean([models.loss(d['qin'][i],d['h'],d['obs'],d['ns'],d['denom'],ids,
            *pars[2,i],a,th,float(LENGTHS[i])) for i in range(20)]))
        rows.append(dict(theta=th,ell_m=ell,alpha=a,alpha_ell_m=product,mean_NSE=1-value,
            delta_NSE_vs_primary=1-value-(1-fits[2]['mean_loss'])))
    csvwrite(OUT/'tables/saturation_product_profile.csv',rows)

if __name__=='__main__':main()
