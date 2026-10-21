from common import *
from prototype import forecast,oracle
import sys,time,h5py
from scipy.optimize import differential_evolution,minimize,minimize_scalar
from scipy.integrate import cumulative_trapezoid
sys.path.insert(0,str(B/'scripts'))
import hec200 as hec

MODES={'RWEC-0':0,'RWEC-L':1,'RWEC-V':2,'RWEC-H':3,'RWEC-1（真实未来水位）':4}
IDS=['E023','E006','E007','E069']
TRAIN=['E023','E007','E069']

def metrics(o,p,t):
    e=p-o;rmse=np.sqrt(np.mean(e*e));sst=np.sum((o-o.mean())**2)
    r=np.corrcoef(o,p)[0,1]
    return dict(NSE=float(1-np.sum(e*e)/sst),RMSE_m3s=float(rmse),KGE=float(1-np.sqrt((r-1)**2+(p.std()/o.std()-1)**2+(p.mean()/o.mean()-1)**2)),
        peak_error_pct=float(100*(p.max()/o.max()-1)),peak_time_error_h=float(t[np.argmax(p)]-t[np.argmax(o)]))

def prepare():
    rows=[]
    for eid in IDS:
        pairs={}
        for kind in ['control','ref']:
            folder=B/'models'/f'{eid}_{kind}_2km_60s'
            quality=read(B/'results'/folder.name/'quality.json')
            assert quality['basic_quality_pass'] and quality['initialization_pass'] and quality['tail_recovered']
            rows.append(dict(event=eid,kind=kind,**{k:v for k,v in quality.items() if isinstance(v,(float,int,bool,str))}))
            with h5py.File(folder/'Backwater200.p01.hdf','r') as f:
                base=f[hec.TS];xs=base['Cross Sections'];tt=base['Time'][:].astype(float)
                t=np.round((tt-tt[0])*24-48,6)
                pairs[kind]=dict(t=t,q=xs['Flow'][:].astype(float)[:,::-1],z=xs['Water Surface'][:].astype(float)[:,::-1])
        c,r=pairs['control'],pairs['ref'];assert np.array_equal(c['t'],r['t'])
        x=np.arange(0.,200.1,2);bed=x*.1;d=c['z']-r['z']
        y=r['z']-bed;da=(500+6*y)*d+3*d*d
        dv=cumulative_trapezoid(da,x*1000,axis=1,initial=0)
        lengths=np.array([extent(v,x)[0] for v in d])
        e=read(B/'config/pilot_frozen.json')['events'][eid]
        inflow=hec.boundary(e,c['t'],'ref')[0]
        np.savez_compressed(OUT/'data'/f'B_{eid}.npz',time=c['t'],x=x,q_control=c['q'],q_ref=r['q'],z_control=c['z'],z_ref=r['z'],lb=lengths,dv=dv,inflow=inflow)
    frame('B_HEC_quality.csv',rows)

def decode(z):return [float(np.exp(z[0])),float(np.exp(z[1])),float(z[2])]

def main():
    dirs();prepare();lookup=loadnp(OUT/'data/hydraulic_lookup.npz')
    data={eid:loadnp(OUT/'data'/f'B_{eid}.npz') for eid in IDS}
    nodes=np.array([20,40,60,80,100,120,0]);cols=(200-nodes)//20
    masks={eid:data[eid]['time']>=0 for eid in IDS}
    def objective(z):
        params=decode(z);scores=[]
        for eid in TRAIN:
            d=data[eid];m=masks[eid];pred=forecast(d['inflow'],1.,params,0.,0,lookup)[0]
            obs=d['q_ref'][m][:,nodes//2];p=pred[m][:,cols]
            scores.append(np.mean(np.sum((p-obs)**2,axis=0)/np.sum((obs-obs.mean(axis=0))**2,axis=0)))
        return float(np.mean(scores))
    pfile=OUT/'fits/core_parameters.json'
    if pfile.exists():params=read(pfile)['params']
    else:
        start=time.perf_counter()
        bounds=[(np.log(.0001),np.log(20)),(np.log(.15),np.log(150)),(0.,8.)]
        opt=differential_evolution(objective,bounds,seed=20261021,popsize=8,maxiter=35,tol=1.e-5,polish=False)
        opt2=minimize(objective,opt.x,method='Nelder-Mead',bounds=bounds,options=dict(maxiter=300,xatol=1.e-6,fatol=1.e-8))
        z=opt2.x if opt2.fun<opt.fun else opt.x;params=decode(z)
        save(pfile,dict(params=params,names=['k','n','tau_per20km_h'],objective=objective(z),train=TRAIN,validation=['E006'],
             fit_target='fixed5 reference flows at six requested sections plus outlet, equal event and section weights',elapsed_s=time.perf_counter()-start,bounds=bounds))
    print('Core',params,flush=True)
    # Separate one-gain calibration for each correction, with core fixed.
    gains={};trials=[]
    for name,mode in MODES.items():
        if mode==0:gains[name]=0.;continue
        def loss(alpha):
            vals=[];fail=0
            for eid in TRAIN:
                d=data[eid];m=masks[eid]
                rr=(oracle(d['inflow'],1.,params,alpha,d['z_control'][:,0],lookup) if mode==4 else forecast(d['inflow'],1.,params,alpha,mode,lookup))
                o=d['q_control'][m][:,nodes//2];p=rr[0][m][:,cols]
                vals.append(np.mean(np.sum((o-p)**2,axis=0)/np.sum((o-o.mean(axis=0))**2,axis=0)))
                fail+=int(rr[6].sum())+rr[9]
            value=float(np.mean(vals))+10*fail
            trials.append(dict(model=name,alpha=float(alpha),loss=value,nonconvergence_or_domain_count=fail))
            return value
        path=OUT/'fits'/f'gain_mode{mode}.json'
        if path.exists():gains[name]=read(path)['alpha'];continue
        start=time.perf_counter();opt=minimize_scalar(loss,bounds=(0.,1.),method='bounded',options=dict(xatol=2.e-4,maxiter=22))
        zero=loss(0.);alpha=0. if zero<=opt.fun else float(opt.x)
        gains[name]=alpha
        save(path,dict(model=name,alpha=alpha,loss=min(zero,float(opt.fun)),core_fixed=True,train=TRAIN,validation=['E006'],elapsed_s=time.perf_counter()-start))
        print(name,'alpha',alpha,flush=True)
    frame('B_gain_calibration_trials.csv',trials)
    save(OUT/'fits/all_parameters.json',dict(core=params,gains=gains,training_background_events=TRAIN,validation_background_events=['E006'],test_opened=0))
    rows=[];lengthrows=[];checkrows=[];sensitivity=[];process=[]
    for eid,d in data.items():
        m=masks[eid];subset='validation' if eid=='E006' else 'train';t=d['time']
        for name,mode in MODES.items():
            start=time.perf_counter()
            rr=(oracle(d['inflow'],1.,params,gains[name],d['z_control'][:,0],lookup) if mode==4 else forecast(d['inflow'],1.,params,gains[name],mode,lookup))
            elapsed=time.perf_counter()-start
            q,h,lb,dz,da,its,failed,mass,fullmass,domain,fixedres=rr
            np.savez_compressed(OUT/'predictions'/f'{eid}_mode{mode}.npz',time=t,q=q,h=h,lb_km=lb/1000,delta_stage=dz,delta_area=da,iterations=its,failed=failed)
            for xc in TARGETS:
                node=(200-xc)//20;o=d['q_control'][m,xc//2];p=q[m,node]
                rows.append(dict(event=eid,subset=subset,model=name,xc_km=int(xc),Lr_km=int(200-xc),**metrics(o,p,t[m])))
                difference=d['q_control'][:,xc//2]-d['q_ref'][:,xc//2]
                process.append(dict(event=eid,subset=subset,xc_km=int(xc),model=name,
                    suppression_min_m3s=float(difference[m].min()),recovery_max_m3s=float(difference[m].max())))
            if mode in (1,2,3):
                err=lb[m]/1000-d['lb'][m]
                lengthrows.append(dict(event=eid,subset=subset,model=name,RMSE_km=float(np.sqrt(np.mean(err*err))),MAE_km=float(np.mean(abs(err))),bias_km=float(err.mean()),
                    reference_peak_km=float(d['lb'][m].max()),predicted_peak_km=float(lb[m].max()/1000)))
            checkrows.append(dict(event=eid,subset=subset,model=name,elapsed_s=elapsed,nonconverged_steps=int(failed.sum()),max_iterations=int(its.max()),
                 mass_residual_max_m3=float(abs(mass).max()),full_balance_max_m3=float(abs(fullmass).max()),lookup_out_of_range_steps=int(domain),max_accepted_fixed_point_residual_m3s=float(fixedres),minimum_predicted_flow_m3s=float(q.min())))
            if mode in (0,1,2,3):
                # Prefix invariance: altering future inflow cannot change past outputs.
                cut=np.flatnonzero(t>=24)[0];altered=d['inflow'].copy();altered[cut+1:]=1000
                other=forecast(altered,1.,params,gains[name],mode,lookup)
                err=float(abs(q[:cut+1]-other[0][:cut+1]).max())
                assert err==0.
                checkrows[-1]['future_input_prefix_difference_m3s']=err
                for dt in [.5,.25]:
                    fine=np.arange(t[0],t[-1]+dt/2,dt);e=read(B/'config/pilot_frozen.json')['events'][eid]
                    inf=hec.boundary(e,fine,'ref')[0];fr=forecast(inf,dt,params,gains[name],mode,lookup)
                    pf=fr[0][::int(1/dt)]
                    for xc in TARGETS:
                        j=(200-xc)//20
                        sensitivity.append(dict(event=eid,subset=subset,model=name,dt_h=dt,xc_km=int(xc),
                            difference_from_1h_RMSE_m3s=float(np.sqrt(np.mean((pf[m,j]-q[m,j])**2))),
                            NSE_vs_HEC=metrics(d['q_control'][m,xc//2],pf[m,j],t[m])['NSE'],nonconverged_steps=int(fr[6].sum())))
    frame('B_event_metrics.csv',rows);frame('B_length_metrics.csv',lengthrows);frame('B_numerical_causality_checks.csv',checkrows)
    frame('B_timestep_sensitivity.csv',sensitivity);frame('B_suppression_recovery.csv',process)
    for src,dest,group in [('B_event_metrics.csv','B_metric_summary.csv',['subset','model','xc_km']),('B_length_metrics.csv','B_length_summary.csv',['subset','model'])]:
        df=pd.read_csv(OUT/'tables'/src);colsnum=[c for c in df.select_dtypes('number').columns if c not in group]
        g=df.groupby(group)[colsnum].agg(['mean','median','min','max']);g.columns=['_'.join(c) for c in g.columns];g.to_csv(OUT/'tables'/dest,encoding='utf-8-sig')
    save(OUT/'checkpoint.json',dict(stage='B_pilot_and_prototype_complete_pending_report',A_completed=True,B_cases=8,test_results_opened=0))
    print('B metrics complete',flush=True)

if __name__=='__main__':main()
