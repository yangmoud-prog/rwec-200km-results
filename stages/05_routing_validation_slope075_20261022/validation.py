from common import *
from routing_core import forecast,oracle,initial_ratio,complementary,interp_history
from numba import njit
from scipy.optimize import differential_evolution,minimize,minimize_scalar
from scipy.integrate import cumulative_trapezoid
import sys,time,traceback,h5py,importlib.util

NODES=np.array([20,40,60,80,100,120,0]);COLS=(200-NODES)//20
MODEL_MODES={'RWEC-0':0,'RWEC-L':1,'RWEC-V':2,'RWEC-H':3,'RWEC-1_oracle':4}

@njit(cache=True)
def fast_base(inflow,dt,k,n,tau):
    nt=len(inflow);q=np.full((nt,11),inflow[0]);q[:,0]=inflow
    s=np.full(10,initial_ratio(k,n)*inflow[0]*dt*3600)
    for i in range(1,nt):
        for j in range(10):
            lag=.5*(interp_history(q[:,j],i,i-tau/dt,q[i,j])+interp_history(q[:,j],i,i-1-tau/dt,q[i,j]))*dt*3600
            sw=s[j]+lag;se=s[j]+k*lag;ss=complementary(sw,se,n)
            q[i,j+1]=(sw-ss)/(dt*3600);s[j]=ss
    return q

def hec_module(root):
    spec=importlib.util.spec_from_file_location('hec_validation',root/'scripts/hec200.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def data_A():
    split=read(P4/'provenance/event_split.json');train=set(split['train']);dev=train|set(split['validation'])
    plan=read(A/'config/formal_plan.json');hec=hec_module(A);groups={};inputs=[]
    for c in plan['formal_cases']:
        if c['event_id'] not in dev:continue
        path=A/'results'/c['case']/'reach_timeseries.npz';d=loadnp(path);t=np.round(d['time_since_flood_h'],6)
        inputs.append(dict(path=str(path.relative_to(WORK)),sha256=sha(path),mtime_ns=path.stat().st_mtime_ns))
        packet=dict(case=c['case'],kind=c['kind'],role=c['role'],time=t,truth=d['flow_m3s'][:,NODES//10],
            true_H=hec.boundary(c['event'],t,c['kind'])[1],event=c['event'],condition='fixed5_reference' if c['kind']=='ref' else ('independent_stage_main' if c['role']=='main' else 'independent_stage_control'))
        groups.setdefault(c['event_id'],dict(subset='train' if c['event_id'] in train else 'validation',cases=[]))['cases'].append(packet)
    for eid,g in groups.items():
        nt=max(len(c['time']) for c in g['cases']);g['time']=np.arange(nt,dtype=float)-48
        g['inflow']=hec.boundary(g['cases'][0]['event'],g['time'],'ref')[0]
    save(OUT/'provenance/A_source_manifest.json',dict(inputs=inputs,development_events=70,scenarios=sum(len(g['cases']) for g in groups.values()),sealed_test_results_opened=0))
    return groups

def data_B():
    frozen=read(B/'config/pilot_frozen.json');hec=hec_module(B);groups={};quality=[];errors=[]
    for eid in frozen['order']:
        pair={}
        for kind in ['control','ref']:
            folder=B/'models'/f'{eid}_{kind}_2km_60s'
            try:
                ql=read(B/'results'/folder.name/'quality.json');quality.append(dict(event=eid,kind=kind,
                  **{k:v for k,v in ql.items() if isinstance(v,(float,int,bool,str))},native_balance_max_fraction=float(max(ql['native_max_balance_residual_fraction_storage_by_length'][1:]))))
                with h5py.File(folder/'Backwater200.p01.hdf','r') as f:
                    b=f[hec.TS];xs=b['Cross Sections'];raw=b['Time'][:].astype(float);t=np.round((raw-raw[0])*24-48,6)
                    pair[kind]=dict(t=t,q=xs['Flow'][:].astype(float)[:,::-1],z=xs['Water Surface'][:].astype(float)[:,::-1])
                # Completed but quality-review cases remain explicitly labelled and evaluated.
                pair[kind]['quality_pass']=ql['basic_quality_pass'] and ql['initialization_pass'] and ql['tail_recovered']
            except Exception as exc:errors.append(dict(event=eid,kind=kind,error=repr(exc),action='retained failure; other events continue'))
        if 'control' not in pair:continue
        c=pair['control'];t=c['t'];e=frozen['events'][eid];inf=hec.boundary(e,t,'ref')[0]
        groups[eid]=dict(subset='train' if eid in frozen['train'] else 'validation',time=t,inflow=inf,cases=[dict(case=f'{eid}_control_2km_60s',kind='control',role='main',condition='fixed_control_dynamic',time=t,truth=c['q'][:,NODES//2],true_H=c['z'][:,0],event=e,quality_pass=c['quality_pass'])])
        arrays=dict(time=t,x=np.arange(0.,200.1,2),q_control=c['q'],z_control=c['z'],inflow=inf)
        if 'ref' in pair:
            r=pair['ref'];assert np.array_equal(t,r['t']);dz=c['z']-r['z'];y=r['z']-np.arange(101)[None,:]*2000*.000075
            da=(500+6*y)*dz+3*dz*dz;dv=cumulative_trapezoid(da,np.arange(101)*2000,axis=1,initial=0)
            arrays.update(q_ref=r['q'],z_ref=r['z'],delta_stage=dz,dv=dv,lb=np.array([extent(v,np.arange(0.,200.1,2)) for v in dz]))
        np.savez_compressed(OUT/'data'/f'B_{eid}.npz',**arrays)
    frame('B_HEC_quality.csv',quality);frame('B_HEC_failures.csv',errors)
    return groups

def weights(g):
    nt=len(g['time']);W=np.zeros((nt,7));U=W.copy();E=W.copy();ns=len(g['cases'])
    for c in g['cases']:
        n=len(c['time']);mask=c['time']>=0;y=c['truth'][mask];den=np.sum((y-y.mean(0))**2,axis=0)
        ids=np.flatnonzero(mask);w=1/den/ns
        W[ids]+=w;U[ids]+=y*w;E[ids]+=y*y*w
    mu=np.divide(U,W,out=np.zeros_like(U),where=W>0);const=float(np.mean(np.sum(E-mu*U,axis=0)))
    g.update(W=W,mu=mu,irreducible=const)

def result(d,params,alpha,mode,lookup,c=None,dt=1.):
    if mode==4:return oracle(d['inflow'],dt,params,alpha,c['true_H'],lookup)
    return forecast(d['inflow'],dt,params,alpha,mode,lookup)

def fit(groups,label,lookup):
    train={e:g for e,g in groups.items() if g['subset']=='train'}
    for g in groups.values():weights(g)
    bounds=[(np.log(.0001),np.log(20)),(np.log(.15),np.log(150)),(0.,8.)]
    def decode(z):return [float(np.exp(z[0])),float(np.exp(z[1])),float(z[2])]
    def loss(z,mode=0,alpha=0.):
        params=decode(z);scores=[]
        for g in train.values():
            if mode==4:
                vals=[]
                for c in g['cases']:
                    gd=dict(inflow=g['inflow'][:len(c['time'])]);r=result(gd,params,alpha,mode,lookup,c)
                    mask=c['time']>=0;o=c['truth'][mask];p=r[0][mask][:,COLS]
                    vals.append(np.mean(np.sum((o-p)**2,0)/np.sum((o-o.mean(0))**2,0)))
                scores.append(np.mean(vals));continue
            if mode==0:p=fast_base(g['inflow'],1.,*params)[:,COLS];pen=0
            else:
                r=result(g,params,alpha,mode,lookup);p=r[0][:,COLS];pen=10*(r[6].sum()+r[9])
            scores.append(np.mean(np.sum(g['W']*(p-g['mu'])**2,0))+g['irreducible']+pen)
        return float(np.mean(scores))
    file=OUT/'fits'/f'{label}_core_final.json'
    if file.exists():core=read(file)['core'];z=np.r_[np.log(core[:2]),core[2]]
    else:
        start=time.perf_counter();opt=differential_evolution(loss,bounds,seed=20261022,popsize=6,maxiter=26,tol=1.e-5,polish=False)
        refine=minimize(loss,opt.x,method='Nelder-Mead',bounds=bounds,options=dict(maxiter=160,xatol=1.e-5,fatol=1.e-7))
        z=refine.x if refine.fun<opt.fun else opt.x;core=decode(z)
        save(file,dict(core=core,names=['k','n','tau_per20km_h'],training_events=list(train),loss=loss(z),bounds=bounds,elapsed_s=time.perf_counter()-start,targets='all available training scenarios event-equally weighted, 6 sections plus outlet'))
    fits=[]
    for name,mode in MODEL_MODES.items():
        gainfile=OUT/'fits'/f'{label}_mode{mode}_final.json'
        if gainfile.exists():f=read(gainfile);fits.append(f);continue
        start=time.perf_counter()
        if mode==0:alpha=0.;val=loss(z)
        else:
            opt=minimize_scalar(lambda a:loss(z,mode,a),bounds=(0.,1.),method='bounded',options=dict(maxiter=20,xatol=5.e-4))
            zero=loss(z,mode,0.);alpha=0. if zero<=opt.fun else float(opt.x);val=min(zero,float(opt.fun))
        f=dict(model=name,mode=mode,core=core,alpha=alpha,loss=val,core_fixed=True,train_events=len(train),elapsed_s=time.perf_counter()-start,
             closure='frozen transplanted control assumption' if label=='A' and mode in (1,2,3) else ('future stage oracle only' if mode==4 else 'known physical B control' if label=='B' else 'none'))
        save(gainfile,f);fits.append(f);print(label,name,alpha,val,flush=True)
    # Equal bounded refinement budget for each model, all using only training data.
    for fixed in fits.copy():
        mode=fixed['mode'];jf=OUT/'fits'/f'{label}_joint_mode{mode}_final.json'
        if jf.exists():fits.append(read(jf));continue
        start=time.perf_counter();z0=np.r_[np.log(fixed['core'][:2]),fixed['core'][2]]
        if mode==0:
            opt=minimize(lambda x:loss(x),z0,method='Nelder-Mead',bounds=bounds,options=dict(maxfev=70,xatol=1.e-5,fatol=1.e-7));cand=dict(core=decode(opt.x),alpha=0.)
        else:
            x0=np.r_[z0,fixed['alpha']];bb=bounds+[(0.,1.)]
            opt=minimize(lambda x:loss(x[:3],mode,x[3]),x0,method='Nelder-Mead',bounds=bb,options=dict(maxfev=70,xatol=1.e-5,fatol=1.e-7));cand=dict(core=decode(opt.x[:3]),alpha=float(opt.x[3]))
        if opt.fun>fixed['loss']:cand=dict(core=fixed['core'],alpha=fixed['alpha'])
        f=dict(fixed,**cand);f.update(model=fixed['model']+'_joint',core_fixed=False,loss=min(fixed['loss'],float(opt.fun)),refinement_evaluations=int(opt.nfev),refinement_converged=bool(opt.success),elapsed_s=time.perf_counter()-start)
        save(jf,f);fits.append(f);print(label,f['model'],f['loss'],flush=True)
    return fits

def evaluate(groups,label,lookup,fits):
    rows=[];numerics=[];lengths=[];process=[];sens=[];fails=[]
    variants=fits.copy()
    if label=='A':
        for name in ['RWEC-L_persist5','RWEC-V_persist5']:
            variants.append(dict(fits[0],model=name,mode=0,alpha=0.,closure='known initial5m persistence, correction exactly zero'))
    for eid,g in groups.items():
        for f in variants:
            name=f['model'];mode=f['mode'];cached=None
            for c in g['cases']:
                try:
                    start=time.perf_counter();nt=len(c['time'])
                    if mode==4:r=result(dict(inflow=g['inflow'][:nt]),f['core'],f['alpha'],mode,lookup,c)
                    else:
                        if cached is None:cached=result(g,f['core'],f['alpha'],mode,lookup)
                        r=cached
                    elapsed=time.perf_counter()-start;mask=c['time']>=0;t=c['time'];pred=r[0][:nt]
                    for xc in TARGETS:
                        j=int(np.where(NODES==xc)[0][0]);node=int((200-xc)//20)
                        rows.append(dict(event=eid,subset=g['subset'],case=c['case'],condition=c['condition'],model=name,xc_km=int(xc),Lr_km=int(200-xc),quality_pass=c.get('quality_pass',True),**metrics(c['truth'][mask,j],pred[mask,node],t[mask])))
                    if mode==4 or cached is r and c is g['cases'][0]:
                        numerics.append(dict(event=eid,subset=g['subset'],model=name,case=c['case'] if mode==4 else 'shared_by_identical_inflow',nonconverged_steps=int(r[6].sum()),max_iterations=int(r[5].max()),min_prediction_m3s=float(r[0].min()),
                          max_mass_residual_m3=float(abs(r[7]).max()),max_full_balance_residual_m3=float(abs(r[8]).max()),out_of_lookup_steps=int(r[9]),fixedpoint_residual_m3s=float(r[10]),runtime_s=elapsed))
                    if mode!=4 and c is g['cases'][0]:
                        np.savez_compressed(OUT/'predictions'/f'{label}_{eid}_{name}.npz',time=g['time'],q=r[0],h=np.full(len(g['time']),5.) if 'persist5' in name else r[1],lb_km=r[2]/1000,delta_stage=r[3],iterations=r[5],failed=r[6])
                except Exception as exc:fails.append(dict(event=eid,model=name,case=c['case'],error=repr(exc),trace=traceback.format_exc()))
            # Dynamic length and suppression / recovery errors against HEC paired targets.
            if cached is not None and mode in (1,2,3):
                if label=='A':d=loadnp(P4/'data'/f'A_{eid}.npz');t=d['time_h'];lb=d['lb_reference_km'];actualq=d['flow_var_m3s'];refq=d['flow_ref_m3s'];obs_h=d['boundary_increment_m']+5
                else:
                    d=loadnp(OUT/'data'/f'B_{eid}.npz');t=d['time'];lb=d.get('lb');actualq=d['q_control'];refq=d.get('q_ref');obs_h=d['z_control'][:,0]
                if lb is not None:
                    n=min(len(lb),len(cached[2]));mask=t[:n]>=0;error=cached[2][:n][mask]/1000-lb[:n][mask]
                    lengths.append(dict(event=eid,subset=g['subset'],model=name,RMSE_km=float(np.sqrt(np.mean(error**2))),MAE_km=float(np.mean(abs(error))),bias_km=float(error.mean()),reference_peak_km=float(lb[:n][mask].max()),predicted_peak_km=float(cached[2][:n][mask].max()/1000),H_RMSE_m=float(np.sqrt(np.mean((cached[1][:n][mask]-obs_h[:n][mask])**2)))))
                if refq is not None:
                    refpred=fast_base(g['inflow'],1.,*f['core'])
                    for xc in TARGETS:
                        col=int(xc//(10 if label=='A' else 2));n=min(len(actualq),len(cached[0]));mask=t[:n]>=0
                        actual=actualq[:n,col]-refq[:n,col];estimate=cached[0][:n,(200-xc)//20]-refpred[:n,(200-xc)//20]
                        process.append(dict(event=eid,subset=g['subset'],model=name,xc_km=int(xc),truth_suppression_min_m3s=float(actual[mask].min()),truth_recovery_max_m3s=float(actual[mask].max()),predicted_suppression_min_m3s=float(estimate[mask].min()),predicted_recovery_max_m3s=float(estimate[mask].max()),difference_process_RMSE_m3s=float(np.sqrt(np.mean((estimate[mask]-actual[mask])**2)))))
        print(label,'evaluated',eid,flush=True)
    frame(f'{label}_event_metrics.csv',rows);frame(f'{label}_numerical_checks.csv',numerics);frame(f'{label}_length_metrics.csv',lengths);frame(f'{label}_suppression_recovery.csv',process);frame(f'{label}_model_failures.csv',fails)
    # All non-oracle primary models: causality and fixed-parameter time-step sensitivity on every validation event.
    causal=[]
    for eid,g in groups.items():
        if g['subset']!='validation':continue
        for f in fits[:5]:
            if f['mode']==4:continue
            base=result(g,f['core'],f['alpha'],f['mode'],lookup);cut=np.flatnonzero(g['time']>=24)[0]
            changed=g['inflow'].copy();changed[cut+1:]=1000
            other=forecast(changed,1.,f['core'],f['alpha'],f['mode'],lookup)
            causal.append(dict(event=eid,model=f['model'],prefix_difference_m3s=float(abs(base[0][:cut+1]-other[0][:cut+1]).max())))
            e=g['cases'][0]['event'];hec=hec_module(A if label=='A' else B)
            for dt in [.5,.25]:
                fine=np.arange(g['time'][0],g['time'][-1]+dt/2,dt);inf=hec.boundary(e,fine,'ref')[0];r=forecast(inf,dt,f['core'],f['alpha'],f['mode'],lookup);p=r[0][::round(1/dt)]
                for c in g['cases']:
                    n=len(c['time']);mask=c['time']>=0
                    for xc in TARGETS:
                        j=int(np.where(NODES==xc)[0][0]);node=(200-xc)//20
                        sens.append(dict(event=eid,model=f['model'],case=c['case'],condition=c['condition'],dt_h=dt,xc_km=int(xc),nonconverged_steps=int(r[6].sum()),difference_from_1h_RMSE_m3s=float(np.sqrt(np.mean((p[:n][mask,node]-base[0][:n][mask,node])**2))),**metrics(c['truth'][mask,j],p[:n][mask,node],c['time'][mask])))
    frame(f'{label}_timestep_sensitivity.csv',sens);frame(f'{label}_causality_checks.csv',causal)
    df=pd.DataFrame(rows)
    for name,keys in [('model_summary',['subset','condition','model']),('section_summary',['subset','condition','model','xc_km'])]:
        s=df.groupby(keys)[['NSE','RMSE_m3s','KGE','peak_error_pct','peak_time_error_h']].agg(['mean','median','min','max']);s.columns=['_'.join(c) for c in s.columns];s.to_csv(OUT/'tables'/f'{label}_{name}.csv',encoding='utf-8-sig')
    frame(f'{label}_parameters.csv',[dict(model=f['model'],k=f['core'][0],n=f['core'][1],tau_per20km_h=f['core'][2],alpha=f['alpha'],core_fixed=f['core_fixed'],training_loss=f['loss'],closure=f['closure']) for f in fits])

def main(label):
    dirs();groups=data_A() if label=='A' else data_B();lookup=loadnp(OUT/'data'/f'{label}_lookup.npz')
    fits=fit(groups,label,lookup);save(OUT/'fits'/f'{label}_all_final.json',fits);evaluate(groups,label,lookup,fits)
    save(OUT/'logs'/f'{label}_validation_checkpoint.json',dict(status='computed',events=len(groups),train=sum(g['subset']=='train' for g in groups.values()),validation=sum(g['subset']=='validation' for g in groups.values()),updated=now(),sealed_test_results_opened=0))
    print(label,'ACTUAL ROUTING VALIDATION COMPLETE',flush=True)
if __name__=='__main__':main(sys.argv[1])
