"""Frozen reference-core calibration and conservative incremental candidates."""
from pathlib import Path
import argparse, importlib.util, json, math, time, sys
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from scipy.optimize import differential_evolution,minimize,minimize_scalar
from prepare import OUT,ROOT,read,save,sha
import models as hy
core=hy.core

def load():
    with np.load(OUT/'data/all.npz') as z:d={k:z[k] for k in z.files}
    return d,read(OUT/'data/metadata.json')

def basefit(joint=False):
    d,meta=load();idx=np.array([i for i,m in enumerate(meta) if m['subset']=='train' and (joint or m['kind']=='ref') and m['role']=='main'],dtype=np.int64)
    assert len(idx)==(98 if joint else 49)
    fits=[None]*10
    def one(l):
        path=OUT/'fits'/f'{"jointbase" if joint else "base"}_L{(l+1)*10:03d}.json'
        if path.exists():return l,read(path)
        calls=0
        def fun(x):
            nonlocal calls
            calls+=1
            return core.loss(d['qin'][l],d['h'],d['obs'],d['ns'],d['denom'],idx,*core.decode(x),0.,0.,float((l+1)*10000))
        pop=np.random.default_rng(20261021+l).uniform(0,1,(24,3))
        for j,(k,n) in enumerate([(.03,2),(.12,8),(.35,16),(1,32),(5,64)]):pop[j]=core.encode(k,n,(l+1)*10/6)
        t=time.perf_counter()
        de=differential_evolution(fun,[(0,1)]*3,init=pop,maxiter=65,seed=20261021+l,polish=False,tol=0,atol=0)
        calls1=calls
        refined=minimize(fun,de.x,method='L-BFGS-B',bounds=[(0,1)]*3,options=dict(maxfun=400,maxiter=100,ftol=1e-12,eps=1e-6))
        x=refined.x if refined.fun<de.fun else de.x
        f=dict(x=x.tolist(),k=core.decode(x)[0],n=core.decode(x)[1],tau_h=core.decode(x)[2],
            loss_before=float(de.fun),loss_after=float(min(de.fun,refined.fun)),de_nfev=calls1,total_nfev=calls,
            optimizer_success=bool(refined.success),message=str(refined.message),elapsed_s=time.perf_counter()-t,
            training='49 floods, both conditions; robustness baseline only' if joint else '49 fixed-stage reference floods; all increments share this core')
        save(path,f);return l,f
    # Warm the compiled objective before launching independent lengths.
    core.loss(d['qin'][0],d['h'],d['obs'],d['ns'],d['denom'],idx,.12,8,1,0.,0.,10000.)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for f in as_completed([pool.submit(one,l) for l in range(10)]):
            l,v=f.result();fits[l]=v;print('基础河长',10*(l+1),'NSE',round(1-v['loss_after'],6),flush=True)
    save(OUT/'fits'/('jointbase.json' if joint else 'base.json'),fits)
    if not joint:makebase(d,meta,fits)

def makebase(d,meta,fits):
    qr=np.zeros_like(d['qin']);ss=qr.copy();ib=qr.copy()
    for l,f in enumerate(fits):
        for j,N in enumerate(d['ns']):
            q,s,i,_=core.route(d['qin'][l,j],d['h'][j],N,f['k'],f['n'],f['tau_h'],0,0,(l+1)*10000.)
            qr[l,j,:N]=q;ss[l,j,:N]=s;ib[l,j,:N]=i
            qr[l,j,N:]=q[-1];ss[l,j,N:]=s[-1];ib[l,j,N:]=i[-1]
    np.savez_compressed(OUT/'predictions/base.npz',q=qr,s=ss,ib=ib)

def arrays(d,meta,tag):
    idx=np.array([i for i,m in enumerate(meta) if m['subset']==tag and m['role']=='main'],dtype=int)
    with np.load(OUT/'predictions/base.npz') as z:b={k:z[k] for k in z.files}
    qr=b['q'][:,idx].reshape(-1,b['q'].shape[-1]);obs=np.tile(d['obs'][idx],(10,1))
    ns=np.tile(d['ns'][idx],10);den=np.tile(d['denom'][idx],10)
    L=np.repeat(np.arange(1,11)*10000,len(idx))
    hh=np.tile(d['h'][idx],(10,1));scale=500*L[:,None]*np.maximum(hh-5,0)/1e6
    infl=b['ib'][:,idx].reshape(qr.shape)
    features=np.stack([np.broadcast_to(L[:,None]/100000,qr.shape),np.maximum(hh-5,0)/6,infl/6000],axis=-1)
    return dict(idx=idx,qr=qr,obs=obs,ns=ns,den=den,L=L,h=hh,scale=scale,features=features)

def mse(q,A):
    return float(np.mean([np.sum((q[j,:n]-A['obs'][j,:n])**2)/A['den'][j] for j,n in enumerate(A['ns'])]))

def increments():
    d,meta=load();A=arrays(d,meta,'train')
    empirical=OUT/'fits/empirical.json'
    if not empirical.exists():
        fits=read(OUT/'fits/base.json');idx=A['idx'];cost=0
        def fun(alpha):
            nonlocal cost
            cost+=1;result=0.
            for l,f in enumerate(fits):
                result+=core.loss(d['qin'][l],d['h'],d['obs'],d['ns'],d['denom'],idx.astype(np.int64),f['k'],f['n'],f['tau_h'],alpha,0,(l+1)*10000.)
            return result/10
        r=minimize_scalar(fun,bounds=(0,2),method='bounded',options=dict(xatol=1e-7,maxiter=100))
        f0=fun(0);alpha=float(r.x) if r.fun<f0 else 0.
        save(empirical,dict(alpha=alpha,loss=float(min(r.fun,f0)),nfev=cost,success=bool(r.success),core_frozen=True))
        print('经验修正',alpha,'训练NSE',1-min(r.fun,f0),flush=True)
    for mode in [1,2]:
        path=OUT/'fits'/f'P{mode}.json'
        if path.exists():continue
        def fun(x):
            alpha=x[0]*2;tb=math.exp(x[1]) if mode==2 else 0.
            q,*_=hy.response(A['qr'],A['scale']*alpha,A['ns'],1.,tb,mode,True)
            return mse(q,A)
        bounds=[(0,1)]+([(math.log(.1),math.log(120))] if mode==2 else [])
        r=differential_evolution(fun,bounds,maxiter=35,popsize=12,seed=20261021+mode,polish=True,tol=1e-8)
        save(path,dict(alpha=float(r.x[0]*2),tb_h=float(math.exp(r.x[1])) if mode==2 else None,
            loss=float(r.fun),nfev=int(r.nfev),success=bool(r.success),message=str(r.message),
            state_constrained=True,raw_original_not_claimed_valid=True,core_frozen=True))
        print('P',mode,'NSE',1-r.fun,'参数',r.x,flush=True)
    save(OUT/'checkpoint.json',dict(stage='physical_candidate_fits_complete',test_read=0))

def neural():
    d,meta=load();A=arrays(d,meta,'train');V=arrays(d,meta,'validation');records=[]
    # Verify analytic chain rule against independent central differences.
    th=np.random.default_rng(17).normal(0,.15,22);th[20]=-3.5;th[21]=math.log(12)
    aa={k:(v[:4].copy() if isinstance(v,np.ndarray) and v.shape[0]==len(A['ns']) else v) for k,v in A.items()}
    fun=lambda x:hy.neural_loss(x,aa['features'],aa['scale'],aa['qr'],aa['obs'],aa['ns'],aa['den'])
    value,g=fun(th);checks=[]
    for k in [0,11,12,17,20,21]:
        e=np.zeros(22);e[k]=1e-5;fd=(fun(th+e)[0]-fun(th-e)[0])/2e-5
        err=abs(fd-g[k]);assert err<1e-5,(k,fd,g[k]);checks.append(dict(parameter=k,analytic=float(g[k]),finite_difference=float(fd),absolute_error=float(err)))
    save(OUT/'provenance/gradient_check.json',dict(passed=True,checks=checks))
    for seed in [20261021,20261022]:
        path=OUT/'fits'/f'neural_{seed}.json'
        if path.exists():records.append(read(path));continue
        rng=np.random.default_rng(seed);theta=rng.normal(0,.15,22);theta[20]=math.log(.05/1.95);theta[21]=math.log(12.)
        bounds=[(-3,3)]*12+[(-5,5)]*4+[(-3,3)]*4+[(-8,3),(math.log(.1),math.log(120))]
        trace=[];count=0;t=time.perf_counter()
        def obj(x):return hy.neural_loss(x,A['features'],A['scale'],A['qr'],A['obs'],A['ns'],A['den'])
        def callback(x):
            nonlocal count
            count+=1
            if count%25==0:
                f,g=obj(x);trace.append(dict(iteration=count,penalized_train_loss=float(f)))
                print('小网络',seed,'迭代',count,'损失',round(f,6),flush=True)
        r=minimize(obj,theta,jac=True,method='L-BFGS-B',bounds=bounds,callback=callback,
            options=dict(maxiter=400,maxfun=520,ftol=1e-10,gtol=1e-7,maxls=20))
        result=dict(seed=seed,theta=r.x.tolist(),loss=float(r.fun),nfev=int(r.nfev),nit=int(r.nit),
            success=bool(r.success),message=str(r.message),trace=trace,elapsed_s=time.perf_counter()-t,
            architecture='3-4-1 tanh/sigmoid capacity coefficient; causal conservative state;21 NN+1 shared response time',
            regularization=1e-5,validation_used_for_optimizer=False)
        save(path,result);records.append(result)
        print('小网络结束',seed,'NSE训练（含惩罚）',1-r.fun,'状态',r.message,flush=True)
    # The gate aggregates lengths within each flood before counting flood improvements.
    report=[];base=mse(V['qr'],V)
    baseline_case=np.array([np.sum((V['qr'][j,:n]-V['obs'][j,:n])**2)/V['den'][j] for j,n in enumerate(V['ns'])]).reshape(10,-1)
    for rec in records:
        th=np.array(rec['theta']);coef,*_=hy.network(th,V['features'])
        q,b,act,raw=hy.response(V['qr'],V['scale']*coef,V['ns'],1,math.exp(th[21]))
        losses=np.array([np.sum((q[j,:n]-V['obs'][j,:n])**2)/V['den'][j] for j,n in enumerate(V['ns'])]).reshape(10,-1)
        byevent={}
        for col,i in enumerate(V['idx']):byevent.setdefault(meta[i]['event_id'],[]).append(float((baseline_case-losses)[:,col].mean()))
        gains=np.array([np.mean(v) for v in byevent.values()]);minimum=float(min(q[j,:n].min() for j,n in enumerate(V['ns'])))
        maxbalance=0.
        for j,n in enumerate(V['ns']):
            rr=np.diff(b[j,:n])*1e6-(V['qr'][j,1:n]-q[j,1:n])*3600
            maxbalance=max(maxbalance,float(np.max(abs(rr))))
        boot=np.random.default_rng(20261021).choice(gains,size=(10000,len(gains)),replace=True).mean(axis=1)
        lower=float(np.quantile(boot,.025));relative=float((baseline_case.mean()-losses.mean())/baseline_case.mean())
        passed=bool(gains.mean()>0 and lower>0 and relative>=.05 and np.median(gains)>0 and (gains>0).mean()>=.6 and minimum>=-1e-9 and maxbalance<.001)
        report.append(dict(seed=rec['seed'],mean_nse_gain=float(gains.mean()),median_nse_gain=float(np.median(gains)),
            min_event_gain=float(gains.min()),max_event_gain=float(gains.max()),improved_events=int((gains>0).sum()),
            validation_events=21,min_prediction_m3s=minimum,max_increment_balance_residual_m3=maxbalance,
            relative_sse_reduction=relative,paired_event_bootstrap_lower_95=lower,
            passed=passed,event_gains={k:float(np.mean(v)) for k,v in byevent.items()}))
    passed=all(r['passed'] for r in report)
    chosen=min(records,key=lambda r:r['loss'])['seed']
    save(OUT/'fits/neural_gate.json',dict(passed=passed,records=report,chosen_seed_by_training_loss=chosen,
        no_symbolic_search_if_false=True,baseline_validation_mean_nse=1-base))
    save(OUT/'checkpoint.json',dict(stage='neural_diagnostic_complete',symbolic_gate_passed=passed,test_read=0))
    print(json.dumps(dict(symbolic_gate=passed,records=report),ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['base','jointbase','increments','neural']);args=p.parse_args()
    {'base':basefit,'jointbase':lambda:basefit(True),'increments':increments,'neural':neural}[args.action]()
