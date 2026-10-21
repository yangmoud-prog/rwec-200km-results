"""Frozen development-only candidate comparison. Raw sources are read-only."""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json,math,os,shutil,sys,time
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from scipy.optimize import differential_evolution,minimize
import models

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
LENGTHS=np.arange(10,201,10)*1000
SEED=20261021

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
    p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8');tmp.replace(p)
def load_np(p):
    with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def csvwrite(p,rows):
    with Path(p).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
def stamp():return datetime.now().astimezone().isoformat()

def prepare():
    for d in ['data','fits','tables','figures','provenance','predictions']:(OUT/d).mkdir(exist_ok=True)
    plan=read(ROOT/'config/formal_plan.json');cp=read(ROOT/'logs/formal_checkpoint.json')
    assert cp['formal_success_count']==232 and not cp['formal_failures']
    mains=[c for c in plan['formal_cases'] if c['role']=='main' and c['kind']=='var']
    cal=sorted([c for c in mains if c['split']=='calibration'],key=lambda c:c['event_id'])
    test=sorted(c['event_id'] for c in mains if c['split']!='calibration')
    assert len(cal)==70 and len(test)==30
    controls=[c for c in plan['formal_cases'] if c['role']=='controlled_addon']
    controls=sorted(controls,key=lambda c:(c['event_id'],c['case']))
    assert len(controls)==32 and {c['event_id'] for c in controls}<={c['event_id'] for c in cal}
    sources=[];metadata={}
    for tag,cases in [('development',cal),('controls',controls)]:
        raw=[];records=[]
        for c in cases:
            path=ROOT/'results'/c['case']/'reach_timeseries.npz';z=load_np(path)
            s=read(ROOT/'models'/c['case']/'settings.json');q=read(ROOT/'results'/c['case']/'quality.json')
            assert q['basic_quality_pass'] and q['initialization_pass'] and q['tail_recovered']
            assert np.allclose(np.diff(z['time_since_flood_h']),1,atol=.02/3600,rtol=0)
            assert np.array_equal(z['x_from_outlet_km'],np.arange(0,201,10))
            raw.append(z);records.append(dict(case=c['case'],event_id=c['event_id'],event=s['event']))
            sources.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path),mtime_ns=path.stat().st_mtime_ns))
        nmax=max(len(z['time_since_flood_h']) for z in raw);n=len(raw)
        qin=np.zeros((20,n,nmax));h=np.zeros((n,nmax));obs=np.zeros((n,nmax));ns=np.zeros(n,dtype=np.int64)
        for j,z in enumerate(raw):
            N=len(z['time_since_flood_h']);ns[j]=N
            qin[:,j,:N]=z['flow_m3s'][:,1:].T;qin[:,j,N:]=z['flow_m3s'][-1,1:,None]
            h[j,:N]=z['stage_m'][:,0];h[j,N:]=h[j,N-1]
            obs[j,:N]=z['flow_m3s'][:,0];obs[j,N:]=obs[j,N-1]
        denom=np.array([np.sum((obs[j,:N]-obs[j,:N].mean())**2) for j,N in enumerate(ns)])
        assert np.all(denom>0) and np.all(np.isfinite(qin))
        np.savez_compressed(OUT/'data'/f'{tag}.npz',qin=qin,h=h,obs=obs,ns=ns,denom=denom)
        metadata[tag]=records
    save(OUT/'data/metadata.json',metadata)
    save(OUT/'provenance/input_manifest.json',dict(created=stamp(),inputs=sources,
        plan_sha256=sha(ROOT/'config/formal_plan.json'),validation_event_ids_not_loaded=test,
        data_selection='70 calibration variable-boundary runs plus32 controlled runs; no reference-target mixing'))
    source=Path(r'E:\01-HUST\Nature\03\Code')
    for old,new in [('routing.py','reference_routing.py'),('routing_aut_complement_3p.py','reference_complement.py')]:
        shutil.copyfile(source/old,OUT/'provenance'/new)
    spec=Path(r'C:\Users\mryan\Downloads\RWEC_%E5%8A%A8%E6%80%81%E5%9B%9E%E6%B0%B4%E5%82%A8%E9%87%8F%E4%BF%AE%E6%AD%A3%E4%B8%8E%E5%88%9D%E6%AD%A5%E6%A8%A1%E5%9E%8B%E5%AF%B9%E6%AF%94_Codex%E6%8C%87%E4%BB%A4.md')
    shutil.copyfile(spec,OUT/'provenance/USER_SPEC.md')
    # Read the archived manuscript without importing an Office application or changing it.
    manuscript=source/'figure code/paper_final_20260831/source_manuscript/基于水能互补方程的河道洪水演算方法_清江笋溪河合并_最新结果版.docx'
    if manuscript.exists():
        import zipfile,xml.etree.ElementTree as ET
        with zipfile.ZipFile(manuscript) as z:tree=ET.fromstring(z.read('word/document.xml'))
        ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
        paras=[''.join(el.text or '' for el in p.iter() if el.tag.rsplit('}',1)[-1]=='t') for p in tree.findall('.//w:p',ns)]
        (OUT/'provenance/manuscript_extraction.txt').write_text('\n'.join(paras),encoding='utf-8')
        save(OUT/'provenance/manuscript_source.json',dict(sha256=sha(manuscript),source_file=str(manuscript),read_only=True))
    save(OUT/'provenance/frozen_protocol.json',dict(created=stamp(),protocol_sha256=sha(OUT/'PROTOCOL.md'),
        model_source_sha256=sha(OUT/'models.py'),base_bounds=dict(k=[1e-6,100],n=[.01,300],tau_h=[0,150]),
        alpha_bounds=[0,10],theta_bounds=[0,200],ell_min_m=1000,linear_limit_at_theta_zero=True,
        common_budget_per_length=1200,cycles=4,baseline_refine_per_length=320,
        candidate_shared_per_cycle=80,candidate_refine_per_length=240,max_reach_evaluations_per_model=49600,
        seed=SEED,parallel_length_workers=2,objective='equal-event mean1-NSE; equal-length mean for global parameters'))
    print('Prepared70development events and32controls;30test events withheld.',flush=True)

def import_reference():
    sys.dont_write_bytecode=True
    for name,filename in [('routing','reference_routing.py'),('reference_complement','reference_complement.py')]:
        spec=importlib.util.spec_from_file_location(name,OUT/'provenance'/filename)
        module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module)
    return sys.modules['reference_complement']

def preflight():
    ref=import_reference();d=load_np(OUT/'data/development.npz');q=d['qin'][9,0,:d['ns'][0]]
    h=d['h'][0,:len(q)];checks=[]
    for k,n,tau in [(0.12,8,12.35),(3,2.3,0),(0.01,.7,4.75),(20,64,48)]:
        old=ref.route_lagged_complement_single_k_solved_s0(q,k,tau,n_comp=n,dt_h=1,pre_lag_mode='steady_initial')
        new=models.route(q,np.full(len(q),5.),len(q),k,n,tau,0,0,100000.)
        err=float(np.max(abs(old['Q_sim']-new[0])))
        serr=float(np.max(abs(old['S_t']-new[1])))
        assert err<1e-5 and serr<1e-7,(err,serr)
        for alpha,theta,L,hh in [(1,0,100000.,np.full(len(q),5.)),(1,4,0,h),(0,4,100000.,h)]:
            test=models.route(q,hh,len(q),k,n,tau,alpha,theta,L)
            assert np.max(abs(test[0]-new[0]))<1e-7
        candidate=models.route(q,h,len(q),k,n,tau,.4,4,100000.)
        assert np.all(candidate[0]>=0) and np.all(candidate[1]>=0)
        assert np.max(abs(candidate[3]))<.001
        checks.append(dict(k=k,n=n,tau_h=tau,max_flow_error_m3s=err,max_storage_error_million_m3=serr,
            max_balance_residual_m3=float(np.max(abs(candidate[3])))))
    idx=np.arange(70,dtype=np.int64)
    models.loss(d['qin'][9],d['h'],d['obs'],d['ns'],d['denom'],idx,.12,8,12,.4,4,100000.)
    start=time.perf_counter()
    for _ in range(30):models.loss(d['qin'][9],d['h'],d['obs'],d['ns'],d['denom'],idx,.12,8,12,.4,4,100000.)
    ms=(time.perf_counter()-start)/30*1000
    save(OUT/'provenance/solver_preflight.json',dict(created=stamp(),pass_all=True,
        baseline_checks=checks,zero_increment_and_zero_length_identity=True,nonnegative_and_complement_balance=True,
        loss_ms_per_length=ms,source_consistency='S0 solved; averaged lagged inflow; old source copied read-only',
        physical_balance_caveat='Model state balances shifted and clipped input; full physical control volume also includes propagation delay volume.'))
    print('Preflight passed; mean objective time',round(ms,2),'ms per length.',flush=True)

class BudgetEnd(Exception):pass

def optimize(fn,x0,cap,seed,anchors=()):
    rng=np.random.default_rng(seed);x0=np.asarray(x0,dtype=float);d=len(x0)
    bestx=x0.copy();bestf=float('inf');calls=0;trace=[]
    def wrapped(x):
        nonlocal calls,bestx,bestf
        if calls>=cap:raise BudgetEnd()
        v=float(fn(np.clip(x,0,1)));calls+=1
        if v<bestf:bestf=v;bestx=np.array(x).copy();trace.append([calls,v])
        return v
    start=time.perf_counter();popn=24 if d==3 else 12
    pop=rng.uniform(0,1,(popn,d));pop[0]=x0
    for i,a in enumerate(anchors):
        if i+1<popn:pop[i+1]=np.clip(a,0,1)
    iters=max(0,int(cap*.75)//popn-1)
    try:
        differential_evolution(wrapped,[(0.,1.)]*d,init=pop,maxiter=iters,tol=0,atol=0,
            seed=seed,polish=False,updating='immediate',workers=1)
        minimize(wrapped,bestx,method='L-BFGS-B',bounds=[(0.,1.)]*d,
            options=dict(maxfun=max(1,cap-calls),maxiter=cap,ftol=1e-11,gtol=1e-6,maxls=12,eps=1e-6))
    except BudgetEnd:pass
    return dict(x=bestx.tolist(),loss=bestf,nfev=calls,elapsed_s=time.perf_counter()-start,trace=trace)

def fit():
    assert read(OUT/'provenance/solver_preflight.json')['pass_all']
    cfg=read(OUT/'provenance/frozen_protocol.json');assert sha(OUT/'PROTOCOL.md')==cfg['protocol_sha256']
    d=load_np(OUT/'data/development.npz');idx=np.arange(70,dtype=np.int64)
    def objective(i,x,shared=(0.,0.)):
        return models.loss(d['qin'][i],d['h'],d['obs'],d['ns'],d['denom'],idx,*models.decode(x),*shared,float(LENGTHS[i]))
    common=[]
    for i,L in enumerate(LENGTHS):
        p=OUT/'fits'/f'common_L{L//1000:03d}.json'
        if p.exists():v=read(p)
        else:
            starts=[models.encode(k,n,min(150,L/1000/6)) for k,n in [(.03,2),(.12,8),(.35,16),(1,32),(5,64)]]
            v=optimize(lambda x:objective(i,x),starts[1],1200,SEED+i,starts)
            save(p,v)
        common.append(v)
        print('COMMON',int(L//1000),'km',round(v['loss'],5),'NFEV',v['nfev'],flush=True)
    commonx=np.array([v['x'] for v in common]);commoncost=sum(v['nfev'] for v in common)
    for model in [0,1,2]:
        final=OUT/'fits'/f'model{model}_final.json'
        if final.exists():continue
        state=dict(model=model,x=commonx.tolist(),shared_x=[] if model==0 else ([0.] if model==1 else [0.,math.log1p(4)/math.log(201)]),
            common_charged_reach_nfev=commoncost,additional_reach_nfev=0,cycles=[],started=stamp())
        for cycle in range(4):
            cp=OUT/'fits'/f'model{model}_cycle{cycle+1}.json'
            if cp.exists():state=read(cp);continue
            xs=np.array(state['x']);sx=state['shared_x']
            shared=(0.,0.) if model==0 else models.decode_shared(sx,model)
            grecord=None
            if model:
                def glob(x):
                    ab=models.decode_shared(x,model)
                    return float(np.mean([objective(i,xs[i],ab) for i in range(20)]))
                anchors=[]
                for alpha,theta in [(.01,0),(.1,1),(.5,4),(1,10),(3,30)]:
                    anchors.append([alpha/10] if model==1 else [alpha/10,math.log1p(theta)/math.log(201)])
                grecord=optimize(glob,sx,80,SEED+10000+model*100+cycle,anchors)
                sx=grecord['x'];shared=models.decode_shared(sx,model)
                state['additional_reach_nfev']+=20*grecord['nfev']
            cap=320 if model==0 else 240
            records=[None]*20
            def one(i):
                return i,optimize(lambda x:objective(i,x,shared),xs[i],cap,SEED+20000+cycle*100+i,
                    [np.clip(xs[i]+np.array([.02,0,0]),0,1),np.clip(xs[i]+np.array([0,.02,0]),0,1)])
            with ThreadPoolExecutor(max_workers=2) as pool:
                for future in as_completed([pool.submit(one,i) for i in range(20)]):
                    i,rr=future.result();records[i]=rr
            state['x']=[rr['x'] for rr in records];state['shared_x']=sx
            state['additional_reach_nfev']+=sum(rr['nfev'] for rr in records)
            value=float(np.mean([rr['loss'] for rr in records]))
            state['cycles'].append(dict(cycle=cycle+1,mean_loss=value,alpha=shared[0],theta=shared[1],
                ell_m=None if shared[1]==0 else 200000/shared[1],global_fit=grecord,length_fits=records))
            state['updated']=stamp();save(cp,state)
            save(OUT/'checkpoint.json',dict(status='calibrating',model=model,cycle=cycle+1,mean_loss=value,updated=stamp()))
            print('MODEL',model,'CYCLE',cycle+1,'meanNSE',round(1-value,6),'shared',shared,flush=True)
        state['finished']=stamp();state['mean_loss']=state['cycles'][-1]['mean_loss']
        state['charged_reach_nfev']=commoncost+state['additional_reach_nfev'];assert state['charged_reach_nfev']<=49600
        save(final,state)
    save(OUT/'checkpoint.json',dict(status='primary_calibration_complete',updated=stamp(),test_results_read=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','preflight','fit']);a=p.parse_args()
    {'prepare':prepare,'preflight':preflight,'fit':fit}[a.action]()
