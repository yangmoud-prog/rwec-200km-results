"""Conditional, training-only PySR distillation. Independent recursion is evaluated later."""
import os,sys,json,math,time,traceback,socket,faulthandler
import urllib.request
from pathlib import Path
OUT=Path(__file__).resolve().parent
sys.stdout.reconfigure(line_buffering=True)
socket.setdefaulttimeout(60)
faulthandler.enable()
faulthandler.dump_traceback_later(180,repeat=True)
sys.path.insert(0,str(OUT/'.deps'))
RUNTIME=OUT.parents[3]/'.runtime'
os.environ['PYTHON_JULIAPKG_PROJECT']=str(RUNTIME/'sr_env')
os.environ['JULIA_DEPOT_PATH']=str(RUNTIME/'julia_depot')
if (RUNTIME/'julia/bin/julia.exe').exists():
    os.environ['PYTHON_JULIAPKG_EXE']=str(RUNTIME/'julia/bin/julia.exe')
os.environ['JULIA_NUM_THREADS']='2'
os.environ['JULIA_NUM_PRECOMPILE_TASKS']='2'
os.environ['JULIA_PKG_SERVER']='https://us-east.pkg.julialang.org'
# Julia's downloader does not automatically inherit Windows Internet Options.
# Reuse the configured system proxy without printing credentials or URLs.
for scheme,proxy in urllib.request.getproxies().items():
    if scheme in ['http','https']:os.environ[scheme+'_proxy']=proxy
from prepare import read,save,frame
from fit import arrays,load,hy
import numpy as np

def run():
    gate=read(OUT/'fits/neural_gate.json');assert gate['passed']
    seed=gate['chosen_seed_by_training_loss'];rec=read(OUT/'fits'/f'neural_{seed}.json');th=np.array(rec['theta'])
    d,meta=load();A=arrays(d,meta,'train');coef,_,_=hy.network(th,A['features'])
    # At zero external-stage increment the coefficient is unidentifiable and
    # multiplied by zero. Distill only its active, training-event domain.
    X=[];y=[];rows=[];rng=np.random.default_rng(20261021)
    for eid in sorted({meta[i]['event_id'] for i in A['idx']}):
        pool=[]
        for j,i in enumerate(A['idx']):
            if meta[i]['event_id']!=eid or meta[i]['kind']!='var':continue
            for l in range(10):
                row=l*len(A['idx'])+j
                for t in np.flatnonzero(A['scale'][row,:A['ns'][row]]>0):pool.append((row,int(t)))
        chosen=rng.choice(len(pool),min(80,len(pool)),replace=False)
        for k in chosen:
            row,t=pool[k];c=coef[row,t]
            X.append(A['features'][row,t]);y.append(math.log(c/(2-c)))
            rows.append(dict(event_id=eid,row=row,time_index=t,l=float(X[-1][0]),h=float(X[-1][1]),i=float(X[-1][2]),logit=float(y[-1])))
    X=np.array(X);y=np.array(y);frame(OUT/'tables/symbolic_training_sample.csv',rows)
    save(OUT/'provenance/symbolic_protocol.json',dict(training_events=49,validation_points=0,samples=len(X),
        per_event_max=80,seed=20261021,target='Neural capacity logit, not HEC target discharge',
        envelope='coefficient=2/(1+exp(-g(l,h,i))); B_eq=500*L*max(Hd-5,0)*coefficient',
        variables=dict(l='Lr/100000 m',h='max(Hd-5m,0)/6m',i='lagged step-average inflow/6000 m3s'),
        maxsize=12,operators=['+','-','*','/'],niterations=60,timeout_search_seconds=120,
        response_time_h=math.exp(th[21]),selection='At most3 training-loss winners at complexity caps4,8,12. No validation selection.'))
    print('PySR初始化；训练采样点',len(X),flush=True)
    start=time.perf_counter()
    from pysr import PySRRegressor
    import pysr,importlib.metadata
    model=PySRRegressor(niterations=60,binary_operators=['+','-','*','/'],unary_operators=[],
        maxsize=12,populations=8,population_size=30,model_selection='accuracy',
        parallelism='serial',deterministic=True,random_state=20261021,
        timeout_in_seconds=120,verbosity=0,progress=False,
        output_directory=str(OUT/'fits/symbolic_search'),run_id='training_only',
        precision=64)
    model.fit(X,y,variable_names=['l','h','i'])
    table=model.equations_;table.to_csv(OUT/'tables/symbolic_pareto.csv',index=False,encoding='utf-8-sig')
    candidates=[];seen=set()
    for cap in [4,8,12]:
        eligible=table[table.complexity<=cap]
        if not len(eligible):continue
        row=eligible.loc[eligible.loss.idxmin()];equation=str(row.sympy_format)
        if equation in seen:continue
        seen.add(equation)
        candidates.append(dict(id=f'SR{len(candidates)+1}',complexity=int(row.complexity),equation_logit=equation,
            offline_logit_loss=float(row.loss),tb_h=math.exp(th[21]),neural_seed=seed,
            nonnegative_capacity_envelope=True,parameters_refitted_to_validation=False))
    assert len(candidates)>0
    save(OUT/'fits/symbolic_candidates.json',dict(status='completed',candidates=candidates,
        pysr_version=pysr.__version__,juliacall_version=importlib.metadata.version('juliacall'),
        elapsed_including_julia_initialization_s=time.perf_counter()-start,
        success_here_means_offline_search_only=True,independent_recurrence_verification_required=True))
    print(json.dumps(candidates,indent=2,ensure_ascii=False),flush=True)

if __name__=='__main__':
    try:run()
    except Exception as e:
        save(OUT/'fits/symbolic_failure.json',dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc(),
            symbolic_equation_not_claimed=True))
        raise
