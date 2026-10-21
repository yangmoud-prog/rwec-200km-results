"""Read-only hydrodynamic screening and event-grouped data preparation."""
from pathlib import Path
import csv, hashlib, json, shutil
from datetime import datetime
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
def read(p): return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x): Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def frame(p,x): pd.DataFrame(x).to_csv(p,index=False,encoding='utf-8-sig')
def phase(t,e):
    if t<e['tb_h']: return '水位上升前'
    if t<e['tb_h']+e['Tr_h']: return '水位上涨'
    if t<e['tb_h']+e['Tr_h']+e['Th_h']: return '高位维持'
    if t<e['Hend_h']: return '水位消落'
    return '边界恢复后'
def negative_duration(t,q):
    total=0.
    for a,b,qa,qb in zip(t[:-1],t[1:],q[:-1],q[1:]):
        if qa<0 and qb<0: total+=b-a
        elif qa<0: total+=(b-a)*(-qa)/(qb-qa)
        elif qb<0: total+=(b-a)*(-qb)/(qa-qb)
    return total

def run():
    for x in ['data','fits','tables','figures','animation','provenance','predictions']:(OUT/x).mkdir(exist_ok=True)
    plan=read(ROOT/'config/formal_plan.json')
    cases=sorted([c for c in plan['formal_cases'] if c['split']=='calibration'],key=lambda c:(c['role'],c['event_id'],c['case']))
    assert len(cases)==172
    eventids=sorted({c['event_id'] for c in cases if c['role']=='main'})
    assert len(eventids)==70
    permutation=np.random.default_rng(20261021).permutation(eventids)
    train=sorted(permutation[:49]); validation=sorted(permutation[49:])
    save(OUT/'provenance/event_split.json',dict(seed=20261021,train=train,validation=validation,
        controls_grouped_by_background_event=True,test_event_count=30,test_results_opened=0,
        created=datetime.now().astimezone().isoformat(),protocol_sha256=sha(OUT/'PROTOCOL.md')))
    manifest=[];raw=[];meta=[];allrows=[];episodes=[]
    for c in cases:
        p=ROOT/'results'/c['case']/'reach_timeseries.npz'
        with np.load(p,allow_pickle=False) as z: d={k:z[k] for k in z.files}
        qual=read(p.parent/'quality.json')
        assert qual['basic_quality_pass'] and qual['complete_duration'] and qual['tail_recovered']
        t=d['time_since_flood_h'];q=d['flow_m3s'];e=c['event']
        assert np.all(np.isfinite(q)) and np.all(np.isfinite(d['stage_m']))
        assert np.allclose(np.diff(t),1,atol=.02/3600,rtol=0)
        assert np.array_equal(d['x_from_outlet_km'],np.arange(0,201,10))
        assert abs(t[0]+48)<1e-5 and t[-1]>=c['end_h']-.02/3600
        for j in range(1,20):
            v=q[:,j];neg=v<0;ii=int(np.argmin(v))
            allrows.append(dict(case=c['case'],event_id=c['event_id'],role=c['role'],kind=c['kind'],
                subset='训练' if c['event_id'] in train else '开发内验证',x_km=j*10,
                n_saved=len(t),time_start_h=float(t[0]),time_end_h=float(t[-1]),
                min_q_m3s=float(v[ii]),min_time_h=float(t[ii]),min_phase=phase(t[ii],e),
                negative_samples=int(neg.sum()),negative_duration_linear_h=negative_duration(t,v),
                first_negative_h=float(t[np.where(neg)[0][0]]) if neg.any() else None,
                last_negative_h=float(t[np.where(neg)[0][-1]]) if neg.any() else None))
            start=np.flatnonzero(neg & ~np.r_[False,neg[:-1]])
            end=np.flatnonzero(neg & ~np.r_[neg[1:],False])
            for a,b in zip(start,end):
                episodes.append(dict(case=c['case'],event_id=c['event_id'],x_km=j*10,
                    first_saved_h=float(t[a]),last_saved_h=float(t[b]),saved_negative_samples=int(b-a+1),
                    min_q_m3s=float(v[a:b+1].min()),first_phase=phase(t[a],e),last_phase=phase(t[b],e)))
        raw.append(d)
        meta.append(dict(case=c['case'],event_id=c['event_id'],kind=c['kind'],role=c['role'],event=e,
            subset='train' if c['event_id'] in train else 'validation'))
        manifest.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns))
    detail=pd.DataFrame(allrows);frame(OUT/'tables/section_case_screen.csv',allrows)
    frame(OUT/'tables/negative_episodes.csv',episodes)
    sections=[]
    for x,g in detail.groupby('x_km'):
        worst=g.loc[g.min_q_m3s.idxmin()]
        sections.append(dict(x_km=int(x),cases=len(g),min_q_m3s=float(g.min_q_m3s.min()),
            negative_samples=int(g.negative_samples.sum()),negative_cases=int((g.negative_samples>0).sum()),
            negative_duration_case_sum_h=float(g.negative_duration_linear_h.sum()),
            worst_case=worst['case'],worst_phase=worst['min_phase'],
            passed=bool((g.negative_samples==0).all()),upstream_available_km=200-int(x)))
    frame(OUT/'tables/section_screen.csv',sections)
    selected=next(s['x_km'] for s in sections if s['x_km'] in [20,50,100] and s['passed'])
    assert selected==100
    lookup={m['case']:i for i,m in enumerate(meta)}
    pairs=[]
    for eid in eventids:
        vi=next(i for i,m in enumerate(meta) if m['event_id']==eid and m['kind']=='var' and m['role']=='main')
        ri=next(i for i,m in enumerate(meta) if m['event_id']==eid and m['kind']=='ref' and m['role']=='main')
        v,r=raw[vi],raw[ri]
        assert v['flow_m3s'].shape==r['flow_m3s'].shape
        assert np.allclose(v['time_since_flood_h'],r['time_since_flood_h'],atol=.02/3600,rtol=0)
        for j in range(1,20):
            dq=v['flow_m3s'][:,j]-r['flow_m3s'][:,j];dz=v['stage_m'][:,j]-r['stage_m'][:,j]
            imin=int(np.argmin(dq));imax=int(np.argmax(dq));iz=int(np.argmax(dz))
            pairs.append(dict(event_id=eid,subset=meta[vi]['subset'],x_km=j*10,
                min_var_minus_ref_q_m3s=float(dq[imin]),min_time_h=float(v['time_since_flood_h'][imin]),
                max_var_minus_ref_q_m3s=float(dq[imax]),max_time_h=float(v['time_since_flood_h'][imax]),
                rms_q_difference_m3s=float(np.sqrt(np.mean(dq*dq))),
                max_stage_increment_m=float(dz[iz]),max_stage_time_h=float(v['time_since_flood_h'][iz]),
                endpoint_q_difference_m3s=float(dq[-1]),max_abs_q_difference_m3s=float(np.max(abs(dq)))))
    frame(OUT/'tables/paired_hydrodynamic_response.csv',pairs)
    p=pd.DataFrame(pairs)
    response=p.groupby('x_km').agg({c:['mean','median','min','max'] for c in ['min_var_minus_ref_q_m3s',
        'max_var_minus_ref_q_m3s','rms_q_difference_m3s','max_stage_increment_m','max_abs_q_difference_m3s']})
    response.columns=['_'.join(c) for c in response.columns];response.reset_index().to_csv(OUT/'tables/response_summary.csv',index=False,encoding='utf-8-sig')
    nmax=max(len(d['time_since_flood_h']) for d in raw);nc=len(raw)
    qin=np.zeros((10,nc,nmax));obs=np.zeros((nc,nmax));h=obs.copy();ns=np.zeros(nc,dtype=np.int64)
    for i,d in enumerate(raw):
        N=len(d['time_since_flood_h']);ns[i]=N
        qin[:,i,:N]=d['flow_m3s'][:,11:].T;qin[:,i,N:]=qin[:,i,N-1,None]
        obs[i,:N]=d['flow_m3s'][:,10];obs[i,N:]=obs[i,N-1]
        h[i,:N]=d['stage_m'][:,0];h[i,N:]=h[i,N-1]
    assert qin.min()>0 and obs.min()>0
    den=np.array([np.sum((obs[i,:n]-obs[i,:n].mean())**2) for i,n in enumerate(ns)])
    np.savez_compressed(OUT/'data/all.npz',qin=qin,obs=obs,h=h,ns=ns,denom=den)
    save(OUT/'data/metadata.json',meta)
    save(OUT/'provenance/input_manifest.json',dict(inputs=manifest,plan_sha256=sha(ROOT/'config/formal_plan.json'),
        test_results_read=0,all_saved_hourly_samples_checked=True,negative_tolerance=0.0,
        no_claim_about_unsaved_solver_steps=True))
    save(OUT/'selection.json',dict(main_xc_km=100,auxiliary_xc_km=150,
        input_xu_km=list(range(110,201,10)),routing_length_km=list(range(10,101,10)),
        external_stage_location_km=0,not_a_dam_location=True,screen_cases=172,
        selection_uses_rwec_scores=False,resolution_h=1,sections=sections,
        selected_response=response.loc[100].to_dict()))
    spec=next(Path(r'C:\Users\mryan\Downloads').glob('RWEC_%E7%AC%AC%E4%BA%8C%E9%98%B6%E6%AE%B5* (1).md'))
    shutil.copyfile(spec,OUT/'provenance/USER_SPEC.md')
    save(OUT/'checkpoint.json',dict(stage='section_screen_complete',selected_xc_km=100,test_read=0))
    print(json.dumps(dict(selected_xc=100,min_q=float(obs.min()),training_events=len(train),validation_events=len(validation),
        response=response.loc[100].to_dict()),ensure_ascii=False,indent=2))

if __name__=='__main__': run()
