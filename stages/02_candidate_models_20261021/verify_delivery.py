"""Numerical/source integrity and export checks; visual review is explicit."""
from pathlib import Path
import argparse,hashlib,json,re,sys
from datetime import datetime
import numpy as np
import pandas as pd
from PIL import Image
from campaign import OUT,ROOT,read,save,sha,load_np,csvwrite
import models

def main(reviewed=False,update_route=False):
    manifest=read(OUT/'provenance/input_manifest.json');changed=[]
    for record in manifest['inputs']:
        p=ROOT/record['path']
        if sha(p)!=record['sha256'] or p.stat().st_mtime_ns!=record['mtime_ns']:changed.append(record['path'])
    assert not changed
    freeze=read(OUT/'provenance/frozen_protocol.json')
    assert sha(OUT/'PROTOCOL.md')==freeze['protocol_sha256']
    assert sha(OUT/'models.py')==freeze['model_source_sha256']
    meta=read(OUT/'data/metadata.json');ids={v['event_id'] for v in meta['development']}
    assert len(ids)==70 and not ids.intersection(manifest['validation_event_ids_not_loaded'])
    assert {v['event_id'] for v in meta['controls']}<=ids
    fits=[read(OUT/'fits'/f'model{m}_final.json') for m in range(3)]
    df=pd.read_csv(OUT/'tables/event_length_metrics.csv');assert len(df)==6120
    assert np.isfinite(df[['NSE','KGE','RMSE']].to_numpy()).all()
    for m in range(3):
        value=df[(df.dataset=='development')&(df.model==f'RWEC-{m}')].NSE.mean()
        assert abs(value-(1-fits[m]['mean_loss']))<1e-10
        assert fits[m]['charged_reach_nfev']<=49600
    assert len(pd.read_csv(OUT/'tables/control_response_metrics.csv'))==2880
    d=load_np(OUT/'data/development.npz');pred=load_np(OUT/'predictions/development.npz')['flow']
    imp=pd.read_csv(OUT/'tables/improvements.csv');chosen=read(OUT/'provenance/representative_events.json')
    representative=[]
    for k,eid in enumerate(chosen['preset']+[chosen['worst_additional']]):
        j=next(j for j,r in enumerate(meta['development']) if r['event_id']==eid);N=int(d['ns'][j])
        L=200 if k<3 else int(imp[(imp.dataset=='development')&(imp.model=='RWEC-2')&(imp.event_id==eid)].sort_values('delta_NSE').iloc[0].L_km)
        li=L//10-1
        for i in range(N):representative.append(dict(event_id=eid,L_km=L,time_h=i-48,
            input_flow_m3s=d['qin'][li,j,i],outlet_stage_m=d['h'][j,i],HEC_RAS=d['obs'][j,i],
            RWEC_0=pred[0,li,j,i],RWEC_1=pred[1,li,j,i],RWEC_2=pred[2,li,j,i]))
    csvwrite(OUT/'tables/representative_hydrographs.csv',representative)
    common=[read(p) for p in (OUT/'fits').glob('common_L*.json')];costs=[]
    for m,f in enumerate(fits):
        elapsed=sum(v['elapsed_s'] for v in common)+sum(v['elapsed_s'] for c in f['cycles'] for v in c['length_fits'])
        elapsed+=sum(c['global_fit']['elapsed_s'] for c in f['cycles'] if c['global_fit'])
        costs.append(dict(model=f'RWEC-{m}',free_parameters=60+m,budget_cap_reach_nfev=49600,
            actual_reach_nfev=f['charged_reach_nfev'],common_charged_reach_nfev=f['common_charged_reach_nfev'],
            additional_reach_nfev=f['additional_reach_nfev'],search_elapsed_sum_s=elapsed,
            timing_note='Sum of search wall times with shared initialization; parallel overlap prevents interpreting as process wall time'))
    csvwrite(OUT/'tables/model_costs.csv',costs)
    j=next(j for j,r in enumerate(meta['development']) if r['event_id']=='E007');N=int(d['ns'][j]);q=d['qin'][19,j,:N];h=d['h'][j,:N]
    causality=[]
    for m in range(3):
        pars=models.decode(fits[m]['x'][19]);shared=(0,0) if m==0 else models.decode_shared(fits[m]['shared_x'],m)
        original=models.route(q,h,N,*pars,*shared,200000.)[0]
        q1=q.copy();h1=h.copy();q1[100:]+=500;h1[100:]+=1
        perturbed=models.route(q1,h1,N,*pars,*shared,200000.)[0]
        error=float(abs(original[:100]-perturbed[:100]).max());assert error==0
        causality.append(dict(model=f'RWEC-{m}',future_input_perturbed_from_index=100,maximum_prior_prediction_difference=error))
    save(OUT/'provenance/causality_check.json',dict(pass_all=True,checks=causality,target_not_an_argument_to_route=True))
    exports=[]
    pdfs=sorted((OUT/'figures').glob('*.pdf'));assert len(pdfs)==9
    for p in pdfs:
        raw=p.read_bytes();fonts=sorted({v.decode('ascii') for v in re.findall(rb'/BaseFont\s*/([^\s/<>]+)',raw)})
        assert fonts and all('TimesNewRoman' in f for f in fonts),fonts
        box=re.search(rb'/MediaBox\s*\[\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)',raw)
        w=(float(box[3])-float(box[1]))*25.4/72;assert abs(w-180)<.01
        with Image.open(p.with_suffix('.tiff')) as im:dpi=[float(v) for v in im.info['dpi']]
        assert min(dpi)>=599
        assert '<text' in p.with_suffix('.svg').read_text(encoding='utf-8')
        exports.append(dict(figure=p.stem,width_mm=w,embedded_fonts=fonts,tiff_dpi=dpi))
    save(OUT/'provenance/export_verification.json',exports)
    preflight=read(OUT/'provenance/source_preflight.json')
    assert not [v for v in preflight['findings'] if v['level']=='FAIL' and v['check_id']!='FONT-FAMILY']
    qa=read(OUT/'provenance/figure_qa.json');assert len(qa['panels'])==33
    qa.update(source_font_override='FONT-FAMILY only: explicit user Times New Roman / SimSun requirement supersedes generic skill default; all embedded PDF fonts independently verified.',
        warning_review='LOG-GUARD: source asserts all k and n values positive before logarithmic axes; data positivity verified.',
        visual_review='33 panels and9 complete figures inspected; labels, units, legends and collisions checked' if reviewed else 'pending',
        source_integrity='All102original result NPZ files unchanged after analysis')
    if reviewed:
        for panel in qa['panels']:panel.update(visual_review='passed',units_legends_and_collisions='checked')
    save(OUT/'provenance/figure_qa.json',qa)
    save(OUT/'provenance/verification.json',dict(all_numerical_checks_pass=True,validation_targets_loaded=0,
        unchanged_sources=102,finite_prediction_metric_rows=6120,response_metric_rows=2880,
        independent_validation_claimed=False,maximum_recurrence_residual_m3=float(df.max_internal_balance_residual_m3.max()),
        timestep_portability_passed=False,primary_optimizer_budget_exhaustion_not_global_convergence=True,
        causal_prediction_check_passed=True,figures=9,visual_review_completed=reviewed))
    completed=datetime.now().astimezone().isoformat()
    status='first_round_complete_waiting_user_decision' if reviewed else 'numerical_complete_visual_review_pending'
    delivery=dict(status=status,completed=completed,development_events=70,test_events_used=0,controls=32,
        candidate_models=['RWEC-0','RWEC-1','RWEC-2'],figures=9,primary_report='REPORT.md',methods='METHODS.md',
        candidate_frozen_for_testing=False,additional_models_created=False,original_data_unchanged=True,
        scientific_conclusion='Linear boundary term improves development fit; saturation advantage not robust or identifiable; recovery and timestep limitations remain.',
        source_sha256={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.suffix in ['.py','.md']})
    save(OUT/'delivery.json',delivery);save(OUT/'checkpoint.json',dict(status=status,updated=completed,test_results_read=False))
    if update_route:
        routepath=ROOT.parents[1]/'CURRENT_EXPERIMENT.json';old=routepath.read_bytes();route=json.loads(old.decode('utf-8-sig'))
        assert Path(route['active_root'])==ROOT
        before=OUT/'provenance/current_experiment_before_stage02.json'
        if not before.exists():before.write_bytes(old)
        route.update(updated=completed,analysis_stage=status,analysis_latest=str(OUT),analysis_findings=str(OUT/'REPORT.md'),
            model_comparison_development_events=70,model_comparison_test_events=0,next_action_requires_user_decision=True)
        save(routepath,route)
    print(json.dumps(delivery,ensure_ascii=False,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--visual-review-confirmed',action='store_true');p.add_argument('--update-route',action='store_true')
    a=p.parse_args();main(a.visual_review_confirmed,a.update_route)
