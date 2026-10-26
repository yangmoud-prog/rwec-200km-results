from pilots_D import *
import hashlib
P=ROOT.parent/'high_stage_hold_2km_60s_v2/analysis/10_new_storage_D_reservoir_20261026'

def reanalyze_pilots():
 rows=[];pairs=[]
 for folder in sorted((ROOT/'models').glob('DP*_v[23]')):
  if (folder/'run.json').exists() and json.loads((folder/'run.json').read_text()).get('exit_code')==0:
   try:
    actual=folder.with_name(folder.name+'_full') if folder.with_name(folder.name+'_full').exists() else folder
    _,r=extract(actual);r['version']=folder.name[-2:];rows.append(r)
   except Exception as ex:rows.append(dict(case=folder.name,analysis_failed=repr(ex)))
 for folder in sorted((ROOT/'models').glob('DP*_reservoir_v[23]')):
  event,geom,*_=folder.name.split('_');ref=ROOT/'models'/f'{event}_{geom}_fixed_v2'
  actual=folder.with_name(folder.name+'_full') if folder.with_name(folder.name+'_full').exists() else folder
  if actual!=folder:ref=ref.with_name(ref.name+'_full')
  if (actual/'extracted.npz').exists() and (ref/'extracted.npz').exists():
   r=pair(actual,ref);r['version']=folder.name[-2:];pairs.append(r)
 import pandas as pd
 pd.DataFrame(rows).to_csv(P/'tables/D_pilot_numerical_quality.csv',index=False,encoding='utf-8-sig');pd.DataFrame(pairs).to_csv(P/'tables/D_pilot_pairs.csv',index=False,encoding='utf-8-sig')
 return rows,pairs

def complete_recessions():
 cfg=json.loads((ROOT/'config/pilot_protocol.json').read_text(encoding='utf-8'));cases={};notes=[]
 for e in sorted(cfg['events'],key=lambda v:v.get('development_split')=='validation'):
  for geom in ['uniform','complex']:
   for version in ['v2','v3']:
    f=ROOT/'models'/f"{e['event_id']}_{geom}_reservoir_{version}"
    if not (f/'run.json').exists():continue
    _,quality=extract(f)
    if quality['tail_Q_max_error_m3s']>1:
     notes.append(dict(original=f.name,tail_error=quality['tail_Q_max_error_m3s'],reason='extend numerical evaluation to include complete recovery; not a parameter or event replacement'))
     new=build_D(f.name+'_full',e,geom,'reservoir',end_h=720,crest=2. if version=='v3' else 3.5);cases[new.name]=new
     ref=build_D(f"{e['event_id']}_{geom}_fixed_v2_full",e,geom,'fixed',end_h=720);cases[ref.name]=ref
 dump(ROOT/'logs/recession_extensions.json',notes)
 def job(f):run(f);return extract(f)[1]
 with ThreadPoolExecutor(max_workers=2) as pool:
  for r in pool.map(job,cases.values()):print(json.dumps(r),flush=True)

def sensitivity(version):
 cfg=json.loads((ROOT/'config/pilot_protocol.json').read_text(encoding='utf-8'));e=cfg['events'][1];crest=2. if version=='v3' else 3.5
 cases=[build_D(f'DS_DP02_{geom}_{tag}_{version}',e,geom,'reservoir',dx=dx,dt=dt,end_h=720,crest=crest) for geom in ['uniform','complex'] for tag,dx,dt in [('dt30',2000,30),('dx1000',1000,60)]]
 def job(f):run(f);return extract(f)[1]
 rows=[]
 with ThreadPoolExecutor(max_workers=2) as pool:
  for row in pool.map(job,cases):rows.append(row);print(json.dumps(row),flush=True)
 metrics=[]
 for f,r in zip(cases,rows):
  ds=np.load(f/'extracted.npz');geom=r['geometry'];ref=np.load(ROOT/'models'/f'DP02_{geom}_reservoir_{version}'/'extracted.npz');nt=len(ref['time']);qq=np.array([np.interp(ref['x'],ds['x'],v) for v in ds['q'][:nt]]);zz=np.array([np.interp(ref['x'],ds['x'],v) for v in ds['z'][:nt]]);errq=np.sqrt(np.mean((qq-ref['q'])**2))/ref['q'].max();errh=abs(zz-ref['z']).max();metrics.append(dict(case=f.name,Q_RMSE_peak_fraction=float(errq),H_max_difference_m=float(errh),numerical_gate=r['numerical_gate'],sensitivity_gate=bool(errq<.01 and errh<.1 and r['numerical_gate'])))
 import pandas as pd
 pd.DataFrame(metrics).to_csv(P/'tables/D_sensitivity_quality.csv',index=False,encoding='utf-8-sig');dump(ROOT/'logs/sensitivity_quality.json',metrics);return all(r['sensitivity_gate'] for r in metrics)

def generate(version):
 # New development sample, fully fixed before formal run; no old sealed event is accessed.
 u=qmc.LatinHypercube(5,seed=2026102601).random(20);valid_indices=set(np.random.default_rng(2026102602).choice(20,6,replace=False).tolist());events=[]
 for i,v in enumerate(u):
  e=event(qp=3000+3500*v[0],duration=96+72*v[1],r=.20+.15*v[2],event_id=f'D{i+1:03d}');e.update(secondary_ratio=(.25+.5*v[3]) if i%2 else 0.,secondary_delay_h=60+84*v[4],development_split='validation' if i in valid_indices else 'train');events.append(e)
 cfg=dict(created=now(),version=version,geometry_frozen=True,physical_selection_only=True,validation_NSE_inspected=False,train=[e['event_id'] for e in events if e['development_split']=='train'],validation=[e['event_id'] for e in events if e['development_split']=='validation'],events=events,geometries=['uniform','complex'],boundaries=['fixed','reservoir'],crest_m=2. if version=='v3' else 3.5,reservoir_area_m2=60000000.,dx_m=2000,dt_s=60,output_h=1,end_h=720,sealed_tests_accessed=0,formal_HEC_cases=80)
 dump(ROOT/'config/D_frozen.json',cfg);dump(P/'provenance/D_frozen.json',cfg);return cfg

def formal(cfg):
 cases=[]
 for e in sorted(cfg['events'],key=lambda v:v['development_split']=='validation'):
  for geom in cfg['geometries']:
   for bc in cfg['boundaries']:cases.append(build_D(f"{e['event_id']}_{geom}_{bc}_{cfg['version']}",e,geom,bc,end_h=cfg['end_h'],crest=cfg['crest_m']))
 dump(ROOT/'config/formal_cases.json',[str(f) for f in cases]);rows=[]
 def job(f):
  rr=run(f)
  if rr.get('exit_code')==0:return extract(f)[1]
  return dict(case=f.name,numerical_gate=False,run=rr)
 with ThreadPoolExecutor(max_workers=2) as pool:
  fs={pool.submit(job,f):f for f in cases}
  for ft in as_completed(fs):
   try:r=ft.result()
   except Exception as ex:r=dict(case=fs[ft].name,numerical_gate=False,error=repr(ex))
   rows.append(r);dump(ROOT/'logs/formal_checkpoint.json',dict(total=80,finished=len(rows),numerical_pass=sum(v.get('numerical_gate',False) for v in rows),rows=rows));print(json.dumps(r),flush=True)
 pairs=[]
 for e in cfg['events']:
  for geom in cfg['geometries']:
   dyn=ROOT/'models'/f"{e['event_id']}_{geom}_reservoir_{cfg['version']}";ref=ROOT/'models'/f"{e['event_id']}_{geom}_fixed_{cfg['version']}"
   try:pairs.append(pair(dyn,ref))
   except Exception as ex:pairs.append(dict(event=e['event_id'],geometry=geom,error=repr(ex)))
 import pandas as pd
 pd.DataFrame(rows).to_csv(P/'tables/D_formal_numerical_quality.csv',index=False,encoding='utf-8-sig');pd.DataFrame(pairs).to_csv(P/'tables/D_formal_pairs.csv',index=False,encoding='utf-8-sig')

def main():
 complete_recessions()
 rows,pairs=reanalyze_pilots()
 candidates=[]
 for version in ['v2','v3']:
  rr=[r for r in rows if r.get('version')==version and r.get('boundary')=='reservoir'];pp=[r for r in pairs if r['version']==version and r['geometry']=='complex']
  fixed=[r for r in rows if r.get('boundary')=='fixed']
  valid=len(rr)==8 and len(pp)==4 and all(r['numerical_gate'] for r in rr+fixed) and max(r['endpoint_range_km'] for r in pp)>20 and max(r['max_suppression_20_80_m3s'] for r in pp)>100 and max(r['max_recovery_20_80_m3s'] for r in pp)>100
  score=np.mean([min(r['max_suppression_20_80_m3s'],r['max_recovery_20_80_m3s']) for r in pp]) if pp else 0
  candidates.append(dict(version=version,eligible=bool(valid),physical_strength=float(score)))
 dump(ROOT/'config/physical_selection.json',dict(candidates=candidates,rule='all pilot numeric and full-process positivity gates; identifiable motion and both storage/release; prefer greater physical suppression/recovery, no NSE used'))
 good=[r for r in candidates if r['eligible']]
 if not good:print('NO_PHYSICALLY_ELIGIBLE_CHANNEL: preserve pilots; formal not started',flush=True);return
 version=max(good,key=lambda r:r['physical_strength'])['version']
 if not sensitivity(version):print('SENSITIVITY_GATE_FAILED: formal not started',flush=True);return
 cfg=generate(version);formal(cfg)
if __name__=='__main__':main()
