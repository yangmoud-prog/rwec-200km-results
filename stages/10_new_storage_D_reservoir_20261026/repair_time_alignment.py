"""Invalidate all pre-validation D fits and reconstruct exact frozen hourly grid.
No HEC computation, parameter bounds, sample split, or optimization budget changes.
"""
from common10 import *
import shutil,datetime
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0,str(D/'scripts'))
from analyze_D import extract,pair

def main():
 archive=P/'provenance/invalid_time_alignment_attempt';archive.mkdir(exist_ok=True)
 moved=[]
 for folder,patterns in [(P/'provenance',['fit_D*.json','D*_calibration_frozen.json','baseline_*.log','baseline_progress.json']),(P,['PARAMETERS_D*.json','D??_completed.json']),(P/'tables',['D??_calibration.csv','D??_event_metrics.csv','D??_solver_quality.csv','D??_dt_sensitivity.csv']),(P/'predictions',['D??_D*.npz'])]:
  for pattern in patterns:
   for f in folder.glob(pattern):
    assert f.resolve().is_relative_to(P.resolve())
    dest=archive/f.relative_to(P);dest.parent.mkdir(parents=True,exist_ok=True)
    if dest.exists():raise RuntimeError('Archive collision: '+str(dest))
    shutil.move(str(f),str(dest));moved.append(str(f.relative_to(P)))
 save(P/'provenance/time_alignment_correction.json',dict(time=datetime.datetime.now().astimezone().isoformat(),reason='Floating HEC days excluded nominal t=0, shifting training truth by one hour. D validation aborted at first training event, before validation metrics.',invalidated_files=moved,old_results_scientifically_invalid=True,HEC_rerun=False,split_changed=False,budget_changed=False,validation_used_to_tune=False))
 folders=[f for f in sorted((D/'models').iterdir()) if (f/'run.json').exists() and read(f/'run.json').get('exit_code')==0]
 with ThreadPoolExecutor(max_workers=2) as pool:
  rows=list(pool.map(lambda f:extract(f)[1],folders))
 cfg=read(D/'config/D_frozen.json');formal=[r for r in rows if r['case'].split('_')[0] in cfg['train']+cfg['validation']]
 save(D/'logs/formal_checkpoint.json',dict(total=80,finished=len(formal),numerical_pass=sum(r['numerical_gate'] for r in formal),rows=formal,time_axis_corrected=True));frame('D_formal_numerical_quality.csv',formal)
 pairs=[]
 for e in cfg['events']:
  for geom in cfg['geometries']:
   pairs.append(pair(D/'models'/f"{e['event_id']}_{geom}_reservoir_{cfg['version']}",D/'models'/f"{e['event_id']}_{geom}_fixed_{cfg['version']}"))
 frame('D_formal_pairs.csv',pairs)
 from prepare_formal import reanalyze_pilots
 reanalyze_pilots()
 from formal_baselines import prepare
 checked=[]
 for domain in ['DUF','DUR','DCF','DCR']:
  for event,d in zip(cfg['train']+cfg['validation'],prepare(domain,cfg['train']+cfg['validation'],cfg)):
   assert len(d['time'])==721 and d['time'][0]==0 and d['time'][-1]==720
   assert len((d['q'][1:,0]+d['q'][:-1,0])/2)==720
   checked.append(dict(domain=domain,event=event,nodes=len(d['time']),intervals=720,initial_inflow=float(d['inflow'][0])))
 frame('D_time_alignment_regression.csv',checked)
 save(P/'provenance/time_alignment_verified.json',dict(cases=len(checked),passed=True,HEC_extractions_rebuilt=len(rows),maximum_jitter_h=max(r['maximum_raw_hourly_time_jitter_h'] for r in rows)))
 print('TIME_ALIGNMENT_VERIFIED',len(checked),flush=True)

if __name__=='__main__':main()
