from analyze_D import *
from concurrent.futures import ThreadPoolExecutor,as_completed

def main():
 es=[event(qp=5000,duration=132,event_id='DP01'),event(qp=6500,duration=96,r=.22,event_id='DP02'),event(qp=3500,duration=168,r=.35,event_id='DP03'),event(qp=5000,duration=168,r=.28,event_id='DP04')];es[-1].update(secondary_ratio=.6,secondary_delay_h=96.)
 cases=[]
 for e in es:
  for geom in ['uniform','complex']:
   for bc in ['fixed','reservoir']:cases.append(build_D(f"{e['event_id']}_{geom}_{bc}_v2",e,geom,bc,end_h=552))
 dump(ROOT/'config/pilot_protocol.json',dict(events=es,cases=[f.name for f in cases],selection='HEC physical gates only; no model validation or NSE is inspected',geometry_source='synthetic engineering stress channel, not a surveyed real river',all_four_factorial=True,initial_Q=1000,initial_H=5))
 def job(f):
  result=run(f)
  if result.get('exit_code')==0:_,qual=extract(f);return qual
  return dict(case=f.name,failed=True,run=result)
 rows=[]
 with ThreadPoolExecutor(max_workers=2) as pool:
  futures={pool.submit(job,f):f for f in cases}
  for ft in as_completed(futures):
   try:row=ft.result()
   except Exception as ex:row=dict(case=futures[ft].name,failed=True,error=repr(ex))
   rows.append(row);dump(ROOT/'logs/pilot_progress.json',rows);print(json.dumps(row),flush=True)
 pairs=[]
 for e in es:
  for geom in ['uniform','complex']:
   dyn=ROOT/'models'/f"{e['event_id']}_{geom}_reservoir_v2";ref=ROOT/'models'/f"{e['event_id']}_{geom}_fixed_v2"
   if (dyn/'extracted.npz').exists() and (ref/'extracted.npz').exists():pairs.append(pair(dyn,ref))
 dump(ROOT/'logs/pilot_pairs.json',pairs);print(json.dumps(pairs),flush=True)
if __name__=='__main__':main()
