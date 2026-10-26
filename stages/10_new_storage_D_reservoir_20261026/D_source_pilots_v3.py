from pilots_D import *
def main():
 cfg=json.loads((ROOT/'config/pilot_protocol.json').read_text(encoding='utf-8'));cases=[]
 dump(ROOT/'config/second_physical_candidate.json',dict(reason='crest 3.5m produced weak backwater response in uniform control; compare crest 2m with the same initial release, storage and floods. Decision uses HEC physics only.',crest_m=2.,validation_NSE_used=False,max_candidates=2))
 for e in cfg['events']:
  for geom in ['uniform','complex']:cases.append(build_D(f"{e['event_id']}_{geom}_reservoir_v3",e,geom,'reservoir',end_h=552,crest=2.))
 def job(f):
  r=run(f)
  if r.get('exit_code')==0:_,quality=extract(f);return quality
  return r
 rows=[]
 with ThreadPoolExecutor(max_workers=2) as pool:
  fs={pool.submit(job,f):f for f in cases}
  for ft in as_completed(fs):
   try:r=ft.result()
   except Exception as ex:r=dict(case=fs[ft].name,error=repr(ex))
   rows.append(r);dump(ROOT/'logs/pilot_v3_progress.json',rows);print(json.dumps(r),flush=True)
 ps=[]
 for e in cfg['events']:
  for geom in ['uniform','complex']:ps.append(pair(ROOT/'models'/f"{e['event_id']}_{geom}_reservoir_v3",ROOT/'models'/f"{e['event_id']}_{geom}_fixed_v2"))
 dump(ROOT/'logs/pilot_v3_pairs.json',ps);print(json.dumps(ps),flush=True)
if __name__=='__main__':main()
