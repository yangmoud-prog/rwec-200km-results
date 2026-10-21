from common import *
import sys,time,concurrent.futures,traceback
sys.path.insert(0,str(B/'scripts'))
import hec200 as h
manifest=read(B/'config/pilot_input_manifest.json')
assert sha(B/'config/pilot_frozen.json')==manifest['frozen_sha256']
progress={}
def run(case):
    folder=B/'models'/case
    try:
        h.run(folder);q=h.extract(folder)
        status='quality_pass' if q['basic_quality_pass'] and q['initialization_pass'] and q['tail_recovered'] else 'quality_review'
        return dict(case=case,status=status,quality=q)
    except Exception as e:return dict(case=case,status='failed_preserved',error=repr(e),trace=traceback.format_exc())
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    futures={pool.submit(run,c['case']):c['case'] for c in manifest['cases']}
    for f in concurrent.futures.as_completed(futures):
        v=f.result();progress[v['case']]=v
        save(B/'logs/pilot_checkpoint.json',dict(cases=progress,finished=len(progress),total=8,updated=h.now(),formal_expansion_authorized=False))
        print(v['case'],v['status'],flush=True)
