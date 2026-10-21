from common import *
import importlib.util,concurrent.futures,traceback
spec=importlib.util.spec_from_file_location('hec075',B/'scripts/hec200.py');hec=importlib.util.module_from_spec(spec);spec.loader.exec_module(hec)
manifest=read(B/'config/pilot_input_manifest.json');assert sha(B/'config/pilot_frozen.json')==manifest['frozen_sha256']
progress={}
def run(c):
    try:
        folder=B/'models'/c['case'];hec.run(folder);q=hec.extract(folder)
        status='quality_pass' if q['basic_quality_pass'] and q['initialization_pass'] and q['tail_recovered'] else 'quality_review'
        return dict(case=c['case'],event_id=c['event_id'],kind=c['kind'],status=status,quality=q)
    except Exception as e:return dict(case=c['case'],event_id=c['event_id'],kind=c['kind'],status='failed_preserved',error=repr(e),trace=traceback.format_exc())
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    fs={pool.submit(run,c):c for c in manifest['cases']}
    for f in concurrent.futures.as_completed(fs):
        v=f.result();progress[v['case']]=v
        save(B/'logs/pilot_checkpoint.json',dict(cases=progress,finished=len(progress),total=40,updated=now(),failed=sum(x['status']=='failed_preserved' for x in progress.values())))
        print(len(progress),v['case'],v['status'],flush=True)
