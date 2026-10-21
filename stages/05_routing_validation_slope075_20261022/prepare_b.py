from common import *
import importlib.util,shutil
dirs()
for d in ['config','scripts','models','results','logs/active_runs']:(B/d).mkdir(parents=True,exist_ok=True)
split=read(P4/'provenance/event_split.json')
events={e['event_id']:e for e in read(A/'config/events.json')['events'] if e['event_id'] in split['train']+split['validation']}
selected=['E023','E006','E007','E069']
ids=sorted(events);vectors=np.array([[events[e][v] for v in ['Qp_m3s','D_h','r']] for e in ids]);vectors=(vectors-vectors.min(0))/(vectors.max(0)-vectors.min(0))
for subset,total in [('train',14),('validation',6)]:
    while sum(e in split[subset] for e in selected)<total:
        candidates=[e for e in split[subset] if e not in selected]
        distances={e:float(np.min(np.sum((vectors[ids.index(e)]-vectors[[ids.index(k) for k in selected]])**2,axis=1))) for e in candidates}
        selected.append(sorted(candidates,key=lambda e:(-distances[e],e))[0])
src=(WORK/'variants/mild_control_2km_60s_pilot_v1/scripts/hec200.py').read_text(encoding='utf-8')
src=src.replace("DESIGN_ID = 'mild_fixed_control_pilot_v1'","DESIGN_ID = 'mild_slope075_twenty_pilot_v1'")
src=src.replace('rs*.0001','rs*.000075').replace('stations*.0001','stations*.000075').replace('bed_slope=.0001','bed_slope=.000075')
src=src.replace('200 km trapezoidal channel. Gamma inflow; paired stage. High-stage hold=0.10D-0.30D.','200 km slope0.075per1000 channel; fixed hydraulic rating or constant5m reference.')
dest=B/'scripts/hec200.py'
if dest.exists():assert dest.read_text(encoding='utf-8')==src
else:dest.write_text(src,encoding='utf-8')
shutil.copyfile(A/'config/plan_template_631.p01',B/'config/plan_template_631.p01')
spec=importlib.util.spec_from_file_location('hec075',dest);hec=importlib.util.module_from_spec(spec);spec.loader.exec_module(hec)
frozen=dict(created_before_runs=True,slope=.000075,length_m=200000,width_m=500,side_slope=3,manning=.025,dx_m=2000,dt_s=60,
    initial_Q_m3s=1000,initial_H_m=5,control='H=3.5+1.5*(Q/1000)^(2/3)',target_stage_rise_m=[1,4],
    events={e:events[e] for e in selected},order=selected,train=[e for e in selected if e in split['train']],validation=[e for e in selected if e in split['validation']],
    n_floods=20,n_paired_HEC_cases=40,sealed_test_results_opened=0,selection='preserved4 then greedy farthest Qp,D,r, fixed14/6 original development partition',
    protocol_sha256=sha(OUT/'PROTOCOL.md'))
freeze=B/'config/pilot_frozen.json'
if freeze.exists():assert read(freeze)==frozen
else:save(freeze,frozen)
cases=[]
for eid in selected:
    for kind in ['control','ref']:
        case=f'{eid}_{kind}_2km_60s';folder=hec.build(case,events[eid],kind)
        cfg=read(folder/'settings.json');assert cfg['bed_slope']==.000075
        cases.append(dict(case=case,event_id=eid,kind=kind,input_hashes=cfg['input_hashes']))
save(B/'config/pilot_input_manifest.json',dict(frozen_sha256=sha(freeze),created=now(),cases=cases))
save(OUT/'provenance/B_frozen_selection.json',frozen)
print('20 paired floods /40 HEC projects frozen',flush=True)
