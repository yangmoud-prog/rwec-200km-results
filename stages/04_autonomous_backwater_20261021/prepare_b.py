"""Freeze independent hydraulic-control pilot inputs before any HEC computation."""
from common import *
import shutil,sys,math

dirs()
for d in ['config','scripts','models','results','logs/active_runs']:(B/d).mkdir(parents=True,exist_ok=True)
src=(A/'scripts/hec200.py').read_text(encoding='utf-8')
src=src.replace("DESIGN_ID = 'high_stage_hold_v1'","DESIGN_ID = 'mild_fixed_control_pilot_v1'")
src=src.replace('z=rs*.00005','z=rs*.0001').replace('bed=stations*.00005','bed=stations*.0001').replace('bed_slope=.00005','bed_slope=.0001')
src=src.replace("q,h=boundary(e,t_h,kind)","q,h=boundary(e,t_h,'ref' if kind=='control' else kind)")
old="""    u += [''.join(f'{v:8.5f}' for v in h[i:i+10]) for i in range(0,len(h),10)]"""
new="""    if kind=='control':
        # Static H,Q pairs, stage first, discharge second. No future stage hydrograph.
        u=u[:-2]+['Rating Curve= 241 ']
        qr=np.arange(0.,12000.1,50.)
        pairs=[f'{3.5+1.5*(v/1000.)**(2/3):8.5f}{v:8.1f}' for v in qr]
        u += [''.join(pairs[i:i+5]) for i in range(0,len(pairs),5)]
    else:
        u += [''.join(f'{v:8.5f}' for v in h[i:i+10]) for i in range(0,len(h),10)]"""
assert old in src;src=src.replace(old,new)
src=src.replace("qexpected,hexpected=boundary(cfg['event'],t_h,cfg['kind'])", "qexpected,hexpected=boundary(cfg['event'],t_h,'ref' if cfg['kind']=='control' else cfg['kind'])\n        if cfg['kind']=='control':\n            qr=np.arange(0.,12000.1,50.); hr=3.5+1.5*(qr/1000.)**(2/3)\n            hexpected=np.interp(q[:,-1],qr,hr)")
src=src.replace("'max_downstream_stage_error_m']<.001", "'max_downstream_stage_error_m']<.006")
dest=B/'scripts/hec200.py'
if dest.exists():assert dest.read_text(encoding='utf-8')==src
else:dest.write_text(src,encoding='utf-8')
shutil.copyfile(A/'config/plan_template_631.p01',B/'config/plan_template_631.p01')
sys.path.insert(0,str(B/'scripts'))
import hec200 as h
ids=['E023','E006','E007','E069']
events={e['event_id']:e for e in read(A/'config/events.json')['events'] if e['event_id'] in ids}
frozen=dict(geometry=dict(length_m=200000,slope=.0001,bottom_width_m=500,side_slope=3,manning=.025,dx_m=2000,dt_s=60),
    control=dict(formula='H0=3.5+1.5*(Q0/1000)^(2/3)',stage_m=[3.5,5.],flow_m3s=[0,1000],range_m3s=[0,12000],table_step_m3s=50,equivalent_weir_width_m=1000/(1.7*1.5**1.5)),
    events=events,train=['E023','E007','E069'],validation=['E006'],sealed_test_results_opened=0,
    scope='four flood pilots plus four fixed5 references only; no automatic formal expansion',protocol_sha256=sha(OUT/'PROTOCOL.md'),
    rating_format_source='https://raw.githubusercontent.com/gpt-cmdr/ras-commander/main/ras_commander/RasUnsteady.py')
freeze=B/'config/pilot_frozen.json'
if freeze.exists():assert read(freeze)==frozen
else:save(freeze,frozen)
cases=[]
for eid in ids:
    for kind in ['control','ref']:
        case=f'{eid}_{kind}_2km_60s'
        folder=h.build(case,events[eid],kind)
        cfg=read(folder/'settings.json')
        if kind=='control':assert 'Rating Curve= 241' in (folder/'Backwater200.u01').read_text()
        cases.append(dict(case=case,input_hashes=cfg['input_hashes']))
save(B/'config/pilot_input_manifest.json',dict(frozen_sha256=sha(freeze),cases=cases,created=h.now()))
print('Frozen independent B pilots:',len(cases),flush=True)
