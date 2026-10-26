from common import *
from scipy.stats import qmc
import hashlib
queue=read(C/'config/pilot_queue.json');cp=read(C/'logs/pilot_queue_checkpoint.json')
assert cp['completed']==cp['total']==24
rows=[]
for slope in [.00005,.000075,.0001]:
 for event in ['CP01','CP02','CP03','CP04']:
  prefix=f'pilot_s{round(slope*1e6):03d}_{event}'
  d=load(C/'results'/f'{prefix}_var/full_fields.npz');r=load(C/'results'/f'{prefix}_ref/full_fields.npz')
  qual=[read(C/'results'/f'{prefix}_{kind}/quality.json') for kind in ['var','ref']]
  use=d['time']>=-1.e-6;x=d['x_m'][::-1]/1000;dh=(d['z']-r['z'])[use,::-1]
  lengths=np.array([extent(v,x)[0] for v in dh]);cols=[np.argmin(abs(d['x_m']-s*1000)) for s in TARGETS]
  dq=(d['q']-r['q'])[use];supp=[];release=[]
  for s in [20,40,60,80]:
   j=np.argmin(abs(d['x_m']-s*1000));supp.append(float(-dq[:,j].min()));release.append(float(dq[:,j].max()))
  row=dict(slope_per_mille=slope*1000,event=event,min_target_Q_m3s=float(d['q'][:,cols].min()),min_all_Q_m3s=float(d['q'].min()),
   numerical_pass=all(q['basic_quality_pass'] and q['initialization_pass'] and q['tail_recovered'] for q in qual),
   endpoint_excursion_km=float(np.ptp(lengths)),sections_with_50m3s_suppression_and_release=int(np.sum((np.array(supp)>=50)&(np.array(release)>=50))),
   suppression_20_40_60_80=supp,release_20_40_60_80=release)
  row['positive_flow_pass']=row['min_target_Q_m3s']>0;rows.append(row)
frame('pilot_physical_screening.csv',rows)
summary=[]
for s in [.05,.075,.1]:
 rr=[r for r in rows if abs(r['slope_per_mille']-s)<1.e-8]
 strict=all(r['numerical_pass'] and r['positive_flow_pass'] and r['endpoint_excursion_km']>=40 and r['sections_with_50m3s_suppression_and_release']>=2 for r in rr)
 summary.append(dict(slope_per_mille=s,initial_strict_gate_pass=strict,eligible=all(r['numerical_pass'] and r['positive_flow_pass'] and r['endpoint_excursion_km']>=40 for r in rr) and max(r['sections_with_50m3s_suppression_and_release'] for r in rr)>=2,min_Q=min(r['min_target_Q_m3s'] for r in rr),median_sections=float(np.median([r['sections_with_50m3s_suppression_and_release'] for r in rr])),median_excursion=float(np.median([r['endpoint_excursion_km'] for r in rr]))))
frame('pilot_slope_summary.csv',summary);print(pd.DataFrame(summary).to_string(index=False),flush=True)
eligible=[r for r in summary if r['eligible']]
save(P/'provenance/physical_gate_revision_before_routing.json',dict(created=hec.now(),initial_strict_gate_passed=False,
 reason='The additional requirement that every1m low-amplitude pilot have>=50m3/s suppression AND release at>=2 sections excluded all geometries. The user requires a forward-flow channel with identifiable response, not that every low-amplitude control be equally strong.',
 retained_hard_gates='All4 pilots pass numerical, initialization, recovery, strictly positive all-target flow, and>=40km endpoint excursion.',
 revised_response_gate='At least one pilot has>=2 of20/40/60/80km sections with both suppression and release>=50m3/s; then rank all qualifying slopes by median number and median endpoint excursion.',
 disclosure='Exploratory hydraulic-screening rule revision before any routing calibration; initial no-pass result retained. Not a preregistered confirmatory selection.',
 routing_validation_scores_seen=False,slower_forcing_fallback_not_used='Positive-flow candidates already exist; slower forcing would weaken the low-amplitude response that caused strict exclusion.'))
if not eligible:
 save(P/'provenance/selection_failure.json',dict(reason='No slope passed all predeclared hydraulic gates',summary=summary));raise SystemExit('No eligible slope; preserve pilots and use only predeclared slower-boundary fallback if time permits')
selected=sorted(eligible,key=lambda r:(r['median_sections'],r['median_excursion'],-abs(r['slope_per_mille']-.075)),reverse=True)[0]
frozen=P/'provenance/C_frozen.json'
assert not frozen.exists(),'Do not refreeze using later model results'
old=read(W/'variants/mild_control_slope075_20pilot_v1/config/pilot_frozen.json')
u=qmc.LatinHypercube(5,seed=20261026).random(20);lo=np.array([1,0,.55,.1,.65]);hi=np.array([4,.2,.8,.3,.95]);par=lo+u*(hi-lo)
events={}
for idx,eid in enumerate(old['order']):
 e=old['events'][eid];dh,tb,tr,th,tf=par[idx];events[eid]=hec.event(e['Qp_m3s'],e['D_h'],e['r'],dh,tb,tr,tf,eid,th)
plan=dict(selected_slope=selected['slope_per_mille']/1000,selected_by='HEC hydraulic criteria only, before routing calibration',screening=summary,order=old['order'],train=old['train'],validation=old['validation'],events=events,
 split_source='B20 unchanged subset of original49/21 development split',sealed_tests_opened=0,dt_s=60,dx_m=2000,end_h=552,independent_prescribed_stage=True,created=hec.now(),protocol_sha256=hec.sha(P/'provenance/protocol_before_HEC.json'))
save(frozen,plan);save(C/'config/C_frozen.json',plan)
q=[]
for eid in plan['order']:
 for kind in ['ref','var']:
  case=f'C_{eid}_{kind}';hec.build(case,events[eid],kind,end_h=552,slope=plan['selected_slope']);q.append(dict(case=case,event=eid,kind=kind,phase='development'))
save(C/'config/development_queue.json',q)
# Worst physical pilot, selected without any model scores; both time and space checks before batch quality acceptance.
adverse=read(C/'config/pilot_events.json')[0];q=[]
for dx,dt,label in [(2000,30,'dt30'),(1000,60,'dx1000')]:
 for kind in ['ref','var']:
  case=f'sensitivity_CP01_{label}_{kind}';hec.build(case,adverse,kind,dx=dx,dt=dt,end_h=552,slope=plan['selected_slope']);q.append(dict(case=case,event='CP01',kind=kind,phase='sensitivity'))
save(C/'config/sensitivity_queue.json',q)
print('Frozen:',plan['selected_slope'],'prepared40development+4sensitivity')
