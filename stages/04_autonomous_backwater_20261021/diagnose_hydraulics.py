"""Given-boundary diagnostics; A uses observed forcing, never autonomous forecast."""
from common import *
from hydraulics import profile
import sys
from scipy.integrate import cumulative_trapezoid
sys.path.insert(0,str(A/'scripts'))
import hec200

def run():
    manifest=read(OUT/'provenance/A_input_manifest.json')
    # Verify every reused result and derived input after model fitting.
    checks=[]
    for row in manifest['inputs']:
        p=A/row['path'];assert sha(p)==row['sha256'] and p.stat().st_mtime_ns==row['mtime_ns']
        checks.append(row['path'])
    save(OUT/'provenance/verification.json',dict(raw_and_derived_sources_unchanged=True,verified_files=len(checks),sealed_test_results_opened=0,
         experiment_B_freeze_sha256=sha(B/'config/pilot_frozen.json'),B_pilot_cases=8))
    events={e['event_id']:e for e in read(A/'config/events.json')['events']}
    split=read(OUT/'provenance/event_split.json');rows=[]
    x=np.arange(0.,200.1,2.)
    y0=profile(np.array([1000.]),np.array([5.]),.00005)[0]
    area=500*y0+3*y0*y0;speed=np.sqrt(9.81*area/(500+6*y0))-1000/area
    delay=cumulative_trapezoid(1/speed,x*1000,initial=0)/3600
    for path in sorted((OUT/'data').glob('A_E*.npz')):
        eid=path.stem[2:];d=loadnp(path);t=d['time_h'];q=hec200.boundary(events[eid],t,'ref')[0];h=d['boundary_increment_m']+5
        yr=profile(q,np.full_like(q,5.),.00005);yc=profile(q,h,.00005);eq=yc-yr
        dynamic=np.column_stack([np.interp(t-delay[j],t,eq[:,j],left=0.) for j in range(101)])
        for name,p in [('quasi_steady',eq),('fixed_initial_characteristic_delay',dynamic)]:
            lb=np.array([extent(v,x)[0] for v in p]);err=lb-d['lb_reference_km']
            rows.append(dict(event=eid,subset='train' if eid in split['train'] else 'validation',method=name,
                RMSE_km=float(np.sqrt(np.mean(err**2))),MAE_km=float(np.mean(abs(err))),bias_km=float(err.mean()),
                future_boundary_stage_used=True,forecast_claimed=False))
        if eid in ['E023','E006','E007','E069']:
            np.savez_compressed(OUT/'data'/f'A_hydraulic_{eid}.npz',time=t,eq=eq,dynamic=dynamic,
                lb_eq=np.array([extent(v,x)[0] for v in eq]),lb_dynamic=np.array([extent(v,x)[0] for v in dynamic]),delay_h=delay)
    frame('A_given_boundary_hydraulic_diagnostic.csv',rows)
    df=pd.DataFrame(rows);s=df.groupby(['subset','method'])[['RMSE_km','MAE_km','bias_km']].agg(['mean','median','min','max']);s.columns=['_'.join(c) for c in s.columns]
    s.to_csv(OUT/'tables/A_hydraulic_diagnostic_summary.csv',encoding='utf-8-sig')
    print('Source preservation and A oracle hydraulic diagnostics complete',flush=True)
if __name__=='__main__':run()
