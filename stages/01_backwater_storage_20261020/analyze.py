"""Read-only diagnostics of the frozen calibration set and controlled experiments."""
from pathlib import Path
import csv
import hashlib
import json
from datetime import datetime
import h5py
import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
TS = 'Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series'
ACCESSED = set()
CACHE = {}

def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))

def dump(p, obj):
    Path(p).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def csvout(name, rows):
    assert rows
    with (OUT / 'tables' / name).open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)

def npz(p):
    p = Path(p); ACCESSED.add(p)
    with np.load(p, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}

def settings(case):
    p = ROOT / 'models' / case / 'settings.json'; ACCESSED.add(p)
    return read(p)

def native(case):
    if case in CACHE:
        return CACHE[case]
    p = ROOT / 'models' / case / 'Backwater200.p01.hdf'; ACCESSED.add(p)
    with h5py.File(p, 'r') as f:
        b = f[TS]; xs = b['Cross Sections']
        x = np.array([float(v.decode().strip()) for v in xs['Cross Section Attributes'][:]['Station']])[::-1]
        t = (b['Time'][:].astype(float) - float(b['Time'][0])) * 24 - 48
        z = xs['Water Surface'][:].astype(float)[:, ::-1]
        q = xs['Flow'][:].astype(float)[:, ::-1]
        flux = xs['Flow Volume Cumulative'][:].astype(float)[:, ::-1]
    assert np.all(np.diff(t) > 0) and np.all(np.diff(x) > 0)
    y = z - x[None, :] * .00005
    area = 500 * y + 3 * y**2
    selected = np.array([int(np.where(x == k)[0][0]) for k in range(0, 200001, 10000)])
    v = cumulative_trapezoid(area, x, axis=1, initial=0)[:, selected]
    d = dict(t=t, x=x, z=z, q=q[:, selected], flux=flux[:, selected], v=v, y=y)
    CACHE[case] = d
    return d

def describe(a):
    a = np.array(a, dtype=float)
    return dict(min=float(a.min()), p10=float(np.quantile(a, .1)), median=float(np.median(a)),
                p90=float(np.quantile(a, .9)), max=float(a.max()))

def stage(e, t):
    u = np.asarray(t) - e['tb_h']; tr, th, tf = e['Tr_h'], e['Th_h'], e['Tf_h']
    return e['H0_m'] + e['deltaH_m'] * np.where(u < 0, 0,
        np.where(u < tr, u / tr, np.where(u < tr + th, 1,
        np.where(u < tr + th + tf, 1 - (u - tr - th) / tf, 0))))

def analyze_pair(var, ref, event_id, role, label, case_code='Main'):
    cfg = settings(var); e = cfg['event']
    a = npz(ROOT / 'results' / var / 'reach_timeseries.npz')
    b = npz(ROOT / 'results' / ref / 'reach_timeseries.npz')
    p = npz(ROOT / 'results' / label / 'paired_response.npz')
    m = npz(ROOT / 'results' / label / 'backwater_masks_full_sections.npz')
    av, bv = native(var), native(ref)
    t, x = p['time_since_flood_h'], p['x_from_outlet_km']
    assert np.allclose(t, a['time_since_flood_h'], rtol=0, atol=.02/3600)
    assert np.allclose(t, b['time_since_flood_h'], rtol=0, atol=.02/3600)
    assert np.allclose(t, av['t'], rtol=0, atol=.02/3600)
    assert np.allclose(t, bv['t'], rtol=0, atol=.02/3600)
    dv = p['delta_volume_m3']; dz = m['delta_stage_m'].astype(float)
    delta_flux = av['flux'] - bv['flux']
    net = (delta_flux - delta_flux[[0], :]) - (delta_flux[:, [0]] - delta_flux[0, 0])
    residual = dv - dv[[0], :] - net
    raw_dz = av['z'] - bv['z']
    # Exact trapezoidal cross-section identity, independently integrated outlet-upstream.
    area_delta = (500 + 6 * bv['y']) * raw_dz + 3 * raw_dz**2
    vcheck = cumulative_trapezoid(area_delta, av['x'], axis=1, initial=0)[:, ::5]
    assert np.max(abs(vcheck-dv)) < .01, 'Geometry and paired volume disagree'
    assert np.max(abs(residual-p['balance_residual_native_m3'])) < .01
    assert np.max(abs(av['v']-a['volume_m3'])) < .01
    shell = np.diff(dv, axis=1)
    assert np.allclose(shell.sum(axis=1), dv[:, -1], rtol=0, atol=.0001)
    peak = dv.max(axis=0); peakidx = dv.argmax(axis=0)
    # No clipping of the signed storage response; thresholds only classify water-level influence.
    for j, L in enumerate(x[1:], 1):
        k = int(peakidx[j]); sj = int(shell[:, j-1].argmax())
        METRICS.append(dict(event_id=event_id, role=role, case=var, case_code=case_code, L_km=int(L),
            peak_delta_volume_m3=float(peak[j]), min_delta_volume_m3=float(dv[:, j].min()),
            time_peak_h=float(t[k]), time_peak_over_D=float(t[k]/e['D_h']),
            peak_minus_previous_length_peak_m3=float(peak[j]-peak[j-1]),
            shell_peak_m3=float(shell[sj, j-1]), time_shell_peak_h=float(t[sj]),
            shell_at_200km_peak_m3=float(shell[int(peakidx[-1]), j-1]),
            peak_fraction_of_200km=float(peak[j]/peak[-1]),
            storage_time_integral_m3_h=float(trapezoid(dv[:, j], t)),
            max_native_balance_residual_m3=float(np.max(abs(residual[:, j])))))
    criterion_stats = []
    for threshold in [.05, .10, .20]:
        pos = raw_dz[:, 1:] > threshold; neg = raw_dz[:, 1:] < -threshold
        interior_x = av['x'][1:] / 1000
        far = np.max(np.where(pos, interior_x[None, :], 0), axis=1)
        contiguous = np.sum(np.cumprod(pos.astype(np.int8), axis=1), axis=1) * 2
        segments = np.sum(np.diff(np.pad(pos.astype(np.int8), ((0,0),(1,1))), axis=1) == 1, axis=1)
        assert np.allclose(far, m[f'extent_{int(threshold*100):02d}cm_km'])
        dt = np.diff(t)
        touch_h = float(np.sum(dt*(pos[:-1,-1].astype(float)+pos[1:,-1].astype(float))/2))
        criterion_stats.append(dict(threshold_m=threshold, max_farthest_km=float(far.max()),
            max_connected_km=float(contiguous.max()), max_segments=int(segments.max()),
            domain_boundary_touched=bool(pos[:,-1].any()), boundary_touch_h=touch_h,
            negative_threshold_samples=int(neg.sum())))
        for k, tt in enumerate(t):
            EXTENT.append(dict(event_id=event_id, role=role, case_code=case_code, threshold_m=threshold,
                time_h=float(tt), time_over_D=float(tt/e['D_h']), farthest_km=float(far[k]),
                connected_km=float(contiguous[k]), segments=int(segments[k]),
                boundary_touched=int(pos[k,-1]), negative_sections=int(neg[k].sum())))
    l90 = float(x[np.flatnonzero(peak >= .9*peak[-1])[0]])
    l95 = float(x[np.flatnonzero(peak >= .95*peak[-1])[0]])
    hreturn = float(e['Hend_h'])
    after = t >= hreturn
    returned_v = float(np.interp(hreturn, t, dv[:,-1]))
    post_t = t[after]; post_dv = dv[after,-1]
    active_after = post_t[post_dv > .01*peak[-1]]
    retention = 0.0 if not len(active_after) else float(active_after[-1]-hreturn)
    event_row = dict(event_id=event_id, role=role, case=var, case_code=case_code,
        D_h=e['D_h'], Qp_m3s=e['Qp_m3s'], deltaH_m=e['deltaH_m'], tb_ratio=e['tb_ratio'],
        Th_ratio=e['Th_ratio'], Tr_ratio=e['Tr_ratio'], Tf_ratio=e['Tf_ratio'],
        peak_200km_m3=float(peak[-1]), peak_time_h=float(t[int(peakidx[-1])]),
        L90_km=l90, L95_km=l95, final10km_peak_difference_m3=float(peak[-1]-peak[-2]),
        monotone_peak_length=bool(np.all(np.diff(peak) >= -.01)),
        min_delta_volume_m3=float(dv.min()), max_native_balance_residual_m3=float(abs(residual).max()),
        max_native_residual_over_peak200=float(abs(residual).max()/peak[-1]),
        geometry_crosscheck_error_m3=float(abs(vcheck-dv).max()),
        max_interior_delta_stage_m=float(raw_dz[:,1:].max()), min_delta_stage_m=float(raw_dz.min()),
        boundary_return_h=hreturn, storage_at_boundary_return_m3=returned_v,
        fraction_peak_at_boundary_return=returned_v/peak[-1],
        last_sample_above_one_pct_peak_after_return_h=retention,
        min_outlet_flow_m3s=float(av['q'][:,0].min()),
        max_paired_inlet_flow_difference_100km_m3s=float(abs(av['q'][:,10]-bv['q'][:,10]).max()),
        native_paired_inlet_volume_at_peak100_m3=float((delta_flux[int(peakidx[10]),10]-delta_flux[0,10])),
        native_paired_outlet_volume_at_peak100_m3=float((delta_flux[int(peakidx[10]),0]-delta_flux[0,0])))
    EVENTS.append(event_row)
    for row in criterion_stats:
        THRESHOLD.append(dict(event_id=event_id, role=role, case_code=case_code, **row))
    for level in [.25,.5,.75]:
        rise=e['tb_h']+level*e['Tr_h']
        fall=e['tb_h']+e['Tr_h']+e['Th_h']+(1-level)*e['Tf_h']
        for L in [10,100,200]:
            j=int(L/10); vr=float(np.interp(rise,t,dv[:,j])); vf=float(np.interp(fall,t,dv[:,j]))
            MATCHED.append(dict(event_id=event_id, role=role, case_code=case_code,L_km=L,
                equal_boundary_increment_m=level*e['deltaH_m'], rise_time_h=rise,fall_time_h=fall,
                rise_storage_m3=vr,fall_storage_m3=vf,fall_minus_rise_m3=vf-vr,
                inlet_flow_rise_m3s=float(np.interp(rise,t,av['q'][:,j])),
                inlet_flow_fall_m3s=float(np.interp(fall,t,av['q'][:,j]))))
    np.savez_compressed(OUT/'derived_cases'/f'{var}.npz', time_h=t,x_km=x,full_x_km=av['x']/1000,
        delta_volume_m3=dv,delta_stage_full_m=raw_dz, shell_m3=shell,
        native_balance_residual_m3=residual,stage_increment_boundary_m=stage(e,t)-e['H0_m'],
        outlet_flow_m3s=av['q'][:,0],reference_outlet_flow_m3s=bv['q'][:,0],
        delta_inlet_flow_m3s=av['q']-bv['q'])
    return dv, t, peak

def main():
    global METRICS, EVENTS, EXTENT, THRESHOLD, MATCHED
    METRICS=[]; EVENTS=[]; EXTENT=[]; THRESHOLD=[]; MATCHED=[]
    for d in ['tables','derived_cases','figures','provenance']:
        (OUT/d).mkdir(exist_ok=True)
    plan=read(ROOT/'config/formal_plan.json'); checkpoint=read(ROOT/'logs/formal_checkpoint.json')
    expected={c['case'] for c in plan['formal_cases']}
    assert len(expected)==232 and expected==set(checkpoint['formal_completed'])
    assert not checkpoint['formal_failures']
    audit=[]
    for c in plan['formal_cases']:
        q=read(ROOT/'results'/c['case']/'quality.json')
        checks=['complete_duration','output_interval_matches_plan','initialization_pass','native_balance_pass','tail_recovered','basic_quality_pass','htab_coverage']
        assert all(q[k] is True for k in checks)
        assert q['nonfinite_count']==0 and q['bank_exceedance_count']==0
        run=read(ROOT/'models'/c['case']/'run.json')
        assert run['status']=='computed' and run['exit_code']==0
        for path in [ROOT/'models'/c['case']/'Backwater200.p01.hdf',ROOT/'results'/c['case']/'reach_timeseries.npz']:
            assert path.stat().st_size > 0
        audit.append(dict(case=c['case'],role=c['role'],event_id=c['event_id'],split=c['split'],quality_pass=True))
    csvout('all232_completeness_audit.csv',audit)
    main_rows=[c for c in plan['formal_cases'] if c['role']=='main']
    byevent={}
    for c in main_rows:
        byevent.setdefault(c['event_id'],{})[c['kind']]=c
    calibration=sorted(eid for eid,g in byevent.items() if g['var']['split']=='calibration')
    validation=sorted(set(byevent)-set(calibration))
    assert len(calibration)==70 and len(validation)==30
    assert set(plan['groups'])<=set(calibration)
    dump(OUT/'provenance/holdout_manifest.json',dict(status='withheld_from_stage01_mechanism_analysis',
        calibration_event_ids=calibration,validation_event_ids=validation,
        validation_results_used=False,all232_quality_metadata_checked=True))
    for eid in calibration:
        g=byevent[eid]
        analyze_pair(g['var']['case'],g['ref']['case'],eid,'calibration_main','PAIR_FORMAL_'+eid)
    contrast=[]; factorial=[]
    signs=np.array([[-1,-1,-1],[-1,-1,1],[-1,1,-1],[-1,1,1],[1,-1,-1],[1,-1,1],[1,1,-1],[1,1,1]])
    for eid,group in sorted(plan['groups'].items()):
        peaks=[]; arrays=[]; times=None
        for i,case in enumerate(group['addon']):
            dv,t,p=analyze_pair(case,group['main']['ref'],eid,'controlled','PAIR_'+case,f'C{i+1:02d}')
            arrays.append(dv); peaks.append(p); times=t
        peaks=np.array(peaks); arrays=np.array(arrays)
        e=settings(group['addon'][0])['event']
        effects={}
        for factor,lohis in {'A':[(0,4),(1,5),(2,6),(3,7)],'B':[(0,2),(1,3),(4,6),(5,7)],'C':[(0,1),(2,3),(4,5),(6,7)]}.items():
            for lo,hi in lohis:
                d=arrays[hi]-arrays[lo]
                for j,L in enumerate(range(10,201,10),1):
                    contrast.append(dict(event_id=eid,factor=factor,low=f'C{lo+1:02d}',high=f'C{hi+1:02d}',L_km=L,
                        difference_individual_peaks_m3=float(peaks[hi,j]-peaks[lo,j]),
                        max_instantaneous_difference_m3=float(d[:,j].max()),
                        min_instantaneous_difference_m3=float(d[:,j].min()),
                        difference_storage_time_integral_m3_h=float(trapezoid(d[:,j],times))))
        for name,cols in [('A',[0]),('B',[1]),('C',[2]),('AB',[0,1]),('AC',[0,2]),('BC',[1,2]),('ABC',[0,1,2])]:
            s=np.prod(signs[:,cols],axis=1)
            effect=np.tensordot(s,arrays,axes=(0,0))/4
            effects['effect_'+name+'_m3']=effect
            pe=s@peaks/4
            for j,L in enumerate(range(10,201,10),1):
                factorial.append(dict(event_id=eid,factor=name,L_km=L,effect_on_individual_peaks_m3=float(pe[j]),
                    max_instantaneous_effect_m3=float(effect[:,j].max()),min_instantaneous_effect_m3=float(effect[:,j].min())))
        np.savez_compressed(OUT/'derived_cases'/f'factorial_{eid}.npz',time_h=times,time_over_D=times/e['D_h'],**effects)
    for name,rows in [('storage_by_length.csv',METRICS),('event_summary.csv',EVENTS),('influence_timeseries.csv',EXTENT),
                      ('threshold_summary.csv',THRESHOLD),('equal_boundary_stage_comparison.csv',MATCHED),
                      ('controlled_contrasts.csv',contrast),('factorial_effects.csv',factorial)]:
        csvout(name,rows)
    cal=[e for e in EVENTS if e['role']=='calibration_main']; ctrl=[e for e in EVENTS if e['role']=='controlled']
    threshold_summary={}
    for th in [.05,.1,.2]:
        rows=[e for e in THRESHOLD if e['role']=='calibration_main' and e['threshold_m']==th]
        threshold_summary[str(th)]=dict(max_extent_km=describe([a['max_farthest_km'] for a in rows]),
            boundary_censored_events=sum(a['domain_boundary_touched'] for a in rows),
            disjoint_events=sum(a['max_segments']>1 for a in rows),
            negative_response_events=sum(a['negative_threshold_samples']>0 for a in rows))
    effects_summary={}
    for f in ['A','B','C']:
        rows=[v for v in contrast if v['factor']==f and v['L_km']==200]
        effects_summary[f]=dict(n_flood_backgrounds=4,n_conditional_contrasts=len(rows),
            peak_difference_m3=describe([v['difference_individual_peaks_m3'] for v in rows]),
            positive_peak_differences=sum(v['difference_individual_peaks_m3']>0 for v in rows),
            integral_difference_m3_h=describe([v['difference_storage_time_integral_m3_h'] for v in rows]),
            by_background={eid:float(np.mean([v['difference_individual_peaks_m3'] for v in rows if v['event_id']==eid])) for eid in sorted(plan['groups'])})
    summary=dict(created=datetime.now().astimezone().isoformat(),stage='01_backwater_storage',
        formal_runs=232,calibration_pairs=70,withheld_validation_pairs=30,controlled_runs=32,
        main_summary={k:describe([e[k] for e in cal]) for k in ['peak_200km_m3','L90_km','L95_km',
            'fraction_peak_at_boundary_return','last_sample_above_one_pct_peak_after_return_h',
            'max_native_residual_over_peak200','geometry_crosscheck_error_m3','final10km_peak_difference_m3']},
        main_nonmonotone_peak_profiles=sum(not e['monotone_peak_length'] for e in cal),
        thresholds=threshold_summary,controlled_effects_at200km=effects_summary,
        minimum_signed_storage_m3=min(e['min_delta_volume_m3'] for e in EVENTS),
        max_geometry_crosscheck_error_m3=max(e['geometry_crosscheck_error_m3'] for e in EVENTS),
        minimum_signed_stage_m=min(e['min_delta_stage_m'] for e in EVENTS),
        equal_stage_200km_50pct_controlled_fall_minus_rise_m3=describe([e['fall_minus_rise_m3'] for e in MATCHED
            if e['role']=='controlled' and e['L_km']==200 and abs(e['equal_boundary_increment_m']-(3 if int(e['case_code'][1:])>=5 else 1.5))<1e-6]),
        caveats=['Finite domain; 200 km threshold hits are right censored.',
            'Same boundary stage does not imply same inlet flow or hydraulic history.',
            'Hold contrast includes delayed start of drawdown.',
            'Nested reaches and factor contexts are not independent flood replicates.',
            'No new RWEC formula fitted and no validation-target features used.'])
    dump(OUT/'summary.json',summary)
    manifest=[]
    for p in sorted(ACCESSED):
        # Metadata-only fingerprints for the native HDFs; portable derived/NPZ sources are SHA256 checked.
        manifest.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,
            sha256=sha(p) if p.suffix!='.hdf' else None,mtime_ns=p.stat().st_mtime_ns))
    dump(OUT/'provenance/input_manifest.json',dict(plan_sha256=sha(ROOT/'config/formal_plan.json'),inputs=manifest,
        hdf_open_mode='read-only',excluded_validation_events=validation,mechanism_analysis_case_count=102))
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
