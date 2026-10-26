"""200-km HEC experiment with an eighth LHS parameter Th/D in[0.10,0.30]."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import shutil
import uuid
import sys
import time
from datetime import datetime, timedelta
import h5py
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import qmc

ROOT = Path(__file__).resolve().parents[1]
BASELINE_ROOT = ROOT.parents[1]
DESIGN_ID = 'forward_dynamic_C_v1'
NAME = 'Backwater200'
TS = 'Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series'
START = datetime(2026, 10, 26)
ORIGIN_H = 48.0
SEED = 20261019

def now():
    return datetime.now().astimezone().isoformat(timespec='seconds')

def dump(p, v):
    p = Path(p)
    tmp = p.with_name(p.name + '.' + uuid.uuid4().hex + '.tmp')
    tmp.write_text(json.dumps(v, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    try:
        for attempt in range(30):
            try:
                tmp.replace(p)
                break
            except PermissionError:
                if attempt==29: raise
                time.sleep(.05 + .01*attempt)
    finally:
        tmp.unlink(missing_ok=True)

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def text(p, v):
    Path(p).write_bytes(v.replace('\r\n', '\n').replace('\n', '\r\n').encode('ascii'))

def csvwrite(p, header, rows):
    with Path(p).open('w', newline='', encoding='utf-8-sig') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)

def event(qp=4250., duration=132., r=.275, dh=4.5, tb=-.05, tr=.45, tf=.6, event_id='PILOT_MID', th=.20):
    a = math.log(.001)/(math.log(1/r)+1-1/r)
    return dict(event_id=event_id, Qb_m3s=1000., Qp_m3s=qp, D_h=duration,
                r=r, Tp_h=r*duration, gamma_a=a, H0_m=5., deltaH_m=dh,
                Hmax_m=5.+dh, tb_ratio=tb, Tr_ratio=tr, Tf_ratio=tf,
                tb_h=tb*duration, Tr_h=tr*duration, Th_ratio=th,Th_h=th*duration, Tf_h=tf*duration,
                Hpeak_h=(tb+tr)*duration, Hfall_start_h=(tb+tr+th)*duration,Hend_h=(tb+tr+th+tf)*duration)

def prepare_design():
    dest = ROOT/'config/events.json'
    if dest.exists():
        saved = json.loads(dest.read_text(encoding='utf-8'))
        assert saved['seed'] == SEED and len(saved['events']) == 100 and saved['design_id']==DESIGN_ID
        return saved['events']
    original=json.loads((BASELINE_ROOT/'config/events.json').read_text(encoding='utf-8'))
    old_unit=np.load(BASELINE_ROOT/'config/lhs_unit_cube.npy')
    hold_unit=qmc.LatinHypercube(d=1,seed=SEED+2).random(100)[:,0]
    u=np.column_stack([old_unit[:,:6],hold_unit,old_unit[:,6]])
    low = [2500.,96.,.20,3.,-.25,.30,.10,.40]
    high = [6000.,168.,.35,6.,.25,.60,.30,.80]
    events = []
    for i, old in enumerate(original['events']):
        e=dict(old);ratio=float(.10+.20*hold_unit[i])
        e.update(Th_ratio=ratio,Th_h=ratio*e['D_h'],Hfall_start_h=e['Hpeak_h']+ratio*e['D_h'],
                 Hend_h=e['Hend_h']+ratio*e['D_h'])
        events.append(e)
    dump(dest, dict(design_id=DESIGN_ID,seed=SEED,hold_seed=SEED+2,split_seed=SEED+1,created=now(),dimensions=8,
                    parent_events_sha256=sha(BASELINE_ROOT/'config/events.json'),
                    first_seven_parameters_and_split_unchanged=True,
                    names=['Qp','D','r','deltaH','tb/D','Tr/D','Th/D','Tf/D'],
                    low=low,high=high,events=events))
    csvwrite(ROOT/'config/event_parameters.csv',list(events[0]),[list(e.values()) for e in events])
    np.save(ROOT/'config/lhs_unit_cube.npy',u)
    assert all(len(np.unique((u[:,j]*100).astype(int)))==100 for j in range(8))
    dump(ROOT/'config/design_audit.json',dict(n_events=100,n_calibration=70,n_validation=30,
        dimensions=8,latin_strata_all_dimensions=True,Th_ratio_range=[.1,.3],seed=SEED,hold_seed=SEED+2,split_seed=SEED+1,
        events_sha256=sha(dest),created=now()))
    return events

def boundary(e, t_h, kind):
    t=np.asarray(t_h,dtype=float)
    q=np.full(t.shape,e['Qb_m3s'],dtype=float)
    mask=t>0; v=t[mask]/e['Tp_h']
    q[mask]+=(e['Qp_m3s']-e['Qb_m3s'])*np.exp(e['gamma_a']*(np.log(v)+1-v))
    if e.get('secondary_ratio',0)>0:
        tt=t-e['secondary_delay_h'];m=tt>0;v2=tt[m]/e['Tp_h']
        q[m]+=e['secondary_ratio']*(e['Qp_m3s']-e['Qb_m3s'])*np.exp(e['gamma_a']*(np.log(v2)+1-v2))
    h=np.full(t.shape,5.,dtype=float)
    if kind=='var':
        tb=e['tb_h']; tr=e['Tr_h']; th=e['Th_h']; tf=e['Tf_h']; dh=e['deltaH_m']
        rising=(t>=tb)&(t<tb+tr)
        holding=(t>=tb+tr)&(t<tb+tr+th)
        falling=(t>=tb+tr+th)&(t<tb+tr+th+tf)
        h[rising]+=dh*(1-np.cos(np.pi*(t[rising]-tb)/tr))/2
        h[holding]+=dh
        h[falling]+=dh*(1+np.cos(np.pi*(t[falling]-tb-tr-th)/tf))/2
    elif kind not in ['ref','base']:
        raise ValueError(kind)
    if kind=='base': q[:]=1000.
    return q,h

def update_keys(s,updates):
    lines=[]; seen=set()
    for line in s.splitlines():
        if line.startswith('Stage Flow Hydrograph='): continue
        key=line.partition('=')[0]
        if key in updates: lines.append(f'{key}={updates[key]}');seen.add(key)
        else: lines.append(line)
    lines += [f'{k}={v}' for k,v in updates.items() if k not in seen]
    return '\n'.join(lines)+'\n'

def build(case,e,kind,dx=2000,dt=60,end_h=None,htab_dh=.05,htab_rows=481,slope=.000075):
    folder=ROOT/'models'/case
    folder.mkdir(parents=True,exist_ok=True)
    if (folder/'settings.json').exists():
        cfg=json.loads((folder/'settings.json').read_text(encoding='utf-8'))
        assert cfg['event']==e and cfg['kind']==kind and cfg['dx_m']==dx and cfg['dt_s']==dt
        assert cfg['bed_slope']==slope
        assert cfg['htab_dh_m']==htab_dh and cfg['htab_rows']==htab_rows
        assert cfg['warmup_steps']==round(86400/dt), 'Use a new case for a full24h warmup at the computation interval'
        if end_h is not None: assert cfg['end_since_flood_h']==end_h
        for name,digest in cfg['input_hashes'].items():
            assert sha(folder/name)==digest, f'Changed input {folder/name}'
        return folder
    assert not (folder/f'{NAME}.p01.hdf').exists(), 'Preserve existing outputs'
    if end_h is None:
        end_h=math.ceil((max(e['D_h'],e['Hend_h'])+168)/24)*24
    duration_h=ORIGIN_H+end_h
    assert duration_h>0 and 200000%dx==0 and dt>0
    # At depth<=12 m this is exactly the requested 500 m / 3H:1V trapezoid.
    # Continue the same side slopes to 24 m solely for numeric coverage.
    # Any depth>12 m is recorded as bank exceedance, never silently accepted.
    g=[f'Geom Title={case}','Program Version=6.31','Viewing Rectangle=0,200000,1000,-1000','',
       'River Reach=Ideal200        ,Main            ','Reach XY= 2 ',
       f'{0:16g}{0:16g}{200000:16g}{0:16g}','Reverse River Text=0','']
    rsvals=np.arange(200000,-1,-dx,dtype=int)
    for rs in rsvals:
        z=rs*slope; length=dx if rs>0 else 0
        pts=[(0,z+24),(36,z+12),(72,z),(572,z),(608,z+12),(644,z+24)]
        g += [f'Type RM Length L Ch R = 1 ,{rs},{length},{length},{length}', '#Sta/Elev= 6 ']
        for i in range(0,6,5):
            g.append(''.join(f'{x:8.2f}{y:8.4f}' for x,y in pts[i:i+5]))
        g += ['#Mann= 3 ,0,0',''.join(f'{x:8g}{.025:8g}{0:8g}' for x in [0,36,608]),
              'Bank Sta=36,608','XS Rating Curve=0,0',
              f'XS HTab Starting El and Incr={z+.01:.4f},{htab_dh:g},{htab_rows}',
              'XS HTab Horizontal Distribution=5,5,5','Exp/Cntr=0,0','']
    g += ['Use User Specified Reach Order=0','GIS Units=METERS','']
    text(folder/f'{NAME}.g01','\n'.join(g))
    cache=ROOT/'config/geometry_cache'/f'{dx}m'
    # Disabled until HEC's interactive cache-validation prompt is eliminated.
    reuse=cache.exists() and (ROOT/'config/ENABLE_GEOMETRY_CACHE').exists()
    if reuse:
        cached=(cache/f'{NAME}.g01').read_text()
        generated=(folder/f'{NAME}.g01').read_text()
        normalize=lambda value:'\n'.join(line for line in value.splitlines() if not line.startswith('Geom Title='))
        assert normalize(cached)==normalize(generated),'Geometry cache does not match requested hydraulic geometry'
        audit=json.loads((cache/'manifest.json').read_text(encoding='utf-8'))
        for name,digest in audit['hashes'].items():
            assert sha(cache/name)==digest,'Geometry cache changed'
            shutil.copyfile(cache/name,folder/name)
    t_h=np.arange(round(duration_h*60)+1)/60-ORIGIN_H
    q,h=boundary(e,t_h,kind)
    csvwrite(folder/'boundary.csv',['t_since_flood_h','elapsed_h','inflow_m3s','stage_m'],
             zip(t_h,t_h+ORIGIN_H,q,h))
    u=[f'Flow Title={case}','Program Version=6.31','Use Restart=0',
       'Initial Flow Loc=Ideal200        ,Main            ,200000,1000',
       'Boundary Location=Ideal200        ,Main            ,200000,        ,                ,                ,                ,                                ,                                ',
       'Interval=1MIN',f'Flow Hydrograph= {len(q)} ']
    u += [''.join(f'{v:8.3f}' for v in q[i:i+10]) for i in range(0,len(q),10)]
    trailer=['Stage Hydrograph TW Check=0','DSS Path=','Use DSS=False','Use Fixed Start Time=False',
             'Fixed Start Date/Time=,','Is Critical Boundary=False','Critical Boundary Flow=']
    u += trailer+['Boundary Location=Ideal200        ,Main            ,0,        ,                ,                ,                ,                                ,                                ',
                  'Interval=1MIN',f'Stage Hydrograph= {len(h)} ']
    u += [''.join(f'{v:8.5f}' for v in h[i:i+10]) for i in range(0,len(h),10)]
    u += trailer+['Met Point Raster Parameters=,,,,','Precipitation Mode=Disable','Wind Mode=No Wind Forces','']
    text(folder/f'{NAME}.u01','\n'.join(u))
    end=START+timedelta(hours=duration_h)
    template=(ROOT/'config/plan_template_631.p01').read_text()
    p=update_keys(template,{'Plan Title':case,'Short Identifier':case[:16],
       'Simulation Date':f'{START:%d%b%Y},{START:%H%M},{end:%d%b%Y},{end:%H%M}'.upper(),
       'Computation Interval':('1MIN' if dt==60 else f'{dt}SEC'),'Output Interval':'1HOUR','Instantaneous Interval':'1HOUR','Mapping Interval':'1HOUR',
       'UNET Theta':'1','UNET Theta Warmup':'1','UNET ZTol':'0.005','UNET ZSATol':'0.005',
       'UNET MxIter':'30','UNET MaxInSteps':str(round(86400/dt)),'UNET DtIC':'0','UNET Use Existing IB Tables':'0',
       'UNET D1 Cores':'1','HDF Write Warmup':'0','Write IC File at Sim End':'0',
       'Run HTab':'0' if reuse else '-1',
       'Run PostProcess':'0','Run RASMapper':'0','HDF Compression':'1'})
    text(folder/f'{NAME}.p01',p)
    prj=(f'Proj Title={case}\nCurrent Plan=p01\nDefault Exp/Contr=0,0\nSI Units\n'
         'Geom File=g01\nUnsteady File=u01\nPlan File=p01\nDSS File=dss\n'
         'BEGIN DESCRIPTION:\n200 km trapezoidal channel. Gamma inflow; paired stage. High-stage hold=0.10D-0.30D.\nEND DESCRIPTION:\n')
    text(folder/f'{NAME}.prj',prj)
    cfg=dict(design_id=DESIGN_ID,case=case,event=e,kind=kind,dx_m=dx,dt_s=dt,n_sections=len(rsvals),
             length_m=200000.,width_m=500.,side_slope=3.,bank_depth_m=12.,geometry_top_depth_m=24.,
             bed_slope=slope,manning_n=.025,output_dt_s=3600,input_dt_s=60,
             start_date=START.isoformat(),flood_origin_h=ORIGIN_H,start_since_flood_h=-ORIGIN_H,
             end_since_flood_h=end_h,duration_h=duration_h,htab_dh_m=htab_dh,htab_rows=htab_rows,
             theta=1.,warmup_steps=round(86400/dt),warmup_dt_s=dt,warmup_interval_rule='default computation interval; UNET DtIC=0',created=now(),
             geometry_cache_reused=reuse,
             geometry_cache_manifest=str(cache/'manifest.json') if reuse else None,
             input_hashes={f'{NAME}{suffix}':sha(folder/f'{NAME}{suffix}') for suffix in ['.g01','.p01','.prj','.u01']})
    dump(folder/'settings.json',cfg)
    return folder

def worker(folder):
    import win32com.client
    start=time.monotonic(); controller=None
    cfg=json.loads((folder/'settings.json').read_text(encoding='utf-8'))
    for name,digest in cfg['input_hashes'].items(): assert sha(folder/name)==digest
    dump(folder/'run.json',dict(status='running',worker_pid=os.getpid(),started=now()))
    try:
        controller=win32com.client.DispatchEx('RAS631.HECRASController')
        version=controller.HECRASVersion()
        controller.Project_Open(str(folder/f'{NAME}.prj'))
        controller.Compute_HideComputationWindow()
        ret=controller.Compute_CurrentPlan(0,None,True)
        dump(folder/'run.json',dict(status='computed' if ret[0] else 'failed',returned=ret,
             version=version,elapsed_s=round(time.monotonic()-start,2),finished=now(),
             hdf_exists=(folder/f'{NAME}.p01.hdf').exists(),exit_code=0 if ret[0] else 1))
        if not ret[0]: raise RuntimeError('HEC computation failed')
    except BaseException as exc:
        current=json.loads((folder/'run.json').read_text(encoding='utf-8'))
        current.update(status='failed',error=repr(exc),elapsed_s=time.monotonic()-start,finished=now(),exit_code=1)
        dump(folder/'run.json',current)
        raise
    finally:
        if controller is not None:
            controller.Project_Save()
            controller.QuitRas()

def run(folder):
    result=folder/'run.json'
    if result.exists():
        saved=json.loads(result.read_text(encoding='utf-8'))
        if saved.get('exit_code')==0 and (folder/f'{NAME}.p01.hdf').exists(): return saved
        if saved.get('status')=='running':
            # A restarted coordinator can adopt its already-running worker.
            # Validate identity; never create a second RAS run for the same case.
            import psutil
            pid=saved['worker_pid']; started=time.monotonic()
            proc=psutil.Process(pid);cmd=proc.cmdline()
            assert 'worker' in cmd and folder.name in cmd and any('hec200.py' in s for s in cmd)
            while True:
                current=json.loads(result.read_text(encoding='utf-8'))
                if current.get('exit_code')==0 and (folder/f'{NAME}.p01.hdf').exists(): return current
                if current.get('status')!='running': raise RuntimeError(f'Adopted run failed: {folder}')
                if not proc.is_running(): raise RuntimeError(f'Worker exited without finished record: {folder}')
                if time.monotonic()-started>1800: raise RuntimeError(f'Adopted worker wait exceeded30min: {folder}')
                time.sleep(5)
        raise RuntimeError(f'Prior running/failed attempt preserved: {folder}')
    with (folder/'compute_stdout.log').open('w') as out,(folder/'compute_stderr.log').open('w') as err:
        p=subprocess.Popen([sys.executable,str(Path(__file__)), 'worker',folder.name],
                           cwd=folder,stdout=out,stderr=err,creationflags=subprocess.CREATE_NO_WINDOW)
        dump(ROOT/'logs/active_runs'/f'{folder.name}.json',dict(case=folder.name,pid=p.pid,status='running',started=now()))
        try:
            code=p.wait(timeout=1800)
        except subprocess.TimeoutExpired:
            # Only the worker tree created by this run is terminated.
            subprocess.run(['taskkill','/PID',str(p.pid),'/T','/F'],capture_output=True)
            saved=json.loads(result.read_text(encoding='utf-8')) if result.exists() else {}
            saved.update(status='timeout_review_required',exit_code=1,finished=now())
            dump(result,saved)
            raise RuntimeError(f'30 min worker timeout; inspect preserved attempt: {folder}')
    if code: raise RuntimeError(f'Worker failed: {folder}')
    return json.loads(result.read_text(encoding='utf-8'))

def extract(folder):
    cfg=json.loads((folder/'settings.json').read_text(encoding='utf-8'))
    out=ROOT/'results'/folder.name;out.mkdir(exist_ok=True)
    with h5py.File(folder/f'{NAME}.p01.hdf','r') as f:
        b=f[TS]; xs=b['Cross Sections']; raw=b['Time'][:].astype(float)
        sec=(raw-raw[0])*86400.; t_h=sec/3600-cfg['flood_origin_h']
        attrs=xs['Cross Section Attributes'][:]
        stations=np.array([float(v.decode().strip().replace('*','')) for v in attrs['Station']])
        assert len(stations)==cfg['n_sections'] and np.all(np.diff(stations)==-cfg['dx_m'])
        bed=stations*cfg['bed_slope']
        q=xs['Flow'][:].astype(float);z=xs['Water Surface'][:].astype(float)
        vel=xs['Velocity Channel'][:].astype(float); depth=z-bed
        area=500*depth+3*depth**2
        xp=200000-stations
        cumarea=cumulative_trapezoid(area,xp,axis=1,initial=0)
        select=np.array([int(np.where(stations==x)[0][0]) for x in range(0,200001,10000)])
        volume=cumarea[:,-1,None]-cumarea[:,select]
        qsel=q[:,select];zsel=z[:,select];ysel=depth[:,select];asel=area[:,select];vsel=vel[:,select]
        slope_down=-np.gradient(z,xp,axis=1)[:,select]
        net5=cumulative_trapezoid(qsel-qsel[:,[0]],sec,axis=0,initial=0)
        native=None
        if 'Flow Volume Cumulative' in xs:
            flux=xs['Flow Volume Cumulative'][:].astype(float)[:,select]
            native=(flux-flux[[0],:])-(flux[:,[0]]-flux[0,0])
        residual=volume-volume[[0],:]-net5
        native_residual=volume-volume[[0],:]-native if native is not None else None
        # Depth>12 is bank exceedance even though geometry is extended to24 m.
        bad=int((~np.isfinite(q)).sum()+(~np.isfinite(z)).sum()+(~np.isfinite(vel)).sum())
        qexpected,hexpected=boundary(cfg['event'],t_h,cfg['kind'])
        t_pre=(t_h>=min(0,cfg['event']['tb_h'])-6)&(t_h<min(0,cfg['event']['tb_h'])-.1)
        tail=t_h>=t_h[-1]-24
        peak_change=np.max(np.abs(volume-volume[[0],:]),axis=0)
        rel=np.max(np.abs(residual),axis=0)/np.maximum(peak_change,1.)
        relnative=(np.max(np.abs(native_residual),axis=0)/np.maximum(peak_change,1.)).tolist() if native is not None else None
        # Verify actual HEC geometry, not merely source file generation.
        geom=f['Geometry/Cross Sections']
        ht=geom['Property Tables']; info=ht['XSEC Info'][:]; ev=ht['XSEC Value']
        names=[v[0].decode().strip() for v in ev.attrs['Variables']]
        ei=names.index('Elevation')
        bottom=np.array([ev[int(s),ei] for s,n,_ in info])
        top=np.array([ev[int(s+n-1),ei] for s,n,_ in info])
        coverage=bool(np.all(z>=bottom) and np.all(z<=top))
        fluxin=float(np.trapz(q[:,0],sec))
        initialflow=float(np.max(abs(q[0]-1000)))
        preflow=float(np.max(abs(q[t_pre]-1000))) if t_pre.any() else None
        metrics=dict(case=folder.name,checked=now(),n_times=len(sec),n_sections=len(stations),
            duration_h=float(sec[-1]/3600),planned_duration_h=cfg['duration_h'],
            complete_duration=bool(abs(sec[-1]/3600-cfg['duration_h'])<.001),
            output_interval_matches_plan=bool(np.max(abs(np.diff(sec)-cfg['output_dt_s']))<.02),
            nonfinite_count=bad,min_flow_m3s=float(q.min()),max_flow_m3s=float(q.max()),
            min_depth_m=float(depth.min()),max_depth_m=float(depth.max()),
            bank_exceedance_count=int((depth>12).sum()),flow_reversal_count=int((q<-.01).sum()),
            max_upstream_boundary_error_m3s=float(np.max(abs(q[:,0]-qexpected))),
            max_downstream_stage_error_m=float(np.max(abs(z[:,-1]-hexpected))),
            initial_flow_error_m3s=initialflow,pre_boundary_flow_error_m3s=preflow,
            final_flow_error_m3s=float(np.max(abs(q[-1]-1000))),
            last24h_flow_error_m3s=float(np.max(abs(q[tail]-1000))),
            final_stage_departure_from_initial_m=float(np.max(abs(z[-1]-z[0]))),
            final_storage_departure_fraction=float(abs(volume[-1,-1]-volume[0,-1])/max(peak_change[-1],1.)),
            maximum_balance_residual_m3=float(np.max(abs(residual))),
            final_balance_residual_pct_total_inflow=float(100*residual[-1,-1]/fluxin),
            max_balance_residual_fraction_storage_by_length=rel.tolist(),
            native_max_balance_residual_fraction_storage_by_length=relnative,
            native_max_balance_residual_m3=float(np.max(abs(native_residual))) if native is not None else None,
            actual_htab_rows_min=int(info[:,1].min()),actual_htab_rows_max=int(info[:,1].max()),
            htab_coverage=coverage,htab_top_min_margin_m=float(np.min(top-z.max(axis=0))),
            geometry_area_definition='A=500*y+3*y^2; same side slopes continued above bank12 for coverage',
            discharge_and_area_sign='x grows upstream; positive Q is downstream; volume is positive',
            native_result_variables=list(xs.keys()))
        baseline=cfg['kind']=='base'
        metrics['initialization_pass']=bool(initialflow<1 and (preflow is None or preflow<1))
        metrics['native_balance_pass']=bool(relnative is not None and max(relnative[1:])<.01)
        metrics['tail_recovered']=bool(metrics['last24h_flow_error_m3s']<1 and metrics['final_stage_departure_from_initial_m']<.005)
        metrics['basic_quality_pass']=bool(metrics['complete_duration'] and metrics['output_interval_matches_plan'] and
            bad==0 and metrics['min_depth_m']>0 and coverage and metrics['max_upstream_boundary_error_m3s']<.02 and
            metrics['max_downstream_stage_error_m']<.001 and metrics['bank_exceedance_count']==0 and
            metrics['native_balance_pass'])
        if baseline: metrics['baseflow_pass']=bool(metrics['basic_quality_pass'] and metrics['last24h_flow_error_m3s']<1)
        np.savez_compressed(out/'reach_timeseries.npz',time_since_flood_h=t_h,elapsed_s=sec,
            x_from_outlet_km=np.arange(0,201,10),flow_m3s=qsel,stage_m=zsel,depth_m=ysel,
            area_m2=asel,velocity_mps=vsel,water_surface_slope_downstream=slope_down,
            volume_m3=volume,storage_change_m3=volume-volume[[0],:],
            balance_residual_sampled_m3=residual,
            balance_residual_native_m3=native_residual if native is not None else np.empty(0))
        # Compact CSVs for immediate review. Full 2001-section results remain in native HDF.
        csvwrite(out/'reach_storage.csv',['t_since_flood_h']+[f'V_{k}km_m3' for k in range(0,201,10)],
                 np.column_stack([t_h,volume]))
        csvwrite(out/'control_sections.csv',['t_since_flood_h','x_from_outlet_km','flow_m3s','stage_m','depth_m','area_m2','velocity_mps','water_surface_slope_downstream'],
            ([t_h[i],j*10,qsel[i,j],zsel[i,j],ysel[i,j],asel[i,j],vsel[i,j],slope_down[i,j]]
             for i in range(len(sec)) for j in range(21)))
    dump(out/'quality.json',metrics)
    return metrics

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','prepare-body','worker','run','extract'])
    ap.add_argument('case',nargs='?'); args=ap.parse_args()
    if args.action=='prepare':
        events=prepare_design()
        mid=event();worst=event(6000,96,.20,6,-.25,.30,.40,'PILOT_ADVERSE',th=.30)
        cases=[('PILOT_mid_var_100m_5s_v3',mid,'var',100,5,None),
               ('PILOT_mid_ref_100m_5s_v3',mid,'ref',100,5,None),
               ('PILOT_adverse_var_100m_5s_v3',worst,'var',100,5,None),
               ('PILOT_adverse_ref_100m_5s_v3',worst,'ref',100,5,None)]
        for c,e,k,dx,dt,end in cases: build(c,e,k,dx,dt,end)
        manifest=[]
        for e in events:
            # Archive authoritative boundary arrays before any result is computed.
            end=math.ceil((max(e['D_h'],e['Hend_h'])+168)/24)*24
            t=np.arange(round((ORIGIN_H+end)*60)+1)/60-ORIGIN_H
            q,h=boundary(e,t,'var')
            np.savez_compressed(ROOT/'boundaries'/f"{e['event_id']}.npz",time_since_flood_h=t,inflow_m3s=q,stage_var_m=h,stage_ref_m=np.full_like(h,5.))
            for kind in ['var','ref']:
                manifest.append(dict(case=f"{e['event_id']}_{kind}_100m_5s",event_id=e['event_id'],kind=kind,split=e['split'],status='pending_quality_gate'))
        dump(ROOT/'config/batch_manifest.json',dict(created=now(),events_sha256=sha(ROOT/'config/events.json'),cases=manifest))
        print(json.dumps(dict(status='prepared',events=len(events),pilot_projects=len(cases),body_cases=len(manifest))))
    elif args.action=='prepare-body':
        events=prepare_design(); count=0
        for e in events:
            for kind in ['var','ref']:
                build(f"{e['event_id']}_{kind}_100m_5s",e,kind)
                count+=1
            if count%20==0: print(json.dumps(dict(prepared_body_projects=count,total=200)),flush=True)
        print(json.dumps(dict(prepared_body_projects=count,total=200)),flush=True)
    elif args.action=='worker': worker(ROOT/'models'/args.case)
    elif args.action=='run':
        folder=ROOT/'models'/args.case;run(folder);print(json.dumps(extract(folder)))
    else: print(json.dumps(extract(ROOT/'models'/args.case)))

if __name__=='__main__':main()
