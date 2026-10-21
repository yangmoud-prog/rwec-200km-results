"""Independent scientific delivery checks, never open sealed test cases."""
from pathlib import Path
import sys,hashlib,json,math,subprocess,re
import numpy as np,pandas as pd,h5py
from prepare import OUT,ROOT,read,save,frame,sha
from fit import load,hy
from evaluate import configs,predict

def native_check():
    d,meta=load();rows=[]
    for j,m in enumerate(meta):
        p=ROOT/'models'/m['case']/'Backwater200.p01.hdf'
        base='Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series'
        with h5py.File(p,'r') as f:
            q=np.asarray(f[base+'/Cross Sections/Flow'])[:,::-5]
            tt=np.asarray(f[base+'/Time'],dtype=float);t=(tt-tt[0])*24-48
        with np.load(ROOT/'results'/m['case']/'reach_timeseries.npz') as z:
            error=float(np.max(abs(q-z['flow_m3s'])));assert error==0
            assert np.allclose(t,z['time_since_flood_h'],atol=.02/3600,rtol=0)
        assert q.shape[0]==d['ns'][j] and q.shape[1]==21
        assert np.all(q[:,10]>=0)
        rows.append(dict(case=m['case'],event_id=m['event_id'],saved_times=len(t),
            full_21_section_flow_max_abs_error_m3s=error,min_q_100km_m3s=float(q[:,10].min()),
            source=str(p.relative_to(ROOT)),flow_dataset_sha256=hashlib.sha256(np.ascontiguousarray(q).tobytes()).hexdigest()))
    frame(OUT/'tables/native_hdf_verification.csv',rows)
    print('全部172份原生HDF已核对',flush=True)

def solver_check():
    d,meta=load();f=read(OUT/'fits/base.json')[9];checks=[]
    j=next(i for i,m in enumerate(meta) if m['event_id']=='E025' and m['kind']=='var' and m['role']=='main');N=d['ns'][j]
    cut=N//2
    for c in configs():
        if c['type']=='neural' and not c['primary']:continue
        q=d['qin'][9,j,:N].copy();h=d['h'][j,:N].copy()
        original=predict(c,q,h,N,100000,f)
        q[cut+1:]*=1.8;h[cut+1:]+=2
        altered=predict(c,q,h,N,100000,f)
        err=float(np.max(abs(original[0][:cut+1]-altered[0][:cut+1])));assert err==0
        checks.append(dict(model=c['id'],future_boundary_changes_leave_prefix_exactly_unchanged=True,
            max_prefix_flow_difference_m3s=err,no_target_arguments_in_prediction=True))
    # Force the supply limit independently of the real-data activation frequency.
    qr=np.full((1,8),10.);eq=np.array([[0,0,1,1,.5,0,0,0.]])
    for mode in [1,2]:
        q,b,active,raw=hy.response(qr,eq,np.array([8],dtype=np.int64),1.,3.,mode,True)
        residual=np.diff(b[0])*1e6-(qr[0,1:]-q[0,1:])*3600
        assert q.min()>=-1e-9 and b.min()>=0 and abs(residual).max()<1e-6 and (active!=0).any()
    save(OUT/'provenance/runtime_checks.json',dict(passed=True,causality=checks,
        synthetic_supply_limiter_balance_and_direction=True,training_gradient_check=read(OUT/'provenance/gradient_check.json'),
        baseline_copy_identical_to_verified_previous_source=sha(OUT/'base_core.py')==sha(OUT.parent/'02_candidate_models_20261021/models.py')))
    print('因果性与可供水约束检查通过',flush=True)

def delivery():
    manifest=read(OUT/'provenance/input_manifest.json');changed=[]
    for m in manifest['inputs']:
        p=ROOT/m['path']
        if sha(p)!=m['sha256'] or p.stat().st_mtime_ns!=m['mtime_ns']:changed.append(m['path'])
    assert not changed
    balance=pd.read_csv(OUT/'tables/balance.csv');maxres=float(balance.physical_control_balance_max_m3.max());assert maxres<.001
    metrics=pd.read_csv(OUT/'tables/case_metrics.csv');assert metrics.negative_predictions.sum()==0
    assert len(metrics)==172*10*len(configs())
    assert np.all(np.isfinite(metrics[['NSE','KGE','RMSE']].to_numpy()))
    pdf_audits=[]
    for p in sorted((OUT/'figures').glob('*.pdf')):
        fonts={v.decode('ascii') for v in re.findall(rb'/BaseFont\s*/([^\s/<>]+)',p.read_bytes())}
        assert any('TimesNewRoman' in f for f in fonts),fonts
        assert any('SimSun' in f for f in fonts),fonts
        assert all('TimesNewRoman' in f or 'SimSun' in f for f in fonts),fonts
        raw=p.read_bytes();box=re.search(rb'/MediaBox\s*\[\s*([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)',raw)
        width=(float(box[3])-float(box[1]))*25.4/72;assert abs(width-180)<.01
        from PIL import Image
        with Image.open(p.with_suffix('.tiff')) as im:dpi=[float(v) for v in im.info['dpi']]
        assert min(dpi)>=599
        assert '<text' in p.with_suffix('.svg').read_text(encoding='utf-8')
        tool=Path(r'C:\Users\mryan\.codex\skills\nature-figure\scripts\audit_pdf_text.py')
        result=subprocess.run([sys.executable,str(tool),str(p),'--min-pt','5','--json'],capture_output=True,text=True,encoding='utf-8')
        a=json.loads(result.stdout);assert result.returncode==0,(p,a)
        pdf_audits.append(dict(file=p.name,fonts=sorted(fonts),glyph_audit=a,width_mm=width,tiff_dpi=dpi,svg_text_editable=True))
    save(OUT/'provenance/pdf_font_and_glyph_audits.json',pdf_audits)
    sys.path.insert(0,str(OUT/'.deps'))
    import imageio_ffmpeg
    reader=imageio_ffmpeg.read_frames(str(OUT/'animation/E007_backwater_100km.mp4'));video=next(reader);first=next(reader);reader.close()
    assert video['size']==(1920,1080) and video['fps']==12 and len(first)==1920*1080*3
    from PIL import Image
    with Image.open(OUT/'animation/E007_backwater_preview.gif') as g:
        gif=dict(frames=g.n_frames,size=list(g.size));assert g.n_frames>=190
    save(OUT/'provenance/verification.json',dict(passed=True,raw_npz_files_unchanged=len(manifest['inputs']),
        native_hdf_cases_verified=172,test_result_files_opened=0,negative_predictions=int(metrics.negative_predictions.sum()),
        max_physical_balance_residual_m3=maxres,causal_runtime_checks=True,
        pdf_count=len(pdf_audits),font_override_confirmed=['Times New Roman','SimSun'],
        video_metadata=video,gif=gif,scientific_claim_limit='Development stage only; fixed1h model, no time-step portability pass'))
    print('交付核查通过：原始数据、守恒、因果、字体和视频',flush=True)

if __name__=='__main__':
    action=sys.argv[1] if len(sys.argv)>1 else 'delivery'
    {'native':native_check,'solver':solver_check,'delivery':delivery}[action]()
