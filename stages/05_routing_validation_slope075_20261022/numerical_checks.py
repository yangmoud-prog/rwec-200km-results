from common import *
import importlib.util,time,concurrent.futures,traceback
spec=importlib.util.spec_from_file_location('hec075_num',B/'scripts/hec200.py');hec=importlib.util.module_from_spec(spec);spec.loader.exec_module(hec)
frozen=read(B/'config/pilot_frozen.json');e=frozen['events']['E023']
cases=[]
for kind in ['control','ref']:
    for dx,dt in [(2000,30),(1000,60)]:
        case=f'CHECK_E023_{kind}_{dx}m_{dt}s';folder=hec.build(case,e,kind,dx=dx,dt=dt)
        cases.append(dict(case=case,dx=dx,dt=dt,kind=kind,folder=folder))
save(B/'config/numerical_checks_frozen.json',dict(cases=[{k:v for k,v in c.items() if k!='folder'} for c in cases],reason='representative E023 time and space checks; physical geometry and control unchanged',prepared=now()))
while True:
    checkpoint=B/'logs/pilot_checkpoint.json'
    if checkpoint.exists() and read(checkpoint).get('finished',0)==40:break
    time.sleep(15)
def run(c):
    try:
        hec.run(c['folder']);q=hec.extract(c['folder'])
        base=loadnp(B/'results'/f"E023_{c['kind']}_2km_60s"/'reach_timeseries.npz');d=loadnp(B/'results'/c['case']/'reach_timeseries.npz')
        assert np.allclose(base['time_since_flood_h'],d['time_since_flood_h'])
        rows=[]
        for xc in TARGETS:
            j=int(xc//10);o=d['flow_m3s'][:,j];p=base['flow_m3s'][:,j];mask=d['time_since_flood_h']>=0
            rows.append(dict(case=c['case'],kind=c['kind'],dx_m=c['dx'],dt_s=c['dt'],xc_km=int(xc),quality_pass=q['basic_quality_pass'],
              base_vs_finer_RMSE_m3s=float(np.sqrt(np.mean((p[mask]-o[mask])**2))),peak_difference_pct=float(100*(p[mask].max()/o[mask].max()-1)),
              max_stage_difference_m=float(abs(base['stage_m'][:,j]-d['stage_m'][:,j]).max())))
        return dict(case=c['case'],status='computed',quality=q,rows=rows)
    except Exception as exc:return dict(case=c['case'],status='failed_preserved',error=repr(exc),trace=traceback.format_exc(),rows=[])
results=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
    for r in pool.map(run,cases):
        results.append(r);print(r['case'],r['status'],flush=True)
save(B/'logs/numerical_checks.json',dict(cases=results,updated=now()))
frame('B_HEC_numerical_sensitivity.csv',[row for r in results for row in r['rows']])
print('HEC time/space checks finished',flush=True)
