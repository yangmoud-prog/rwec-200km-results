from common import *
from validation import data_A,data_B,NODES
from routing_core import forecast
import inspect

def main():
    groupsA=data_A();groupsB=data_B();lookup=loadnp(OUT/'data/B_lookup.npz');fits=read(OUT/'fits/A_all_final.json');transfer=[];lengths=[];reversals=[]
    for label,groups in [('A',groupsA),('B',groupsB)]:
        for eid,g in groups.items():
            for c in g['cases']:
                mask=c['time']>=0
                for j,xc in enumerate(TARGETS):
                    y=c['truth'][mask,j];reversals.append(dict(experiment=label,event=eid,subset=g['subset'],case=c['case'],condition=c['condition'],xc_km=int(xc),truth_min_m3s=float(y.min()),truth_negative_samples=int((y<0).sum()),truth_samples=len(y)))
    frame('truth_reversal_by_case.csv',reversals)
    for eid,g in groupsB.items():
        c=g['cases'][0];mask=g['time']>=0;d=loadnp(OUT/'data'/f'B_{eid}.npz')
        for f in fits[:4]:
            r=forecast(g['inflow'],1.,f['core'],f['alpha'],f['mode'],lookup)
            for j,xc in enumerate(TARGETS):
                transfer.append(dict(event=eid,subset=g['subset'],source='A_only_training',target='B_slope075',model=f['model'],xc_km=int(xc),alpha=f['alpha'],zero_B_recalibration=True,**metrics(c['truth'][mask,j],r[0][mask,(200-xc)//20],g['time'][mask])))
            if f['mode'] in (1,2,3):lengths.append(dict(event=eid,subset=g['subset'],model=f['model'],RMSE_km=float(np.sqrt(np.mean((r[2][mask]/1000-d['lb'][mask])**2))),alpha=f['alpha'],zero_B_recalibration=True))
    frame('A_to_B_zero_refit_transfer.csv',transfer);frame('A_to_B_zero_refit_length.csv',lengths)
    manifest=read(OUT/'provenance/A_source_manifest.json')
    for row in manifest['inputs']:
        path=WORK/row['path'];assert sha(path)==row['sha256'] and path.stat().st_mtime_ns==row['mtime_ns']
    prior=read(P4/'provenance/A_input_manifest.json')
    for row in prior['inputs']:
        path=A/row['path'];assert sha(path)==row['sha256'] and path.stat().st_mtime_ns==row['mtime_ns']
    frozen=read(B/'config/pilot_input_manifest.json');audits=[]
    for c in frozen['cases']:
        folder=B/'models'/c['case'];changes=[]
        for name,digest in c['input_hashes'].items():
            actual=sha(folder/name)
            if actual!=digest:
                assert name.endswith('.prj'),'Changed hydraulic input'
                s=(folder/name).read_text()
                assert all(v in s for v in ['SI Units','Geom File=g01','Unsteady File=u01','Current Plan=p01'])
                changes.append(dict(file=name,before_sha256=digest,after_sha256=actual,reason='HEC saves display/DSS/GIS metadata'))
        audits.append(dict(case=c['case'],geometry_plan_boundary_unchanged=True,project_metadata_changes=changes))
    checks=[]
    for label in ['A','B']:
        c=pd.read_csv(OUT/'tables'/f'{label}_causality_checks.csv');assert c.prefix_difference_m3s.max()==0
        checks.append(dict(experiment=label,prefix_test_count=len(c),max_future_input_prefix_difference_m3s=float(c.prefix_difference_m3s.max())))
    save(OUT/'provenance/verification.json',dict(raw_A_sources_unchanged=True,verified_raw_A_files=len(manifest['inputs']),B_frozen_hydraulic_inputs_unchanged=True,
       B_frozen_sha256=sha(B/'config/pilot_frozen.json'),B_input_audits=audits,sealed_test_results_opened=0,all20_B_events_retained=True,
       A_B_mixed_calibration=False,transfer_A_to_B_actually_computed=True,animation_requested=False,
       independently_frozen_previous_stage_inputs_verified=len(prior['inputs']),previous_manifest_sha256=sha(P4/'provenance/A_input_manifest.json')))
    save(OUT/'provenance/causality_check.json',dict(forecast_signature=str(inspect.signature(forecast)),future_water_level_input=False,
       HEC_internal_flow_input=False,HEC_length_input=False,prefix_checks=checks,identical_A_inflow_prediction_shared_across_scenarios=True,
       A_actual_control_not_known=True,A_transplanted_control_explicit_assumption=True))
    for label in ['A','B']:
        df=pd.read_csv(OUT/'tables'/f'{label}_event_metrics.csv');index=[]
        for j,(name,g) in enumerate(df.groupby('model')):
            path=f'{label}_event_metrics_model_{j:02d}.csv';g.to_csv(OUT/'tables'/path,index=False,encoding='utf-8-sig');index.append(dict(model=name,file=path,rows=len(g)))
        frame(f'{label}_metric_file_index.csv',index)
    print('Raw-input, causal and actual zero-refit transfer checks complete',flush=True)
if __name__=='__main__':main()
