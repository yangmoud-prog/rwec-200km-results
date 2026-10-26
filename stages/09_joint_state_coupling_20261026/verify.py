from common import *
from joint import simulate
from datetime import datetime

def main():
 cfg=read(P/'provenance/frozen_protocol.json');fit=read(P/'PARAMETERS.json')['Joint'];d=load(P/'data/B_E023.npz');dt=fit['dt_h'];par=np.exp(fit['theta']);t=np.arange(0,288+dt/2,dt);inp=np.interp(t,GRID,d['inflow']);zero=np.zeros(len(t))
 original=simulate(inp,dt,*par,d['yinit'],zero,0);changed=inp.copy();changed[t>96]+=300
 future=simulate(changed,dt,*par,d['yinit'],zero,0)
 assert np.array_equal(original[0][:,t<=96],future[0][:,t<=96])
 # A hostile future length array must have no influence in autonomous mode.
 hostile=simulate(inp,dt,*par,d['yinit'],np.ones(len(t))*199,0)
 assert np.array_equal(original[0],hostile[0])
 steady=simulate(np.ones(len(t))*1000,dt,*par,d['yinit'],zero,0)
 assert steady[-2]==0 and np.max(abs(steady[0]-1000))<1e-7
 q=pd.read_csv(P/'tables/solver_quality.csv');sens=pd.read_csv(P/'tables/fixed_parameter_dt_sensitivity.csv');fixed=pd.read_csv(P/'tables/fixed_section_event_metrics.csv');manifest=read(P/'provenance/source_manifest.json')
 assert len(fixed)==20*7*7
 assert len(set(cfg['train'])&set(cfg['validation']))==0
 for source in manifest['sources']:assert sha(W/source['path'])==source['sha256']
 causal=dict(future_upstream_suffix_changes_past_outputs=False,autonomous_ignores_future_endpoint_array=True,autonomous_future_HEC_stage_or_storage_inputs=False,initial_state_only_HEC_input=True,serial8_external_preprocessor=True,oracle_separately_labelled=True,C_not_used_to_fit_B=True,sealed_test_opened=0)
 save(P/'provenance/causality_check.json',causal)
 checks=dict(checked=datetime.now().isoformat(),status='first_development_round_complete_candidate_not_supported',events=20,train=14,validation=6,models=7,fixed_section_event_metric_rows=len(fixed),base_completed_runs=int((q.failed_step==0).sum()),base_total_runs=len(q),base_negative_reconstructed_flow_count=int(q.negative_count.sum()),base_max_balance_flow_m3s=float(q.max_balance_flow_m3s.max()),joint_validation_all_complete=bool((q.query("model=='Joint' and subset=='validation'").failed_step==0).all()),sensitivity_failed_event_model_dt_runs=int(len(sens[~sens.complete][['event','model','dt_factor']].drop_duplicates())),steady_flow_max_error_m3s=float(np.max(abs(steady[0]-1000))),source_hashes_rechecked=True,causality=causal,no_new_HEC=True,no_animation=True,uncertainty='six development validation events; sealed30 not used; no independent final-test claim')
 save(P/'provenance/verification.json',checks);print(json.dumps(checks,ensure_ascii=False,indent=2))
 print(sens.groupby(['model','dt_factor']).agg(NSE=('NSE','mean'),completed=('complete','sum'),total=('complete','count')).to_string())
if __name__=='__main__':main()
