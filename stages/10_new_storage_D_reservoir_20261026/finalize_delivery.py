from common10 import *
import datetime,platform,importlib.metadata
def main():
 summary=read(P/'SUMMARY.json');qa=read(P/'provenance/figure_qa.json');assert qa['all_machine_checks_passed']
 env=dict(python=sys.version,executable=sys.executable,platform=platform.platform(),HEC_RAS='6.3.1',packages={k:importlib.metadata.version(k) for k in ['numpy','scipy','pandas','numba','h5py','matplotlib','pywin32','psutil']})
 save(P/'EXECUTION_ENVIRONMENT.json',env)
 cfg=read(P/'provenance/frozen_protocol.json');sources=[]
 for domain in ['B','C']:
  for e in cfg['train']+cfg['validation']:
   f=(P9/'data'/f'B_{e}.npz') if domain=='B' else P8/'data'/f'{e}.npz';sources.append(dict(domain=domain,event=e,path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),access='existing development data only'))
 save(P/'SOURCE_DATA_MANIFEST.json',dict(existing=sources,D_cases_manifest=str(D/'config/formal_cases.json'),D_input_audit='D_INPUT_AUDIT.json',old_sealed_test_access=False))
 quality=read(P/'provenance/verification.json');quality.update(stage10_delivery=summary['status'],D_checks=summary['checks'],causality=read(P/'provenance/causality_check.json'),HEC_formal_cases=80,HEC_passed=80,oracle_closure_failed=16,autonomous_coupling_not_validated=True,old_sealed_test_access=False,time_alignment=read(P/'provenance/time_alignment_verified.json'),figure_machine_QA_passed=True,LaTeX_compiled_with_desktop_editor=True,storage_and_interface_audit=read(P/'STORAGE_AUDIT.json'))
 save(P/'provenance/verification.json',quality)
 delivery=read(P/'DELIVERY.json');delivery.update(completed=True,complete_means='initial formula, HEC, baseline and physical-gate result package delivered; autonomous coupling remains scientifically unvalidated',full_research_complete=False);save(P/'DELIVERY.json',delivery)
 save(P/'checkpoint.json',dict(stage=10,status=summary['status'],HEC_cases=80,baseline_domains=4,baseline_models=4,validation_events=6,train_events=14,old_sealed_tests_read=0,oracle_profile_candidates=2,oracle_cases=16,oracle_passed=0,autonomous_coupling_run=False,feedback_ablation_run=False,remaining='Replace failed quasi-steady one-shape closure with physically justified causal shape/storage memory before autonomous joint forecasting.',no_new_automation=True))
 current=read(W/'CURRENT_EXPERIMENT.json');current.update(updated=datetime.datetime.now().astimezone().isoformat(),analysis_stage='stage10_formula_D_baselines_delivered_physical_coupling_gate_failed',analysis_findings=str(P/'REPORT.md'),analysis_latest=str(P),analysis_latest_status='stage_result_delivery_local_complete_repository_sync_pending',analysis_candidate_status='new_SE_not_consistently_better_two_profile_closures_failed',stage10_summary=str(P/'SUMMARY.json'),stage10_D_root=str(D),stage10_HEC_formal_cases=80,stage10_HEC_passed=80,stage10_autonomous_coupling_validated=False,requested_computations_complete=False,full_research_complete=False,forecast_window_h=[0,720],future_downstream_levels_used_as_predictor=False,next_action_requires_user_decision=False);save(W/'CURRENT_EXPERIMENT.json',current)
 print('LOCAL_DELIVERY_VERIFIED')
if __name__=='__main__':main()
