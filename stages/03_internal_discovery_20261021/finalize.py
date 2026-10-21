"""Register completion only after numerical and rendered-figure review."""
import argparse,datetime,json
from pathlib import Path
import pandas as pd
from prepare import OUT,ROOT,read,save,sha

def run(visual=False):
    assert visual,'Rendered panels must be reviewed before finalization.'
    verification=read(OUT/'provenance/verification.json');assert verification['passed'] and verification['pdf_count']==8
    assert verification['raw_npz_files_unchanged']==172 and verification['test_result_files_opened']==0
    audit=read(OUT/'provenance/source_figure_audit.json')
    assert not [x for x in audit['findings'] if x['level']=='FAIL' and x['check_id']!='FONT-FAMILY']
    contracts=read(OUT/'provenance/figure_contracts.json');assert len(contracts)==8 and sum(c['panels'] for c in contracts)==28
    panels=[]
    for c in contracts:
        for letter in 'abcd'[:c['panels']]:
            panels.append(dict(figure=c['figure'],panel=letter,visual_review='passed',
                labels_units_legend_and_collisions='checked at rendered final dimensions',
                evidence_claim=c['claim'],source_data_traceable=True))
    save(OUT/'provenance/figure_qa.json',dict(visual_review_completed=True,panels=panels,
        source_font_exception='Only FONT-FAMILY overridden by explicit user Times New Roman / SimSun requirement; actual embedded PDF fonts checked',
        warning_review='Chinese legend text is not subject to English title-case rules; no non-font source FAIL remains',
        uncertainty='Per-event IQR or full observed range defined consistently in corresponding caption; no pseudo-independent length tests',
        animation_review='Beginning, rise/hold, maximum suppression, maximum recovery and final frames inspected; axes fixed and shared physical time verified'))
    report=(OUT/'REPORT.md').read_text(encoding='utf-8')
    assert not any(ord(c)<32 and c not in '\n\r\t' for c in report)
    assert '\\frac' in report and '\\tau' in report and '{{' not in report
    assert read(OUT/'fits/symbolic_candidates.json')['status']=='completed'
    completed=datetime.datetime.now().astimezone().isoformat()
    state='second_stage_complete_waiting_researcher_review'
    sourcehash={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*.py') if '.deps' not in p.parts and '.julia_env' not in p.parts and '__pycache__' not in p.parts}
    sourcehash.update({p.name:sha(p) for p in OUT.glob('*.md')})
    delivery=dict(status=state,completed=completed,stage='internal-section screening and conservative equation discovery',
        main_internal_xc_km=100,auxiliary_hydrodynamic_only_xc_km=150,screened_cases=172,
        training_events=49,development_internal_validation_events=21,sealed_test_events_used=0,new_hec_runs=0,
        primary_models=['RWEC0','RWEC1','P1','P2','NN20261021','SR1','SR2','SR3'],
        additional_robustness=['RWEC0J','NN20261022'],pysr_successful=True,
        static_figures=8,static_panels=28,video='animation/E007_backwater_100km.mp4',gif='animation/E007_backwater_preview.gif',
        primary_report='REPORT.md',reproduction='REPRODUCE.md',runtime_verification='provenance/verification.json',
        scientific_conclusion='100km saved series direction pass; NN and SR2/SR3 improve development responses; no single equation dominates all diagnostics; fixed-parameter timestep transfer remains unresolved',
        candidate_frozen_for_external_test=False,next_action_requires_user_decision=True,
        original_extracted_data_unchanged=True,source_sha256=sourcehash)
    save(OUT/'delivery.json',delivery)
    save(OUT/'checkpoint.json',dict(status=state,updated=completed,test_results_opened=0,raw_sources_unchanged=True,
        scientific_next_step='Researcher review before any sealed testing or formula extension',auto_continue_to_test=False))
    routepath=ROOT.parents[1]/'CURRENT_EXPERIMENT.json'
    before=OUT/'provenance/current_experiment_before_stage03_completion.json'
    if not before.exists():before.write_bytes(routepath.read_bytes())
    route=read(routepath);assert Path(route['active_root'])==ROOT
    route.update(updated=completed,analysis_stage=state,analysis_latest=str(OUT),analysis_findings=str(OUT/'REPORT.md'),
        internal_prediction_x_km=100,internal_training_events=49,internal_validation_events=21,
        model_comparison_development_events=70,model_comparison_test_events=0,next_action_requires_user_decision=True,
        analysis_candidate_status='development_hypotheses_with_unresolved_timestep_transfer')
    save(routepath,route)
    print(json.dumps({k:v for k,v in delivery.items() if k!='source_sha256'},ensure_ascii=False,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--visual-review-confirmed',action='store_true');a=p.parse_args();run(a.visual_review_confirmed)
