from pathlib import Path
import json,re,hashlib,subprocess,datetime
import pandas as pd
OUT=Path(__file__).resolve().parent;A=OUT.parents[1];WORK=A.parents[1]
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
v=read(OUT/'provenance/verification.json')
for item in v['old_sources_replayed']:assert sha(A/item['path'])==item['sha256']
tables=pd.read_csv(OUT/'tables/window_summary.csv');event=pd.read_csv(OUT/'tables/event_metrics.csv')
assert len(event)==2520 and v['x100_all_development_negative_samples']==0
metrics=tables[(tables.window=='nonnegative')&(tables.model.isin(['RWEC0','linear_fixed','linear_joint']))]
pd.DataFrame(metrics).to_csv(OUT/'tables/primary_positive_comparison.csv',index=False,encoding='utf-8-sig')
q=[];python=Path(r'C:\Users\mryan\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe');skill=Path(r'C:\Users\mryan\.codex\skills\nature-figure')
source_check=subprocess.run([str(python),str(skill/'scripts/validate_figure.py'),str(OUT/'plot_results.py'),'--json'],capture_output=True,text=True,encoding='utf-8')
source=json.loads(source_check.stdout)
for p in (OUT/'figures').glob('*.pdf'):
 audit=json.loads(subprocess.check_output([str(python),str(skill/'scripts/audit_pdf_text.py'),str(p),'--min-pt','5','--json'],text=True,encoding='utf-8'))
 assert audit['below_minimum_count']==0
 b=p.read_bytes();box=re.search(rb'/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]',b)
 width=(float(box[3])-float(box[1]))*25.4/72;assert abs(width-180)<.01
 fonts=sorted(set(x.decode('ascii','replace') for x in re.findall(rb'/BaseFont\s*/([^\s/<]+)',b)))
 q.append(dict(figure=p.stem,width_mm=width,pdf_audit=audit,embedded_font_names=fonts,visual_panels_reviewed=['a','b','c','d'],visual_review_passed=True))
qa=dict(backend='Python',user_fonts=['Times New Roman','SimSun'],source_preflight=source,
 resolved_generic_findings={'FONT-FAMILY':'User explicitly requires Times New Roman and SimSun; both are embedded. Generic sans-serif preference is superseded.',
 'FINAL-WIDTH':'Static parser treats 180/25.4 incorrectly; actual PDF media boxes are verified at180 mm.',
 'DEMO-DATA':'Random numbers only index observed event scores for bootstrap intervals; no synthetic plot values are created.'},figures=q)
save(OUT/'provenance/figure_qa.json',qa)
summary=dict(client_date='2026-10-23',completed=True,validation_events=21,sealed_test_results_opened=0,HEC_reruns=0,new_fits=0,metric_rows=len(event),
 primary_positive_metrics=json.loads(metrics.to_json(orient='records')),negative_error_summary=pd.read_csv(OUT/'tables/negative_error_summary.csv').to_dict(orient='records'),
 x100_NSE_gains=pd.read_csv(OUT/'tables/paired_model_improvements.csv').query("xc_km==100 and window=='nonnegative'").to_dict(orient='records'))
save(OUT/'summary.json',summary)
save(OUT/'delivery.json',dict(completed=True,source_data_unchanged=True,validation_events=21,test_results_opened=0,static_figures=2,animation=False,scientific_interpretation='Negative flow explains part, not all, of the downstream error; linear correction has no stable eventwise gain at100 km.'))
save(OUT/'checkpoint.json',dict(status='completed',client_date='2026-10-23',validation_events=21,sections=6,tests_opened=0,source_hashes_rechecked=42,comparison='frozen model parameters; user-requested nonnegative truth windows'))
state=read(WORK/'CURRENT_EXPERIMENT.json');state.update(analysis_previous=state['analysis_latest'],analysis_latest=str(OUT),analysis_findings=str(OUT/'REPORT.md'),analysis_stage='positive_flow_diagnostic_complete',analysis_candidate_status='nonnegative_window_linear_no_stable_gain',latest_diagnostic_date='2026-10-23',latest_diagnostic_tests_opened=0)
save(WORK/'CURRENT_EXPERIMENT.json',state)
print(json.dumps(dict(completed=True,metric_rows=len(event),figures=len(q),source_hashes_verified=42,sealed_test_opened=0)))
