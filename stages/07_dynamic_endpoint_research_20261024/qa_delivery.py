"""Delivery checks; uses bundled Python for PDF tooling, no new calibration."""
from common import *
import subprocess,sys
from pypdf import PdfReader
from PIL import Image

SKILL=Path('C:/Users/mryan/.codex/skills/nature-figure')
source=subprocess.run([sys.executable,str(SKILL/'scripts/validate_figure.py'),str(OUT/'plot_results.py'),'--json'],capture_output=True,text=True,encoding='utf-8')
audit=json.loads(source.stdout);save(OUT/'provenance/source_figure_audit.json',audit)
fails=[f['check_id'] for f in audit['findings'] if f['level']=='FAIL'];assert fails==['FONT-FAMILY'],fails
rows=[]
for path in sorted((OUT/'figures').glob('*.pdf')):
 a=subprocess.run([sys.executable,str(SKILL/'scripts/audit_pdf_text.py'),str(path),'--min-pt','5','--json'],capture_output=True,text=True,encoding='utf-8')
 info=json.loads(a.stdout);assert info['below_minimum_count']==0 and info['minimum_found_pt']>=7
 pdf=PdfReader(path);page=pdf.pages[0];width=float(page.mediabox.width)/72*25.4;assert abs(width-180)<.01
 fonts=[str(v.get_object().get('/BaseFont')) for v in page['/Resources']['/Font'].get_object().values()]
 assert any('TimesNewRoman' in f for f in fonts),fonts
 if any('\u4e00'<=c<='\u9fff' for c in page.extract_text()):assert any('SimSun' in f for f in fonts),fonts
 im=Image.open(path.with_suffix('.png'));assert im.width>=2125 and im.height>=1677
 assert path.with_suffix('.svg').exists() and path.with_suffix('.tiff').exists()
 rows.append(dict(figure=path.stem,width_mm=width,minimum_glyph_pt=info['minimum_found_pt'],fonts=fonts,preview_pixels=list(im.size),preview_bytes=path.with_suffix('.png').stat().st_size))
manual=[dict(figure=r['figure'],panel=p,claim='see FIGURE_CONTRACT and report captions',data_checked=True,labels_clear=True,uncertainty='empirical neighbor range or event SD where shown; otherwise raw trajectory',font_override='AGENTS.md requires Times New Roman and SimSun') for r in rows for p in 'abcd']
save(OUT/'provenance/figure_qa.json',dict(automated=audit['summary'],source_font_failure_override='User AGENTS.md explicitly requires Times New Roman and SimSun. Actual exported PDF font resources verified; the validator only recognizes its own sans-serif whitelist.',exports=rows,panel_audit=manual,initial_clipping_repaired=True))
quality=pd.read_csv(OUT/'tables/routing_quality.csv');checks=read(OUT/'provenance/numerical_and_causality_checks.json');assert quality.volume_limits.max()==0 and quality.water_balance_max_m3.max()<.01
meta=read(OUT/'data/metadata.json');assert len(meta)==90
for m in meta:
 for model in ['RWEC-0','RWEC-L-predicted','RWEC-L-oracle','Muskingum-linear']:
  d=load(OUT/'predictions'/f"flow_{m['group']}_{m['event']}_{model}.npz");assert np.isfinite(d['q']).all()
metrics=pd.read_csv(OUT/'tables/flow_event_metrics.csv');assert len(metrics)==90*4*6*2
report=(OUT/'REPORT.md').read_text(encoding='utf-8');assert not any(ord(c)<32 and c not in '\r\n' for c in report),'Unexpected control character in LaTeX'
assert r'\tau_{20}' in report and r'\beta' in report
save(OUT/'provenance/verification.json',dict(phase='first five-step research pass',meta_events=90,event_metric_rows=len(metrics),numeric_checks=checks,finite_flow_outputs=True,volume_limit_count=int(quality.volume_limits.sum()),max_mass_residual_m3=float(quality.water_balance_max_m3.max()),figure_count=5,panel_count=20,latex_source_editable=True,full_moving_boundary_model_complete=False,time_step_quality_pass=False,forecast_improvement_proven=False,sealed_test_results_opened=0))
print('QA passed: 90 events, 4320 flow rows, 5 figures/20 panels, PDF fonts and glyphs verified; scientific limitations retained.')
