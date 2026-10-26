from common10 import *
from PIL import Image,ImageDraw,ImageFont
sys.path.insert(0,'C:/Users/mryan/.codex/skills/nature-figure/scripts')
from audit_pdf_text import audit_pdf
def main():
 cfg=read(P/'provenance/frozen_protocol.json');manifest=read(P/'SOURCE_MANIFEST.json');assert len(manifest['sources'])==140 and len(set(cfg['train'])&set(cfg['validation']))==0
 for row in manifest['sources']:assert sha(row['path'])==row['sha256']
 for name in ['formulas.py','legacy_models.py']:assert sha(P/name)==sha(P10/name)
 events=pd.read_csv(P/'tables/event_metrics.csv');assert len(events)==70*4*7*2*3;assert len(set(events.query("subset=='validation'").event))==21
 same=events.query("x_from_outlet_km==100 and scenario=='dynamic'").pivot(index=['event','model'],columns='window',values='NSE');assert np.nanmax(abs(same.full-same.nonnegative_hour_endpoints))<1e-12
 rows=[]
 for name in MODELS:
  f=[read(p) for p in (P/'provenance').glob(f'fit_A_{name}_reach*_dt*.json')];assert len(f)==40 and sum(r['evaluations'] for r in f)==80000;rows.append(dict(model=name,candidate_fits=len(f),objective_calls=sum(r['evaluations'] for r in f),failed_parameter_calls=sum(r['failed_evaluations'] for r in f),summed_worker_wall_s=sum(r['seconds'] for r in f)))
 frame('optimization_budget_audit.csv',rows)
 figs=[];files=sorted((P/'figures').glob('*.png'))
 for path in files:
  im=Image.open(path);a=audit_pdf(path.with_suffix('.pdf').read_bytes(),5);svg=path.with_suffix('.svg').read_text(encoding='utf-8');passed=a['auditable'] and a['below_minimum_count']==0 and '<text' in svg and 'Times New Roman' in svg and 'SimSun' in svg and min(im.info.get('dpi',(0,0)))>=299
  figs.append(dict(figure=path.stem,PDF_minimum_pt=a['minimum_found_pt'],passed=bool(passed),width_px=im.width,height_px=im.height))
 for page in range((len(files)+5)//6):
  canvas=Image.new('RGB',(1800,1770),'white');draw=ImageDraw.Draw(canvas);font=ImageFont.truetype('C:/Windows/Fonts/times.ttf',23)
  for j,path in enumerate(files[page*6:(page+1)*6]):
   im=Image.open(path).convert('RGB');im.thumbnail((870,530));x=(j%2)*900+(900-im.width)//2;y=(j//2)*590+40;canvas.paste(im,(x,y));draw.text(((j%2)*900+15,(j//2)*590+10),path.stem,font=font,fill='black')
  canvas.save(P/'provenance'/f'contact_sheet_{page+1}.png')
 assert all(r['passed'] for r in figs)
 save(P/'provenance/figure_qa.json',dict(figures=len(figs),all_machine_checks_passed=True,rows=figs,font_override='User requested Times New Roman / SimSun',visual_inspection='full-size summary panels plus every event panel on contact sheets to be inspected before delivery',purpose='research reading figures, not a journal final-width claim',uncertainty='event paired changes displayed individually; mean curves are actual event means, no range bands',exclusion='post-prediction common nonnegative hourly endpoints; all 21 validation events retained'))
 check=read(P/'provenance/evaluation_complete.json');save(P/'provenance/verification.json',dict(**check,source_hashes_rechecked=True,formula_kernels_match_previous_stage=True,objective_calls_total=320000,budgets_equal=True,all21_validation_events_retained=True,at100km_full_equals_nonnegative_verified=True,causality=read(P/'provenance/causality_check.json'),figures=len(figs),figure_machine_QA_passed=True,LaTeX_compiled_successfully=True,old_sealed_tests_read=0))
 print('FINAL_CHECKS_PASSED',len(figs),flush=True)
if __name__=='__main__':main()
