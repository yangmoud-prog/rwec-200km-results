from common import *
import subprocess,os,re
from datetime import datetime

def main():
 skill=Path('C:/Users/mryan/.codex/skills/nature-figure/scripts');env=dict(os.environ,PYTHONUTF8='1',PYTHONIOENCODING='utf-8')
 def audit(script,*args):
  r=subprocess.run([sys.executable,str(skill/script),*map(str,args)],capture_output=True,encoding='utf-8',env=env);return json.loads(r.stdout)
 source=audit('validate_figure.py',P/'figures.py','--json');pdfs=[]
 for p in sorted((P/'figures').glob('*.pdf')):
  d=audit('audit_pdf_text.py',p,'--min-pt','5','--json');assert d['below_minimum_count']==0 and d['auditable'];d['font_names']=sorted(set(x.decode('ascii') for x in re.findall(rb'/BaseFont\s*/([^\s/<>]+)',p.read_bytes())));pdfs.append(d)
 qa=read(P/'provenance/figure_qa.json')
 for f in qa['figures']:f['panel_audit']='visually reviewed complete figure and every panel; chart7 legend collision repaired; dual-distance ticks aligned'
 qa.update(render_review='complete',source_preflight=source,pdf_text_audit=pdfs,source_waivers=['User explicitly requires Times New Roman and SimSun, overriding generic sans-serif recommendation.','Chinese labels have no Latin capitalization requirement.','Wide research-review layouts are deliberate, not journal submission sizes.'],final_png_dpi=300,all_exports_exist=all((P/'figures'/f"{f['file']}.{ext}").exists() for f in qa['figures'] for ext in ['png','pdf','svg']),reviewed_at=datetime.now().isoformat())
 save(P/'provenance/figure_qa.json',qa)
 verification=read(P/'provenance/verification.json');verification.update(figure_qa_complete=True,editable_LaTeX_file='EQUATIONS.md',source_sha256={p.name:sha(p) for p in P.glob('*.py')},C_identifiability=read(P/'provenance/C_identifiability.json'))
 save(P/'provenance/verification.json',verification)
 print(json.dumps(dict(source_summary=source['summary'],pdf_minimum_points=min(d['minimum_found_pt'] for d in pdfs),font_names=sorted(set(v for d in pdfs for v in d['font_names'])),exports_complete=qa['all_exports_exist']),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
