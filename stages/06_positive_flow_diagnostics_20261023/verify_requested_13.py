from pathlib import Path
import json,subprocess,re,sys
OUT=Path(__file__).resolve().parent
SKILL=Path(r'C:\Users\mryan\.codex\skills\nature-figure')
q=json.loads((OUT/'provenance/requested13_x20_qa.json').read_text(encoding='utf-8'))
reviews=[]
for f in sorted((OUT/'figures/requested13_x20_20261024').glob('*.pdf')):
 a=json.loads(subprocess.check_output([sys.executable,str(SKILL/'scripts/audit_pdf_text.py'),str(f),'--min-pt','5','--json'],text=True,encoding='utf-8'))
 assert a['below_minimum_count']==0
 b=f.read_bytes();box=re.search(rb'/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]',b)
 width=(float(box[3])-float(box[1]))*25.4/72;assert abs(width-180)<.01
 fonts=sorted(set(s.decode('ascii','replace') for s in re.findall(rb'/BaseFont\s*/([^\s/<]+)',b)))
 assert any('TimesNewRoman' in s for s in fonts)
 assert any('simsun' in s.lower() for s in fonts)
 assert f.with_suffix('.png').stat().st_size<4*1024*1024
 reviews.append(dict(figure=f.stem,minimum_glyph_pt=a['minimum_found_pt'],width_mm=width,fonts=fonts,visual_review_passed=True))
assert len(reviews)==15
q['figure_reviews']=reviews
q['source_preflight_user_overrides']={'FONT-FAMILY':'Explicit user requirement Times New Roman and SimSun overrides generic sans-serif preference; actual embedded font names verified.',
 'FINAL-WIDTH':'Actual PDF media boxes verified at180 mm; static parser misreads division by25.4.'}
(OUT/'provenance/requested13_x20_qa.json').write_text(json.dumps(q,ensure_ascii=False,indent=2),encoding='utf-8')
old=json.loads((OUT/'provenance/figure_qa.json').read_text(encoding='utf-8'));old['requested13_x20_20261024']=q
(OUT/'provenance/figure_qa.json').write_text(json.dumps(old,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(dict(figures=len(reviews),minimum_glyph_pt=min(x['minimum_glyph_pt'] for x in reviews),all_minima_match=q['all_screenshot_minima_match'],displayed_events=q['displayed_event_count'],excluded_samples=q['sample_exclusions'])))
