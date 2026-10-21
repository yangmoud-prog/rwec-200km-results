from common import *
import sys,re,subprocess,ast
from PIL import Image,ImageDraw,ImageFont

def main():
    skill=Path(r'C:\Users\mryan\.codex\skills\nature-figure/scripts')
    source=subprocess.run([sys.executable,str(skill/'validate_figure.py'),str(OUT/'figures.py'),'--json'],capture_output=True,text=True,encoding='utf-8')
    preflight=__import__('json').loads(source.stdout)
    assert [x['check_id'] for x in preflight['findings'] if x['level']=='FAIL']==['FONT-FAMILY']
    audits=[]
    for p in sorted((OUT/'figures').glob('*.pdf')):
        b=p.read_bytes();names=sorted(set(x.decode() for x in re.findall(rb'/BaseFont\s*/([^\s/>]+)',b)))
        box=re.search(rb'/MediaBox\s*\[([^]]+)\]',b).group(1).split();width=float(box[2])*25.4/72
        assert abs(width-180)<.02 and all('TimesNewRoman' in x or 'SimSun' in x for x in names)
        result=subprocess.run([sys.executable,str(skill/'audit_pdf_text.py'),str(p),'--min-pt','5','--json'],capture_output=True,text=True,encoding='utf-8')
        assert result.returncode==0,result.stdout
        audits.append(dict(file=p.name,width_mm=width,base_fonts=names,glyph_audit=__import__('json').loads(result.stdout)))
    assert len(audits)==48
    for p in OUT.glob('*.py'):ast.parse(p.read_text(encoding='utf-8'))
    for name in ['REPORT.md','EQUATIONS.md']:
        p=OUT/name
        if p.exists():
            text=p.read_text(encoding='utf-8');assert '\t' not in text and text.count('$$')%2==0
            assert not any(ord(c)<32 and c not in '\n\r' for c in text)
    for label in ['A','B']:
        df=pd.read_csv(OUT/'tables'/f'{label}_event_metrics.csv');assert np.isfinite(df[['NSE','RMSE_m3s','KGE']]).all().all()
    # Contact sheets for inspecting every gallery page without changing source figures.
    paths=sorted((OUT/'figures').glob('A_validation_*.png'))+sorted((OUT/'figures').glob('B_flood_*.png'))
    font=ImageFont.truetype(r'C:\Windows\Fonts\times.ttf',16)
    for start in range(0,len(paths),15):
        tile=Image.new('RGB',(1200,1500),'white');draw=ImageDraw.Draw(tile)
        for i,p in enumerate(paths[start:start+15]):
            im=Image.open(p).convert('RGB');im.thumbnail((390,265));x=(i%3)*400;y=(i//3)*300
            tile.paste(im,(x,y+25));draw.text((x+5,y+4),p.stem,fill='black',font=font)
        tile.save(OUT/'data'/f'gallery_QA_{start//15}.png')
    save(OUT/'provenance/figure_qa.json',dict(backend='Python',source_preflight=preflight,
        user_font_override='Times New Roman English/numbers/symbols; SimSun Chinese, overriding generic sans-serif validator',
        warning_width='source parser misreads180/25.4; actual MediaBoxes all180mm',Chinese_capitalization_warning_inapplicable=True,
        PDF_audits=audits,static_main_figures=7,A_validation_gallery_events=21,B_gallery_events=20,animation_requested=False,
        uncertainty_definition='event interquartile range; no stochastic confidence intervals',all_events_retained=True,raw_negative_flows_not_clipped=True))
    print('48 PDFs and all static source/data checks passed, with explicit requested-font override',flush=True)
if __name__=='__main__':main()
