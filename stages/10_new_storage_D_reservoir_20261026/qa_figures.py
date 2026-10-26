from common10 import *
from PIL import Image,ImageOps,ImageDraw,ImageFont
sys.path.insert(0,'C:/Users/mryan/.codex/skills/nature-figure/scripts')
from audit_pdf_text import audit_pdf
def main():
 rows=[];files=sorted((P/'figures').glob('*.png'))
 for f in files:
  im=Image.open(f);a=audit_pdf(f.with_suffix('.pdf').read_bytes(),5);svg=f.with_suffix('.svg').read_text(encoding='utf-8');row=dict(figure=f.stem,width_px=im.width,height_px=im.height,dpi=im.info.get('dpi'),PDF_minimum_pt=a['minimum_found_pt'],PDF_text_auditable=a['auditable'],PDF_below_5pt=a['below_minimum_count'],SVG_editable_text='<text' in svg,Times_New_Roman_declared='Times New Roman' in svg,SimSun_declared='SimSun' in svg)
  row['passed']=bool(a['auditable'] and not a['below_minimum_count'] and row['SVG_editable_text'] and row['Times_New_Roman_declared'] and row['SimSun_declared'] and min(row['dpi'])>=299)
  rows.append(row)
 for page in range((len(files)+5)//6):
  sheet=Image.new('RGB',(1800,1770),'white');draw=ImageDraw.Draw(sheet);font=ImageFont.truetype('C:/Windows/Fonts/times.ttf',23)
  for j,f in enumerate(files[page*6:(page+1)*6]):
   im=Image.open(f).convert('RGB');im.thumbnail((870,535));x=(j%2)*900+(900-im.width)//2;y=(j//2)*590+40;sheet.paste(im,(x,y));draw.text(((j%2)*900+15,(j//2)*590+10),f.stem,font=font,fill='black')
  sheet.save(P/'provenance'/f'contact_sheet_{page+1}.png')
 save(P/'provenance/figure_qa.json',dict(figures=len(rows),all_machine_checks_passed=all(r['passed'] for r in rows),font_override='User requested Times New Roman / SimSun; supersedes generic journal Arial suggestion.',output_purpose='scientific research reading figures; PNG 300 dpi, editable SVG, vector PDF; not final journal-size layout',rows=rows,visual_inspection='contact sheets and full-size key figures inspected by agent; see delivery notes'))
 print('FIGURE_QA',len(rows),all(r['passed'] for r in rows))
if __name__=='__main__':main()
