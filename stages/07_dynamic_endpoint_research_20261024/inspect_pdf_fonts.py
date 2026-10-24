from pathlib import Path
from pypdf import PdfReader
import sys
sys.stdout.reconfigure(encoding='utf-8')
p=Path(__file__).parent/'figures/01_endpoint_and_storage.pdf';f=PdfReader(p);page=f.pages[0]
print('text extract',page.extract_text()[:1500]);print('resources',page['/Resources'])
for key,v in page['/Resources'].items():
 if key=='/XObject':
  ob=v.get_object();print('Xobjects',list(ob)[:10])
  for k in list(ob)[:3]:print(k,str(ob[k].get_object())[:500])
for path in sorted(p.parent.glob('*.pdf')):
 pg=PdfReader(path).pages[0];print(path.name,[str(v.get_object().get('/BaseFont')) for v in pg['/Resources']['/Font'].get_object().values()])
