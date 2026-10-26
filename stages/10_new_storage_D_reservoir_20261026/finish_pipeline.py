"""Current-turn bounded result pipeline; not a scheduler or future quota continuation."""
from common10 import *
import subprocess,time
def main():
 start=time.monotonic()
 while not all((P/f'{d}_completed.json').exists() for d in ['DUF','DUR','DCF','DCR']):
  if (P/'provenance/baseline_progress.json').exists():
   progress=read(P/'provenance/baseline_progress.json')
   if any(not r['complete'] for r in progress):raise RuntimeError('Baseline worker failed; preserve files and inspect')
  if time.monotonic()-start>3600:raise TimeoutError('Baseline not completed within current pipeline window')
  time.sleep(10)
 for script in ['aggregate_formal.py','audit_final_models.py','storage_audit.py','formal_figures.py','storage_figure.py','qa_figures.py','deliver_stage10.py','finalize_delivery.py']:
  print('START',script,flush=True)
  with (P/'provenance'/('final_'+script+'.log')).open('w',encoding='utf-8') as log:
   subprocess.run([sys.executable,str(P/script)],cwd=P,stdout=log,stderr=subprocess.STDOUT,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
  print('DONE',script,flush=True)
 print('LOCAL_PIPELINE_COMPLETE',flush=True)
if __name__=='__main__':main()
