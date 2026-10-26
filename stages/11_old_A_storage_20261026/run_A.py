from common10 import *
import subprocess
from concurrent.futures import ThreadPoolExecutor,as_completed
def job(name):
 if not (P/'provenance'/f'{name}_frozen.json').exists():
  with (P/'logs'/f'fit_{name}.log').open('a',encoding='utf-8') as f:subprocess.run([sys.executable,str(P/'fit_A.py'),name],cwd=P,stdout=f,stderr=subprocess.STDOUT,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
 return dict(model=name,complete=True)
def main():
 rows=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for ft in as_completed([pool.submit(job,n) for n in MODELS]):
   r=ft.result();rows.append(r);save(P/'provenance/training_progress.json',rows);print(r,flush=True)
 save(P/'provenance/all_parameters_frozen.json',dict(complete=True,models=MODELS,validation_metrics_inspected=False))
 print('ALL_A_PARAMETERS_FROZEN',flush=True)
if __name__=='__main__':main()
