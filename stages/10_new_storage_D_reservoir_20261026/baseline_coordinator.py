from common10 import *
import subprocess,time,os,psutil
from concurrent.futures import ThreadPoolExecutor,as_completed

def ready(domain,names,cfg):
 geom,bc={'DUF':('uniform','fixed'),'DUR':('uniform','reservoir'),'DCF':('complex','fixed'),'DCR':('complex','reservoir')}[domain]
 for e in names:
  folder=D/'models'/f"{e}_{geom}_{bc}_{cfg['version']}"
  if not (folder/'quality.json').exists():return False
  q=read(folder/'quality.json')
  if not q['numerical_gate']:raise RuntimeError('Numerical quality failure '+folder.name)
 return True

def job(domain):
 start=time.monotonic()
 while not (D/'config/D_frozen.json').exists():
  if time.monotonic()-start>7200:raise TimeoutError('No frozen D design')
  time.sleep(10)
 cfg=read(D/'config/D_frozen.json')
 while not ready(domain,cfg['train'],cfg):
  if time.monotonic()-start>10800:raise TimeoutError('Training HEC not ready')
  time.sleep(10)
 with (P/'provenance'/f'baseline_{domain}.log').open('a',encoding='utf-8') as f:
  # Adopt already-running own calibration workers if a coordinator is restarted.
  while True:
   active=[]
   for proc in psutil.process_iter(['pid','cmdline']):
    cmd=proc.info.get('cmdline') or []
    if proc.pid!=os.getpid() and any(str(P/'formal_baselines.py').lower()==str(arg).lower() for arg in cmd) and domain in cmd and 'fit' in cmd:active.append(proc.pid)
   if not active:break
   if time.monotonic()-start>14400:raise TimeoutError('Adopted fit exceeded coordinator timeout')
   time.sleep(10)
  if not (P/'provenance'/f'{domain}_calibration_frozen.json').exists():
   r=subprocess.run([sys.executable,str(P/'formal_baselines.py'),domain,'fit'],cwd=P,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
   if r.returncode:raise RuntimeError('Training failed '+domain)
  while not ready(domain,cfg['validation'],cfg):
   if time.monotonic()-start>14400:raise TimeoutError('Validation HEC not ready')
   time.sleep(10)
  r=subprocess.run([sys.executable,str(P/'formal_baselines.py'),domain,'evaluate'],cwd=P,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
  if r.returncode:raise RuntimeError('Evaluation failed '+domain)
 return dict(domain=domain,complete=True,elapsed_s=time.monotonic()-start)

if __name__=='__main__':
 rows=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  fs={pool.submit(job,d):d for d in ['DUF','DUR','DCF','DCR']}
  for f in as_completed(fs):
   try:r=f.result()
   except Exception as ex:r=dict(domain=fs[f],complete=False,error=repr(ex))
   rows.append(r);save(P/'provenance/baseline_progress.json',rows);print(json.dumps(r),flush=True)
