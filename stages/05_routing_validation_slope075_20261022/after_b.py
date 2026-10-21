from common import *
import time,subprocess,sys
while True:
    p=B/'logs/pilot_checkpoint.json'
    if p.exists() and read(p).get('finished',0)==40:break
    time.sleep(15)
for script,args in [('validation.py',['B']),('diagnostics.py',['B'])]:
    result=subprocess.run([sys.executable,str(OUT/script),*args],cwd=WORK)
    if result.returncode:
        save(OUT/'logs/B_postprocess_failure.json',dict(script=script,returncode=result.returncode,updated=now()))
        print('B postprocessing error recorded; HEC/numerical jobs continue',flush=True)
        # Continue the second independent postprocessing attempt if usable fits exist.
print('B postprocessing coordinator finished',flush=True)
