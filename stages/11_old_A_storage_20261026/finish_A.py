from common10 import *
import subprocess,time
def main():
 started=time.monotonic()
 while not (P/'provenance/all_parameters_frozen.json').exists():
  if time.monotonic()-started>7200:raise TimeoutError('Current-turn training did not finish')
  time.sleep(10)
 for script in ['evaluate_A.py','historical_reference.py','figures_A.py','report_A.py','final_checks.py']:
  print('START',script,flush=True)
  with (P/'logs'/('final_'+script+'.log')).open('w',encoding='utf-8') as f:subprocess.run([sys.executable,str(P/script)],cwd=P,stdout=f,stderr=subprocess.STDOUT,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
  print('DONE',script,flush=True)
 print('A_COMPARISON_AND_LOCAL_DELIVERY_COMPLETE',flush=True)
if __name__=='__main__':main()
