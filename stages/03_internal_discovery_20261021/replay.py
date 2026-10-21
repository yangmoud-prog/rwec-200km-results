"""Replay saved candidates without HEC runs, refitting or symbolic searching."""
from pathlib import Path
import os,subprocess,sys
OUT=Path(__file__).resolve().parent
os.environ.setdefault('PYTHONIOENCODING','utf-8')
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
for args in [ ['verify.py','native'],['evaluate.py','full'],['evaluate.py','sensitivity'],['verify.py','solver'],
    ['plot_results.py'],['animation/make_animation.py'],['write_report.py'],['verify.py','delivery']]:
    subprocess.run([sys.executable,str(OUT/args[0]),*args[1:]],cwd=OUT,check=True)
