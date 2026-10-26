from common10 import *
from formulas import equilibrium,simulate
import time

def main():
 rows=[];checks=[]
 for domain in ['DUF','DUR','DCF','DCR']:
  for name in ['original','new','general','NLM']:
   for f in read(P/f'PARAMETERS_{domain}_{name}.json')['reaches']:
    t=f['theta'];row=dict(domain=domain,model=name,reach=f['reach'],x_from_outlet_km=f['outlet_x_km'],distance_from_inlet_km=200-f['outlet_x_km'],dt_h=f['dt_h'],training_NSE=f['training_NSE'],boundary_parameter_count=len(f['boundary_parameters']))
    if name=='NLM':row.update(tref_h=float(np.exp(t[0])),X=t[1],m=t[2])
    else:
     a,b,n,tau=(1.,np.exp(t[0]),np.exp(t[1]),t[2]) if name=='original' else (np.exp(t[0]),1.,np.exp(t[1]),t[2]) if name=='new' else (np.exp(t[0]),np.exp(t[1]),np.exp(t[2]),t[3])
     r=equilibrium(a,b,n);rho=(r/(r+1))**(n+1)+a*(r/(a*r+b))**(n+1)
     row.update(a=float(a),b=float(b),n=float(n),tau_h=tau,steady_S_per_20km_m3=float(r*1000*f['dt_h']*3600),local_memory_factor=float(rho),local_decay_h=float(-f['dt_h']/np.log(rho)))
    rows.append(row)
    inp=np.full(int(96/f['dt_h'])+1,1000.);q,s,v,ok,step,bal=simulate(np.array(t),f['model_id'],inp,f['dt_h'],np.array([1e8]),1)
    err=float(abs(q-1000).max());passed=bool(ok and err<1e-5 and s.min()>=0 and q.min()>=0)
    checks.append(dict(domain=domain,model=name,reach=f['reach'],constant_flow_error_m3s=err,storage_min_m3=float(s.min()),flow_min_m3s=float(q.min()),balance_residual_m3=float(bal),passed=passed))
 frame('D_parameters_interpreted.csv',rows);frame('D_final_parameter_steady_checks.csv',checks)
 qual=pd.read_csv(P/'tables/D_solver_quality.csv');metricsdf=pd.read_csv(P/'tables/D_event_metrics.csv');sens=pd.read_csv(P/'tables/D_fixed_parameter_dt_sensitivity.csv')
 frame('D_evaluation_completeness.csv',metricsdf.groupby(['domain','subset','model']).agg(rows=('complete','size'),completed=('complete','sum'),finite_NSE=('NSE','count')).reset_index())
 from formal_baselines import prepare,predict
 cfg=read(D/'config/D_frozen.json');runtime=[];budget=[]
 for domain in ['DUF','DUR','DCF','DCR']:
  d=prepare(domain,[cfg['validation'][0]],cfg)[0]
  for name in ['original','new','general','NLM']:
   fit=read(P/f'PARAMETERS_{domain}_{name}.json')['reaches'];predict(d,fit,cfg['end_h']);timings=[]
   for repeat in range(3):
    start=time.perf_counter();predict(d,fit,cfg['end_h']);timings.append(time.perf_counter()-start)
   runtime.append(dict(domain=domain,model=name,event=cfg['validation'][0],forecast_h=720,median_forecast_s=float(np.median(timings)),includes_IO=False,excludes_first_compilation=True))
   trials=[read(f) for f in (P/'provenance').glob(f'fit_{domain}_{name}_reach*_dt*.json')];calls=sum(f['evaluations'] for f in trials);assert len(trials)==40 and calls==80000
   budget.append(dict(domain=domain,model=name,candidate_fits=len(trials),objective_calls=calls,failed_evaluations=sum(f['failed_evaluations'] for f in trials),sum_candidate_wall_s=sum(f['seconds'] for f in trials)))
 frame('D_runtime.csv',runtime);frame('D_calibration_budget_audit.csv',budget)
 save(P/'provenance/final_model_checks.json',dict(steady_checks=len(checks),steady_all_passed=all(r['passed'] for r in checks),minimum_full_run_Q_m3s=float(qual.min_Q.min()),maximum_step_balance_m3s=float(qual.balance_m3s.max()),failed_reach_runs=int((~qual.ok).sum()),sensitivity_failed_target_rows=int((~sens.complete).sum()),validation_target_rows=int((metricsdf.subset=='validation').sum()),all_calibrated_parameters_reported=True))
 print(read(P/'provenance/final_model_checks.json'))
if __name__=='__main__':main()
