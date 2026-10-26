from common10 import *
from formulas import *
def main():
 rows=[]
 for a,b,n in [(1,.1,1),(1,3,1),(1,20,2),(.1,1,1),(3,1,1),(3,1,2),(.3,4,.7),(2,5,1.7)]:
  r=equilibrium(a,b,n);root=(r/(r+1))**n+(r/(a*r+b))**n-1
  analytic=(a-1+np.sqrt((a-1)**2+4*b))/2 if n==1 else r
  rho=(r/(r+1))**(n+1)+a*(r/(a*r+b))**(n+1)
  q,s,v,ok,step,bal=rwec(np.ones(1001)*1000,.5,a,b,n,1.3,np.ones(10)*1e8)
  err=np.max(abs(q-1000));passed=ok and err<1e-6 and abs(root)<1e-10 and abs(r-analytic)<1e-7 and 0<rho<1
  assert passed,(a,b,n,r,err)
  rows.append(dict(a=a,b=b,n=n,r=r,rho=rho,constant_Q_error_m3s=err,root_residual=root,mass_balance_max_m3=bal,passed=bool(passed)))
 # Causal pulse, no future-state dependency, all formula families.
 t=np.arange(1001)*.5;inp=1000+4500*np.exp(-((t-75)/25)**2);inp[t==0]=1000
 for a,b in [(1,3),(3,1),(2,4)]:
  q,s,v,ok,step,bal=rwec(inp,.5,a,b,1.5,1.3,np.ones(10)*1e8);changed=inp.copy();changed[500:]*=1.5
  qc,*_=rwec(changed,.5,a,b,1.5,1.3,np.ones(10)*1e8)
  assert ok and q.min()>=0 and s.min()>=0 and np.max(abs(q[:,:500]-qc[:,:500]))==0
 frame('formula_unit_checks.csv',rows);save(P/'provenance/verification.json',dict(steady_cases_passed=len(rows),nonnegative_causal_pulse_cases_passed=3,math='unique stable positive equilibrium for a,b,n>0',numerical_condition_limit_r=1e10))
 print('FORMULA_CHECKS_PASS',flush=True)
if __name__=='__main__':main()
