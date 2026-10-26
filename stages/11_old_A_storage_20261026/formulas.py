"""Stage 10 exact SE alternatives; independent equilibrium, causal nonnegative routing."""
import numpy as np
from numba import njit
from legacy_models import at,muskingum

@njit(cache=True)
def equilibrium(a,b,n):
 # log residual avoids cancellation when r/(r+1) is near one.
 lo=-35.;hi=np.log(1.e10)
 for k in range(100):
  u=(lo+hi)/2;inv=np.exp(-u)
  g=np.log(-np.expm1(-n*np.log1p(inv)))+n*np.log(a+b*inv)
  if g>0:lo=u
  else:hi=u
 r=np.exp((lo+hi)/2)
 if r>1.e10*(1-1.e-8):return np.nan
 return r

@njit(cache=True)
def rwec(inflow,dt,a,b,n,tau,v0,ncells=10):
 nt=len(inflow);sec=dt*3600.;q=np.full((ncells+1,nt),inflow[0]);q[0]=inflow
 s=np.zeros((ncells,nt));transit=np.zeros((ncells,nt));vol=np.zeros((ncells,nt))
 r=equilibrium(a,b,n)
 if not np.isfinite(r):return q,s,vol,False,-1,np.inf
 s0=r*inflow[0]*sec;s[:,0]=s0;transit[:,0]=inflow[0]*tau*3600;vol[:,0]=v0
 maxbal=0.
 for i in range(1,nt):
  for j in range(ncells):
   raw=.5*(q[j,i-1]+q[j,i])*sec
   v=.5*(at(q[j,:i+1],i-1-tau/dt)+at(q[j,:i+1],i-tau/dt))*sec
   transit[j,i]=transit[j,i-1]+raw-v
   w=s[j,i-1]+v;e=a*s[j,i-1]+b*v
   if w<=0 or e<=0 or transit[j,i]<-1.e-4:return q,s,vol,False,i,np.inf
   logratio=n*np.log(w/e)
   exponent=(max(0.,logratio)+np.log1p(np.exp(-abs(logratio))))/n
   out=w*(-np.expm1(-exponent))/sec
   s[j,i]=w-out*sec;q[j+1,i]=out
   vol[j,i]=v0[j]+s[j,i]-s0+transit[j,i]-transit[j,0]+.5*sec*(out-inflow[0])
   residual=vol[j,i]-vol[j,i-1]-(raw-.5*(q[j+1,i-1]+out)*sec)
   maxbal=max(maxbal,abs(residual))
   if not np.isfinite(out) or out<0 or out>1.e6 or vol[j,i]<=0 or s[j,i]<0:return q,s,vol,False,i,maxbal
 return q,s,vol,True,0,maxbal

@njit(cache=True)
def simulate(theta,model,inp,dt,v0,ncells=10):
 if model==3:
  q,vol,ok,step,bal=muskingum(inp,dt,np.exp(theta[0]),theta[1],theta[2],ncells=ncells);return q,vol,vol,ok,step,bal
 if model==0:a=1.;b=np.exp(theta[0]);n=np.exp(theta[1]);tau=theta[2]
 elif model==1:a=np.exp(theta[0]);b=1.;n=np.exp(theta[1]);tau=theta[2]
 else:a=np.exp(theta[0]);b=np.exp(theta[1]);n=np.exp(theta[2]);tau=theta[3]
 return rwec(inp,dt,a,b,n,tau,v0,ncells)

@njit(cache=True)
def hourly(q,dt,hours):
 # Exact integration of the piecewise-linear nodal discharge, including dt=2h.
 out=np.zeros((q.shape[0],hours));integ=np.zeros_like(q)
 for i in range(1,q.shape[1]):integ[:,i]=integ[:,i-1]+.5*(q[:,i]+q[:,i-1])*dt
 prev=np.zeros(q.shape[0])
 for h in range(1,hours+1):
  k=min(int(h/dt),q.shape[1]-2);rem=h-k*dt
  cur=integ[:,k]+rem*q[:,k]+.5*rem*rem*(q[:,k+1]-q[:,k])/dt
  out[:,h-1]=cur-prev;prev=cur
 return out

@njit(cache=True)
def objective(theta,model,dt,inputs,truth,den,v0):
 total=0.
 for ev in range(len(inputs)):
  q,s,vol,ok,step,bal=simulate(theta,model,inputs[ev],dt,v0[ev])
  if not ok:return 1000.+(len(inputs)-ev)/len(inputs)
  pred=hourly(q,dt,truth.shape[2])[4:11]
  total+=np.sum(np.sum((pred-truth[ev])**2,axis=1)/den[ev])
 return total/(len(inputs)*7)
