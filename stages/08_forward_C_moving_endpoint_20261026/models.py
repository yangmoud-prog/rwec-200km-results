"""Conservative routing. SI storage, hours for time parameters, x upstream from outlet."""
import numpy as np
from numba import njit

@njit(cache=True)
def complement(w,e,n):
 lo=min(w,e);hi=max(w,e)
 if lo<=0:return np.nan
 return lo/np.exp(np.log1p(np.exp(n*np.log(lo/hi)))/n)

@njit(cache=True)
def equilibrium_ratio(k,n):
 lo=0.;hi=1.e10
 for _ in range(100):
  r=(lo+hi)/2
  if (r/(r+1))**n+(r/(r+k))**n>1:hi=r
  else:lo=r
 return (lo+hi)/2

@njit(cache=True)
def at(a,z):
 if z<=0:return a[0]
 k=int(z)
 if k>=len(a)-1:return a[-1]
 return a[k]+(z-k)*(a[k+1]-a[k])

@njit(cache=True)
def geom_extra(b_km,alpha,bed_slope,y0):
 # Fixed geometry closure: delta_y = alpha S0 [b-x]+.
 # Integrate quadratic trapezoidal area exactly over each linear subinterval.
 out=np.zeros(10);b=b_km*1000.
 for j in range(10):
  lo=(9-j)*20000.;hi=(10-j)*20000.;top=min(hi,b)
  if top<=lo:continue
  h1=alpha*bed_slope*(b-lo);h2=alpha*bed_slope*(b-top)
  out[j]=(top-lo)*((500+6*y0[j])*(h1+h2)/2+(h1*h1+h1*h2+h2*h2))
 return out

@njit(cache=True)
def rwec(inflow,b,dt,k,n,tau,alpha,relax_h,bed_slope,y0,initial_volumes):
 nt=len(inflow);sec=dt*3600.;q=np.full((11,nt),inflow[0]);q[0]=inflow
 ss=np.zeros((10,nt));tt=np.zeros((10,nt));bb=np.zeros((10,nt));vol=np.zeros((10,nt))
 r0=equilibrium_ratio(k,n)*inflow[0]*sec
 ss[:,0]=r0;tt[:,0]=inflow[0]*tau*3600.;bb[:,0]=geom_extra(b[0],alpha,bed_slope,y0)
 vol[:,0]=initial_volumes
 maxbal=0.;maxlocal=0.;valid=True;failstep=0
 for i in range(1,nt):
  target=geom_extra(b[i],alpha,bed_slope,y0)
  fac=1. if relax_h<=0 else -np.expm1(-dt/relax_h)
  for j in range(10):
   raw=.5*(q[j,i-1]+q[j,i])*sec
   lagged=.5*(at(q[j,:i+1],i-1-tau/dt)+at(q[j,:i+1],i-tau/dt))*sec
   tt[j,i]=tt[j,i-1]+raw-lagged
   bb[j,i]=bb[j,i-1]+fac*(target[j]-bb[j,i-1])
   effective=lagged-(bb[j,i]-bb[j,i-1])
   w=ss[j,i-1]+effective;e=ss[j,i-1]+k*effective
   if w<=0 or e<=0 or tt[j,i]<-1.e-4:
    return q,ss,tt,bb,vol,False,i,np.inf,np.inf
   ss[j,i]=complement(w,e,n);q[j+1,i]=(w-ss[j,i])/sec
   if not np.isfinite(q[j+1,i]) or q[j+1,i]<0 or q[j+1,i]>1.e6:
    return q,ss,tt,bb,vol,False,i,np.inf,np.inf
   # J converts the native RWEC end-step mean outflow to a trapezoidal flux convention.
   vol[j,i]=initial_volumes[j]+ss[j,i]-r0+tt[j,i]-tt[j,0]+bb[j,i]-bb[j,0]+sec*.5*(q[j+1,i]-q[j+1,0])
   residual=vol[j,i]-vol[j,i-1]-(raw-.5*(q[j+1,i-1]+q[j+1,i])*sec)
   maxlocal=max(maxlocal,abs(residual))
   if vol[j,i]<=0:return q,ss,tt,bb,vol,False,i,np.inf,np.inf
  residual=(vol[:,i]-vol[:,i-1]).sum()-.5*(q[0,i-1]+q[0,i]-q[10,i-1]-q[10,i])*sec
  maxbal=max(maxbal,abs(residual))
 return q,ss,tt,bb,vol,valid,failstep,maxbal,maxlocal

@njit(cache=True)
def muskingum(inflow,dt,tref,x,m,qref=1000.):
 nt=len(inflow);q=np.full((11,nt),np.nan);q[0]=inflow
 kp=tref/(m*qref**(m-1));vol=np.zeros((10,nt));maxbal=0.
 for j in range(10):
  u=q[j];out=q[j+1];out[0]=u[0];s=kp*u[0]**m;vol[j,0]=s*3600
  for i in range(1,nt):
   available=s+.5*dt*(u[i-1]+u[i]-out[i-1])
   zero=kp*max(0.,x*u[i])**m
   if available-zero < -1.e-10*max(1.,abs(available),abs(zero)):return q,vol,False,i,np.inf
   lower=0.;upper=max(u[i],out[i-1],2*available/dt,1.);cur=min(out[i-1],upper)
   for iteration in range(45):
    w=max(1.e-12,x*u[i]+(1-x)*cur);v=kp*w**m+.5*dt*cur-available
    if abs(v)<=1.e-12*max(1.,abs(available)):break
    if v<0:lower=cur
    else:upper=cur
    derivative=kp*m*(1-x)*w**(m-1)+.5*dt;candidate=cur-v/derivative
    if not np.isfinite(candidate) or candidate<=lower or candidate>=upper:candidate=(lower+upper)/2
    cur=candidate
   snew=kp*max(0.,x*u[i]+(1-x)*cur)**m
   err=(snew-s-.5*dt*(u[i-1]+u[i]-out[i-1]-cur))*3600.
   if not np.isfinite(cur) or cur<0 or abs(err)>1.e-5*max(1.,snew*3600):return q,vol,False,i,np.inf
   maxbal=max(maxbal,abs(err));out[i]=cur;s=snew;vol[j,i]=s*3600
 return q,vol,True,0,maxbal

@njit(cache=True)
def loss(theta,model,dt,inputs,lengths,truth,denoms,bed_slope,y0,initial_volumes):
 total=0.;nodes=np.array([4,5,6,7,8,9,10]);nttruth=truth.shape[2]
 for ev in range(len(inputs)):
  if model in (0,3,4):
   alpha=0. if model==0 else theta[3];relax=0. if model==0 else theta[4]
   q,ss,tt,bb,vol,ok,step,err,local=rwec(inputs[ev],lengths[ev],dt,np.exp(theta[0]),np.exp(theta[1]),theta[2],alpha,relax,bed_slope,y0,initial_volumes)
  else:
   m=1. if model==1 else theta[2]
   q,vol,ok,step,err=muskingum(inputs[ev],dt,np.exp(theta[0]),theta[1],m)
  if not ok:return 1000.+(len(inputs)-ev)/len(inputs)
  for j in range(7):
   sse=0.
   for t in range(nttruth):
    d=at(q[nodes[j]],t/dt)-truth[ev,j,t];sse+=d*d
   total+=sse/denoms[ev,j]
 return total/(len(inputs)*7)
