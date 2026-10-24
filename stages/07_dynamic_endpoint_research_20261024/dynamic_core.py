"""Minimal conservative dynamic-delay complementary routing, SI volumes."""
import numpy as np
from numba import njit

@njit(cache=True)
def complement(a,b,n):
 lo=min(a,b);hi=max(a,b)
 return lo/np.exp(np.log1p(np.exp(n*np.log(lo/hi)))/n)

@njit(cache=True)
def initial_ratio(k,n):
 lo=0.;hi=1.e8
 for _ in range(100):
  r=.5*(lo+hi)
  if (r/(r+1))**n+(r/(r+k))**n>1:hi=r
  else:lo=r
 return .5*(lo+hi)

@njit(cache=True)
def interp(q,j,i,z,u):
 if z<=0:return q[0,j]
 k=int(np.floor(z));w=z-k
 if k>=i:return u
 b=u if k+1>=i else q[k+1,j]
 return q[k,j]+w*(b-q[k,j])

@njit(cache=True)
def route(inflow,length,dt,k,n,tau,beta):
 nt=len(inflow);sec=dt*3600.;q=np.full((nt,11),inflow[0]);q[:,0]=inflow
 s=np.full(10,initial_ratio(k,n)*inflow[0]*sec);transit=np.empty(10);clock=np.empty(10)
 for j in range(10):
  r=max(0.,min(length[0],(10-j)*20.)-(9-j)*20.)/20.
  lag=tau*(1+beta*r);clock[j]=-lag/dt;transit[j]=inflow[0]*lag*3600.
 initial=s.sum()+transit.sum()+sec*.5*q[0,1:10].sum();totalin=0.;totalout=0.;maxres=0.;limits=0;mintransit=transit.min()
 ss=np.zeros((nt,10));tt=np.zeros((nt,10));ss[0]=s;tt[0]=transit
 for i in range(1,nt):
  for j in range(10):
   r=max(0.,min(length[i],(10-j)*20.)-(9-j)*20.)/20.
   lag=tau*(1+beta*r);newclock=max(clock[j],i-lag/dt);newclock=min(float(i),newclock)
   a=interp(q,j,i,clock[j],q[i,j]);b=interp(q,j,i,newclock,q[i,j])
   v=.5*(a+b)*(newclock-clock[j])*sec
   raw=.5*(q[i-1,j]+q[i,j])*sec;available=transit[j]+raw
   if v>available+1.e-7:limits+=1;v=max(0.,available)
   sw=s[j]+v;se=s[j]+k*v;new= complement(sw,se,n)
   q[i,j+1]=(sw-new)/sec;s[j]=new;transit[j]=available-v;clock[j]=newclock
  ss[i]=s;tt[i]=transit;mintransit=min(mintransit,transit.min())
  totalin+=.5*(inflow[i-1]+inflow[i])*sec;totalout+=q[i,10]*sec
  state=s.sum()+transit.sum()+sec*.5*q[i,1:10].sum()
  maxres=max(maxres,abs(state-initial-totalin+totalout))
 return q,ss,tt,maxres,limits,mintransit

@njit(cache=True)
def muskingum(inflow,dt,K,X):
 nt=len(inflow);q=np.full((nt,11),inflow[0]);q[:,0]=inflow
 D=2*K*(1-X)+dt;c0=(dt-2*K*X)/D;c1=(dt+2*K*X)/D;c2=(2*K*(1-X)-dt)/D
 for i in range(1,nt):
  for j in range(10):q[i,j+1]=c0*q[i,j]+c1*q[i-1,j]+c2*q[i-1,j+1]
 return q

@njit(cache=True)
def batch_loss(theta,inputs,lengths,truth,masks,denoms,kind):
 total=0.
 if kind==3:K=np.exp(theta[0]);X=theta[1]*.5/K
 else:k=np.exp(theta[0]);n=np.exp(theta[1]);tau=theta[2];beta=0. if kind==0 else theta[3]
 nodes=np.array([9,8,7,6,5,4])
 for ev in range(len(inputs)):
  if kind==3:pred=muskingum(inputs[ev],1.,K,X)
  else:pred=route(inputs[ev],lengths[ev],1.,k,n,tau,beta)[0]
  for j in range(6):
   v=0.
   for i in range(49,289):
    if masks[ev,i,j]:v+=(pred[i,nodes[j]]-truth[ev,i,j])**2
   total+=v/denoms[ev,j]
 return total/(len(inputs)*6)
