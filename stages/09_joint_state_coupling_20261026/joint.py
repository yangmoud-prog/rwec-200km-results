"""Experimental two-zone moving-domain solver. SI units; positive x upstream.

Finite-volume remapping is an assumption, not evidence of physical adequacy.
The forecast interface accepts only inflow, known initial depth and parameters.
The optional endpoint array is confined to the labelled oracle experiment.
"""
import numpy as np
from numba import njit

@njit(cache=True)
def normal_depth(q):
 y=max(.2,2.8552*(q/1000.)**.6)
 for _ in range(12):
  a=500*y+3*y*y;t=500+6*y;p=500+2*np.sqrt(10.)*y;r=a/p
  flow=a*r**(2/3)*np.sqrt(.000075)/.025
  der=flow*(t/a+2/3*(t/a-2*np.sqrt(10.)/p))
  delta=(flow-q)/der;y=max(.01,y-delta)
  if abs(delta)<1e-9:break
 return y

@njit(cache=True)
def length_scale(q,mult):
 y=normal_depth(q);a=500*y+3*y*y;t=500+6*y;p=500+2*np.sqrt(10.)*y
 fr2=q*q*t/(9.81*a*a*a)
 derivative=2*.000075*(t/a+2/3*(t/a-2*np.sqrt(10.)/p))
 return mult*(1-fr2)/derivative

@njit(cache=True)
def rating(h):
 # Same rounded 241-pair, piecewise-linear HEC downstream control.
 q=max(0.,1000*((h-3.5)/1.5)**1.5) if h>3.5 else 0.
 j=min(239,int(q/50));a=round(3.5+1.5*(j*50/1000)**(2/3),5);b=round(3.5+1.5*((j+1)*50/1000)**(2/3),5)
 return j*50+50*(h-a)/(b-a)

@njit(cache=True)
def shape(q,h,mult,yinit):
 lam=length_scale(q,mult);yn=normal_depth(q)-normal_depth(1000.)
 a=np.zeros(101);cum=np.zeros(101)
 for j in range(101):
  decay=np.exp(-j*2000/lam);y=yinit[j]+yn*(1-decay)+(h-5)*decay
  a[j]=500*y+3*y*y
  if j>0:cum[j]=cum[j-1]+1000*(a[j]+a[j-1])
 b=0. if h<=5.1 else min(200000.,lam*np.log((h-5)/.1))
 return b,a,cum

@njit(cache=True)
def atx(a,x):
 if x<=0:return a[0]
 if x>=200000:return a[-1]
 j=int(x/2000);return a[j]+(x/2000-j)*(a[j+1]-a[j])

@njit(cache=True)
def exponent(k,r):
 lo=.001;hi=100000.
 for _ in range(48):
  n=(lo+hi)/2
  val=np.exp(n*np.log(r/(r+1)))+np.exp(n*np.log(r/(r+k)))
  if val>1:lo=n
  else:hi=n
 return (lo+hi)/2

@njit(cache=True)
def comp(w,e,n):
 if w<=0 or e<=0:return np.nan
 lo=min(w,e);hi=max(w,e)
 return lo/np.exp(np.log1p(np.exp(n*np.log(lo/hi)))/n)

@njit(cache=True)
def route(qin,b,oldb,oldu,oldcum,initcum,k,Ksec,dt,ntable,offflows,mode):
 u=np.zeros(10);wtotal=0.;qq=qin;valid=True
 for j in range(10):
  hi=(10-j)*20000.;lo=hi-20000.;oldlo=max(lo,min(hi,oldb));newlo=max(lo,min(hi,b));oldlen=hi-oldlo;newlen=hi-newlo
  if newlo>oldlo:
   swept=oldu[j]*(newlo-oldlo)/oldlen if oldlen>0 else 0.
  else:swept=-(atx(oldcum,oldlo)-atx(oldcum,newlo))
  wtotal+=swept
  if newlen<1e-6:continue
  f=newlen/20000.;base=atx(initcum,hi)-atx(initcum,newlo);r0=1000*Ksec*f
  r=oldu[j]-swept-base+r0
  if mode==2:
   # Feedback off: one-way flux from an independent fixed-domain RWEC.
   qnext=offflows[j]+f*(offflows[j+1]-offflows[j])
   u[j]=oldu[j]-swept+dt*(qq-qnext)
  else:
   idx=min(999,int(f*1000));n=ntable[1000] if f>=1 else (exponent(k,Ksec*f/dt) if f<.001 else ntable[idx]+(f*1000-idx)*(ntable[idx+1]-ntable[idx]))
   w=r+qq*dt;e=r+k*qq*dt;rn=comp(w,e,n);qnext=(w-rn)/dt;u[j]=base+rn-r0
  if not np.isfinite(qnext) or qnext<0 or u[j]<0:valid=False
  qq=qnext
 return qq,u,wtotal,valid

@njit(cache=True)
def residual(q,h,qin,oldb,oldu,oldcum,oldsb,initcum,k,Ksec,mult,dt,ntable,yinit,forcedb,offflows,mode):
 b,area,cum=shape(q,h,mult,yinit)
 if mode==1:b=forcedb
 qq,u,w,ok=route(qin,b,oldb,oldu,oldcum,initcum,k,Ksec,dt,ntable,offflows,mode)
 sb=atx(cum,b);f1=q-qq;f2=(sb-oldsb-w)/dt-q+rating(h)
 return np.array([f1,f2]),b,u,cum,sb,w,ok

@njit(cache=True)
def cumulative(x,b,u,cum):
 total=atx(cum,min(x,b))
 if x<=b:return total
 for j in range(10):
  hi=(10-j)*20000.;lo=max(hi-20000.,min(hi,b));length=hi-lo
  if length>0:total+=u[j]*max(0.,min(x,hi)-lo)/length
 return total

@njit(cache=True)
def simulate(inflow,dt_h,k,K_h,mult,yinit,oracle,mode):
 dt=dt_h*3600;nt=len(inflow);Ksec=K_h*3600
 result=np.full((11,nt),np.nan);bb=np.full(nt,np.nan);hh=bb.copy();qb=bb.copy();sbout=bb.copy();suout=bb.copy();sweep=bb.copy();iters=np.zeros(nt);errout=np.zeros(nt)
 ntable=np.zeros(1001)
 for i in range(1,1001):ntable[i]=exponent(k,Ksec*(i/1000)/dt)
 b,area,cum=shape(1000.,5.,mult,yinit);initcum=cum.copy();oldcum=cum.copy();u=np.array([atx(cum,(10-j)*20000)-atx(cum,(9-j)*20000) for j in range(10)])
 h=5.;q=1000.;sb=0.;oldb=0.;oldsb=0.;fullR=np.ones(10)*1000*Ksec
 result[:,0]=1000;bb[0]=0;hh[0]=5;qb[0]=1000;sbout[0]=0;suout[0]=u.sum();sweep[0]=0;maxmass=0.;failure=0
 for ti in range(1,nt):
  qin=.5*(inflow[ti-1]+inflow[ti]);off=np.ones(11)*qin
  if mode==2:
   for j in range(10):
    w=fullR[j]+off[j]*dt;e=fullR[j]+k*off[j]*dt;rn=comp(w,e,ntable[1000]);off[j+1]=(w-rn)/dt;fullR[j]=rn
  oldstate=u.copy();oldb=b;oldsb=sb;oldcum=cum.copy();forced=oracle[ti]*1000
  args=(qin,oldb,oldstate,oldcum,oldsb,initcum,k,Ksec,mult,dt,ntable,yinit,forced,off,mode)
  success=False
  for iteration in range(35):
   f,bnew,unew,cnew,snew,sw,ok=residual(q,h,*args);norm=max(abs(f[0]),abs(f[1]))
   if ok and norm<1e-5:success=True;break
   if not np.all(np.isfinite(f)):break
   dq=max(.01,q*1e-5);dh=1e-5
   fq=residual(q+dq,h,*args)[0];fh=residual(q,h+dh,*args)[0]
   a=(fq[0]-f[0])/dq;bbj=(fh[0]-f[0])/dh;c=(fq[1]-f[1])/dq;d=(fh[1]-f[1])/dh;det=a*d-bbj*c
   if not np.isfinite(det) or abs(det)<1e-12:break
   stepq=(-d*f[0]+bbj*f[1])/det;steph=(c*f[0]-a*f[1])/det
   scale=min(1.,max(10.,.8*q)/max(abs(stepq),1e-12),2./max(abs(steph),1e-12));accepted=False
   for ls in range(14):
    nq=q+scale*stepq;nh=h+scale*steph
    if nq>1 and 3.5<nh<30:
     ff=residual(nq,nh,*args)[0]
     if np.all(np.isfinite(ff)) and max(abs(ff[0]),abs(ff[1]))<norm:
      q=nq;h=nh;accepted=True;break
    scale*=.5
   if not accepted:break
  if not success:
   failure=ti;errout[ti]=norm;iters[ti]=iteration+1;break
  b=bnew;u=unew;cum=cnew;sb=snew;q0=rating(h)
  for j in range(11):
   x=(10-j)*20000.
   result[j,ti]=q0+(cumulative(x,b,u,cum)-cumulative(x,oldb,oldstate,oldcum))/dt
  bb[ti]=b/1000;hh[ti]=h;qb[ti]=q;sbout[ti]=sb;suout[ti]=u.sum();sweep[ti]=sw/dt;iters[ti]=iteration+1;errout[ti]=norm
  balance=(u.sum()+sb-oldstate.sum()-oldsb)/dt-qin+q0;maxmass=max(maxmass,abs(balance))
 return result,bb,hh,qb,sbout,suout,sweep,iters,errout,failure,maxmass
