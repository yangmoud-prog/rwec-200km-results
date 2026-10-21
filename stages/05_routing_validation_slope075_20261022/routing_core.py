"""Causal RWEC cascade. Forecast function has no HEC target input.

The oracle comparator is an explicitly separate wrapper, carrying true future
boundary level through oracle_h. This input is forbidden for causal modes.
"""
import numpy as np
from numba import njit

@njit(cache=True)
def complementary(a,b,n):
    if a<=0 or b<=0:return 0.
    lo=min(a,b);hi=max(a,b)
    return lo/np.exp(np.log1p(np.exp(n*np.log(lo/hi)))/n)

@njit(cache=True)
def initial_ratio(k,n):
    lo=0.;hi=1.e8
    for _ in range(100):
        r=(lo+hi)/2
        if (r/(r+1))**n+(r/(r+k))**n>1:hi=r
        else:lo=r
    return (lo+hi)/2

@njit(cache=True)
def interp_history(series,i,z,current):
    if z<=0:return series[0]
    j=int(np.floor(z));w=z-j
    if j>=i:return current
    a=series[j];b=current if j+1>=i else series[j+1]
    return a+(b-a)*w

@njit(cache=True)
def table_value(tab,qq,j):
    pos=(qq-500)/50
    k=int(np.floor(pos));w=pos-k
    if k<0:return tab[0,j]
    if k>=tab.shape[0]-1:return tab[-1,j]
    return tab[k,j]+w*(tab[k+1,j]-tab[k,j])

@njit(cache=True)
def surfaces(i,dt,qhist,guess,yref,dz,delay):
    d=np.empty(101);ar=np.empty(101)
    for j in range(101):
        past=interp_history(qhist,i,i-delay[j]/dt,guess)
        d[j]=table_value(dz,past,j)
        y=table_value(yref,guess,j)
        ar[j]=(500+6*y)*d[j]+3*d[j]*d[j]
    lb=0.
    for j in range(101):
        if d[j]>.1:
            lb=j*2000.
            if j<100 and d[j+1]<=.1:lb+=2000*(d[j]-.1)/(d[j]-d[j+1])
    return d,ar,lb

@njit(cache=True)
def engine(inflow,dt,k,n,tau,alpha,mode,yref,dz,delay,oracle_h):
    # nodes: 200,180,...,0 km. One global triplet per20km segment.
    nt=len(inflow);q=np.full((nt,11),inflow[0]);q[:,0]=inflow
    s=np.full(10,initial_ratio(k,n)*inflow[0]*dt*3600)
    h=np.full(nt,5.);lbs=np.zeros(nt);vol=np.zeros((nt,101));delta=np.zeros((nt,101))
    iterations=np.zeros(nt,np.int32);failed=np.zeros(nt,np.int32)
    mass=np.zeros((nt,10));transit=np.zeros(10);fullmass=np.zeros((nt,10));sinit=s.copy()
    sumin=np.zeros(10);sumout=np.zeros(10)
    qdomain=0;history=q[:,10].copy();maxfixedres=0.
    for i in range(1,nt):
        guess=q[i-1,10];oldh=3.5+1.5*(guess/1000)**(2/3);oldlb=lbs[i-1]
        niter=1 if mode==0 or mode==4 else 80
        ss=np.empty(10);candidate=np.empty(11);candidate[0]=inflow[i]
        lag=np.empty(10);raw=np.empty(10);storedlast=np.empty(10)
        converged=False
        for it in range(niter):
            dh=0.;ar=np.zeros(101);d=np.zeros(101);lb=0.
            if mode in (1,2,3):
                d,ar,lb=surfaces(i,dt,history,guess,yref,dz,delay)
                hh=3.5+1.5*(guess/1000)**(2/3);dh=max(0.,hh-5.)
            elif mode==4:hh=oracle_h[i];dh=max(0.,hh-5.)
            else:hh=3.5+1.5*(guess/1000)**(2/3)
            for j in range(10):
                u=candidate[j]
                a=interp_history(q[:,j],i,i-tau/dt,u)
                b=interp_history(q[:,j],i,i-1-tau/dt,u)
                lag[j]=.5*(a+b)*dt*3600
                raw[j]=.5*(q[i-1,j]+u)*dt*3600
                extra=0.
                xd=(9-j)*20000.;xu=xd+20000
                if mode==1:extra=alpha*500*dh*max(0.,min(lb,xu)-xd)
                elif mode==2:
                    jl=(9-j)*10
                    for m in range(jl,jl+10):extra+=alpha*1000*(ar[m]+ar[m+1])
                elif mode==3 or mode==4:extra=alpha*500*dh*20000
                sw=s[j]+lag[j];se=s[j]+k*lag[j]+extra
                ss[j]=complementary(sw,se,n)
                candidate[j+1]=(sw-ss[j])/(dt*3600)
            new=candidate[10]
            residual=abs(new-guess)
            if niter==1 or (residual<1.e-4 and abs(hh-oldh)<1.e-6 and abs(lb-oldlb)<1.):
                converged=True;break
            oldh=hh;oldlb=lb;guess=.5*(guess+new)
        iterations[i]=it+1;failed[i]=0 if converged else 1
        if niter>1:maxfixedres=max(maxfixedres,residual)
        q[i]=candidate;history[i]=candidate[10]
        h[i]=3.5+1.5*(candidate[10]/1000)**(2/3)
        lbs[i]=lb;vol[i]=ar;delta[i]=d
        if candidate[10]<500 or candidate[10]>12000:qdomain+=1
        for j in range(10):
            mass[i,j]=ss[j]-s[j]-lag[j]+candidate[j+1]*dt*3600
            transit[j]+=raw[j]-lag[j]
            sumin[j]+=raw[j];sumout[j]+=candidate[j+1]*dt*3600
            fullmass[i,j]=ss[j]-sinit[j]+transit[j]-sumin[j]+sumout[j]
        s=ss.copy()
    return q,h,lbs,delta,vol,iterations,failed,mass,fullmass,qdomain,maxfixedres

def forecast(inflow,dt,params,alpha,mode,lookup):
    if mode not in (0,1,2,3):raise ValueError('Only causal modes allowed')
    return engine(np.asarray(inflow),dt,*params,alpha,mode,lookup['yref'],lookup['dz'],lookup['delay_h'],np.zeros(len(inflow)))

def oracle(inflow,dt,params,alpha,future_boundary_h,lookup):
    return engine(np.asarray(inflow),dt,*params,alpha,4,lookup['yref'],lookup['dz'],lookup['delay_h'],np.asarray(future_boundary_h))
