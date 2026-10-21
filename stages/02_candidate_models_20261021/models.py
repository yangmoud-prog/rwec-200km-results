"""Document-specified lagged complementary routing, in million cubic metres."""
import math
import numpy as np
from numba import njit

@njit(cache=True, nogil=True)
def storage(sw, se, n):
    if sw<=0 or se<=0 or n<=0:
        return 0.0
    lo=min(sw,se);hi=max(sw,se)
    return lo / math.exp(math.log1p(math.exp(n*math.log(lo/hi)))/n)

@njit(cache=True, nogil=True)
def initial_ratio(k,n):
    lo=0.0;hi=max(1.0,k)
    while (hi/(hi+1))**n+(hi/(hi+k))**n<1:
        hi*=2
    for _ in range(80):
        m=(lo+hi)/2
        if (m/(m+1))**n+(m/(m+k))**n>1:
            hi=m
        else:
            lo=m
    return (lo+hi)/2

@njit(cache=True, nogil=True)
def sample(q,idx,N):
    if idx<=0:return q[0]
    if idx>=N-1:return q[N-1]
    j=int(math.floor(idx));f=idx-j
    return q[j]*(1-f)+q[j+1]*f

@njit(cache=True, nogil=True)
def length_factor(L,theta):
    if theta<1e-10:return L
    ell=200000/theta
    return -ell*math.expm1(-L/ell)

@njit(cache=True, nogil=True)
def route(q,h,N,k,n,tau,alpha,theta,L,dt_h=1.0):
    out=np.empty(N);states=np.empty(N);ibars=np.empty(N);residual=np.zeros(N)
    dt=dt_h*3600.0
    prev=initial_ratio(k,n)*max(0,q[0])*dt/1e6
    out[0]=max(0,q[0]);states[0]=prev;ibars[0]=max(0,q[0])
    lf=length_factor(L,theta)
    for i in range(1,N):
        raw=.5*(sample(q,i-tau/dt_h,N)+sample(q,i-1-tau/dt_h,N))
        ib=max(0.0,raw);ibars[i]=ib
        vin=ib*dt/1e6
        sw=max(1e-12,prev+vin)
        se=max(1e-12,prev+k*vin+alpha*500*lf*max(h[i]-5,0)/1e6)
        s=storage(sw,se,n)
        qo=max(0.0,(sw-s)*1e6/dt)
        residual[i]=(s-prev)*1e6-(ib-qo)*dt
        out[i]=qo;states[i]=s;prev=s
    return out,states,ibars,residual

@njit(cache=True, nogil=True)
def loss(q,h,obs,ns,denom,indices,k,n,tau,alpha,theta,L):
    # No target data enter initialization, inflow sampling or state propagation.
    result=0.0;dt=3600.0;ir=initial_ratio(k,n);lf=length_factor(L,theta)
    for jj in range(len(indices)):
        j=indices[jj];N=ns[j]
        prev=ir*max(0,q[j,0])*dt/1e6
        sse=(max(0,q[j,0])-obs[j,0])**2
        for i in range(1,N):
            ib=max(0.0,.5*(sample(q[j],i-tau,N)+sample(q[j],i-1-tau,N)))
            vin=ib*dt/1e6;sw=max(1e-12,prev+vin)
            se=max(1e-12,prev+k*vin+alpha*500*lf*max(h[j,i]-5,0)/1e6)
            s=storage(sw,se,n);qo=max(0.0,(sw-s)*1e6/dt)
            sse+=(qo-obs[j,i])**2;prev=s
        result+=sse/denom[j]
    return result/len(indices)

def decode(x):
    return (math.exp(math.log(1e-6)+float(x[0])*math.log(1e8)),
            math.exp(math.log(.01)+float(x[1])*math.log(30000)),float(x[2])*150)

def encode(k,n,tau):
    return np.array([(math.log(k)-math.log(1e-6))/math.log(1e8),
        (math.log(n)-math.log(.01))/math.log(30000),tau/150])

def decode_shared(x,model):
    alpha=float(x[0])*10
    theta=0.0 if model==1 else math.expm1(float(x[1])*math.log(201))
    return alpha,theta
