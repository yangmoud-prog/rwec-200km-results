"""Incremental states in million m3; flow m3/s and time hours."""
import math
from pathlib import Path
import numpy as np
from numba import njit
import base_core as core

@njit(cache=False,nogil=True)
def response(qr,eq,ns,dt_h,tb_h,mode=2,limit=True):
    """Exact linear-forcing relaxation, with state supply constraint."""
    nc,nt=qr.shape;pred=np.zeros_like(qr);bs=np.zeros_like(qr);active=np.zeros_like(qr,dtype=np.int64)
    raw_min=np.full(nc,1e30);conv=1e6/(dt_h*3600)
    a=0.;w0=0.;w1=1.
    if mode==2:
        a=math.exp(-dt_h/tb_h);r=tb_h/dt_h*(-math.expm1(-dt_h/tb_h));w0=r-a;w1=1-r
    for j in range(nc):
        bs[j,0]=eq[j,0];pred[j,0]=qr[j,0]
        for i in range(1,ns[j]):
            b=a*bs[j,i-1]+w0*eq[j,i-1]+w1*eq[j,i]
            rawq=qr[j,i]-(b-bs[j,i-1])*conv;raw_min[j]=min(raw_min[j],rawq)
            if limit:
                top=bs[j,i-1]+qr[j,i]/conv
                if b>top:b=top;active[j,i]=1
                elif b<0:b=0.;active[j,i]=-1
            bs[j,i]=b;pred[j,i]=qr[j,i]-(b-bs[j,i-1])*conv
        pred[j,ns[j]:]=pred[j,ns[j]-1];bs[j,ns[j]:]=bs[j,ns[j]-1]
    return pred,bs,active,raw_min

@njit(cache=False,nogil=True)
def response_gradient(qr,eq,obs,ns,den,dt_h,tb_h):
    pred,b,active,rawmin=response(qr,eq,ns,dt_h,tb_h)
    ge=np.zeros_like(eq);glogt=0.;loss=0.;nc=len(ns);conv=1e6/(dt_h*3600)
    a=math.exp(-dt_h/tb_h);r=tb_h/dt_h*(-math.expm1(-dt_h/tb_h));w0=r-a;w1=1-r
    da=a*dt_h/tb_h;dw0=r-a-da;dw1=-r+a
    for j in range(nc):
        direct=np.zeros(pred.shape[1])
        for i in range(ns[j]):
            dif=pred[j,i]-obs[j,i];loss+=dif*dif/den[j]/nc
            if i:
                gq=2*dif/den[j]/nc;direct[i]-=gq*conv;direct[i-1]+=gq*conv
        future=0.
        for i in range(ns[j]-1,0,-1):
            g=direct[i]+future
            if active[j,i]==0:
                ge[j,i]+=g*w1;ge[j,i-1]+=g*w0
                glogt+=g*(da*b[j,i-1]+dw0*eq[j,i-1]+dw1*eq[j,i])
                future=g*a
            elif active[j,i]==1: future=g
            else:future=0.
        ge[j,0]+=direct[0]+future
    return loss,ge,glogt

def network(theta,features):
    w=theta[:12].reshape(4,3);b=theta[12:16];v=theta[16:20];c=theta[20]
    z=np.tanh(features@w.T+b)
    logit=z@v+c;s=1/(1+np.exp(-np.clip(logit,-40,40)))
    return 2*s,z,s

def neural_loss(theta,features,scale,qr,obs,ns,den,regularization=1e-5):
    coef,z,s=network(theta,features);eq=scale*coef
    loss,ge,gt=response_gradient(qr,eq,obs,ns,den,1.,math.exp(theta[21]))
    g=ge*scale*2*s*(1-s)
    vf=theta[16:20];gz=g[...,None]*vf*(1-z*z)
    grad=np.empty(22);grad[:12]=np.einsum('jth,jtf->hf',gz,features).ravel()
    grad[12:16]=gz.sum(axis=(0,1));grad[16:20]=np.einsum('jt,jth->h',g,z)
    grad[20]=g.sum();grad[21]=gt
    loss+=regularization*np.sum(theta[:20]**2);grad[:20]+=2*regularization*theta[:20]
    return loss,grad

def balance(qin,states,ibars,b,pred,dt_h,tau_h):
    dt=dt_h*3600;raw=.5*(qin[1:]+qin[:-1])
    residual=np.diff(states*1e6+b*1e6)-(ibars[1:]-pred[1:])*dt
    transit=np.r_[qin[0]*tau_h*3600,qin[0]*tau_h*3600+np.cumsum((raw-ibars[1:])*dt)]
    physical=np.diff(states*1e6+b*1e6+transit)-(raw-pred[1:])*dt
    return residual,physical,transit
