"""Geometry-only quasi-steady lookup. Contains no HEC result dependency."""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from common import OUT,save

def gradient(y,q,slope):
    a=500*y+3*y*y;t=500+6*y;p=500+2*y*np.sqrt(10)
    fr2=q*q*t/(9.81*a**3)
    sf=.025**2*q*q/(a*a*(a/p)**(4/3))
    if np.any(y<=0) or np.any(fr2>=.95):raise ValueError('GVF left subcritical wet domain')
    return (sf-slope)/(1-fr2)

def profile(q,h,slope=.0001):
    x=np.arange(0.,200000.1,500.)
    y=np.empty((len(q),len(x)));y[:,0]=h
    for j in range(len(x)-1):
        z=y[:,j];dx=500.
        a=gradient(z,q,slope);b=gradient(z+dx*a/2,q,slope)
        c=gradient(z+dx*b/2,q,slope);d=gradient(z+dx*c,q,slope)
        y[:,j+1]=z+dx*(a+2*b+2*c+d)/6
    return y[:,::4]

def build():
    q=np.arange(500.,12000.1,50.);x=np.arange(0.,200000.1,2000.)
    yc=profile(q,3.5+1.5*(q/1000)**(2/3));yr=profile(q,np.full_like(q,5.))
    initial=yr[np.flatnonzero(q==1000)[0]]
    area=500*initial+3*initial**2;top=500+6*initial
    c=np.sqrt(9.81*area/top)-1000/area
    assert np.all(c>0)
    delay=cumulative_trapezoid(1/c,x,initial=0)/3600
    np.savez_compressed(OUT/'data/hydraulic_lookup.npz',q=q,x=x,yref=yr,dz=yc-yr,delay_h=delay)
    save(OUT/'provenance/hydraulic_lookup.json',dict(q_min=500,q_max=12000,q_step=50,dx_integration_m=500,dx_saved_m=2000,
         source='fixed channel geometry, control relation and initial Q1000/H5 only',slope=.0001,
         delay_200km_h=float(delay[-1]),max_reference_fr2=float(np.max(q[:,None]**2*(500+6*yr)/(9.81*(500*yr+3*yr**2)**3))),
         future_HEC_results_used=False))
    return q,x,yr,yc-yr,delay
if __name__=='__main__':build()
