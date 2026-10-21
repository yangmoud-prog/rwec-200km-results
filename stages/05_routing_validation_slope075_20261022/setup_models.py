from common import *
import shutil
dirs()
for source,dest in [('prototype.py','routing_core.py'),('hydraulics.py','hydraulic_core.py')]:
    shutil.copyfile(P4/source,OUT/dest)
from hydraulic_core import profile
from scipy.integrate import cumulative_trapezoid
for label,slope in [('A',.00005),('B',.000075)]:
    q=np.arange(500.,12000.1,50.);x=np.arange(0.,200000.1,2000.)
    yr=profile(q,np.full_like(q,5.),slope);yc=profile(q,3.5+1.5*(q/1000)**(2/3),slope)
    y0=yr[np.flatnonzero(q==1000)[0]];a=500*y0+3*y0*y0
    c=np.sqrt(9.81*a/(500+6*y0))-1000/a;assert np.all(c>0)
    delay=cumulative_trapezoid(1/c,x,initial=0)/3600
    np.savez_compressed(OUT/'data'/f'{label}_lookup.npz',q=q,x=x,yref=yr,dz=yc-yr,delay_h=delay)
save(OUT/'provenance/model_source.json',dict(original_complement_source_sha256=sha(P4/'prototype.py'),
    copied_core_sha256=sha(OUT/'routing_core.py'),hydraulic_source_sha256=sha(OUT/'hydraulic_core.py'),no_future_HEC_inputs=True))
print('A/B fixed hydraulic lookups ready',flush=True)
