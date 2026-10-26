import numpy as np

def fractions(b):
 return np.clip((np.arange(200,0,-20)[:,None]-np.asarray(b)[None,:])/20,0,1)

def sample_flow(q,b):
 # q: upstream->downstream faces; b: km upstream of outlet.
 return np.array([np.interp(200-bb,np.arange(0,201,20),q[:,i]) for i,bb in enumerate(b)])

def audit(q,vol,b,dt_h):
 f=fractions(b);sm=np.sum(f*vol,axis=0);sec=dt_h*3600
 dq=.5*(q[:,1:]+q[:,:-1]);fb=.5*(f[:,1:]+f[:,:-1]);vb=.5*(vol[:,1:]+vol[:,:-1])
 sweep=-np.sum(vb*np.diff(f,axis=1),axis=0)
 # Exact bilinear space/time quadrature; separate from instantaneous physical Q(b,t).
 qquad=dq[0]-np.sum(fb*(dq[:-1]-dq[1:]),axis=0)
 residual=np.diff(sm)-(dq[0]-qquad)*sec+sweep
 physical=sample_flow(q,b);relative=qquad+sweep/sec
 sampled_residual=np.diff(sm)-(dq[0]-.5*(physical[:-1]+physical[1:]))*sec+sweep
 no_sweep=np.diff(sm)-(dq[0]-qquad)*sec
 return dict(sm=sm,qb_physical=physical,qb_quadrature=qquad,relative_flux=relative,swept_volume=sweep,
  balance_residual=residual,sampled_balance_residual=sampled_residual,omitting_sweep_residual=no_sweep,
  max_balance_m3=float(np.max(abs(residual))),rms_balance_m3s=float(np.sqrt(np.mean((residual/sec)**2))),
  rms_sampled_balance_m3s=float(np.sqrt(np.mean((sampled_residual/sec)**2))),rms_no_sweep_m3s=float(np.sqrt(np.mean((no_sweep/sec)**2))),
  jump_steps_over20km=int(np.sum(abs(np.diff(b))>20)))
