"""Training-only diagnostic: six requested sections, excluding outlet loss.

Triggered by training alpha grids; not a replacement of frozen primary target.
Validation targets are never used to choose coefficients or search budgets.
"""
from common import *
from validation import data_A,data_B,weights,fast_base,NODES,COLS
from routing_core import forecast
from scipy.optimize import differential_evolution,minimize,minimize_scalar
import time
def run(label):
    groups=data_A() if label=='A' else data_B();train={e:g for e,g in groups.items() if g['subset']=='train'};lookup=loadnp(OUT/'data'/f'{label}_lookup.npz')
    for g in groups.values():weights(g)
    bounds=[(np.log(.0001),np.log(20)),(np.log(.15),np.log(150)),(0.,8.)]
    def decode(z):return [float(np.exp(z[0])),float(np.exp(z[1])),float(z[2])]
    def loss(z,mode=0,alpha=0.):
        p=decode(z);v=[]
        for g in train.values():
            q=fast_base(g['inflow'],1.,*p) if mode==0 else forecast(g['inflow'],1.,p,alpha,mode,lookup)[0]
            v.append(np.mean(np.sum(g['W'][:,:6]*(q[:,COLS[:6]]-g['mu'][:,:6])**2,axis=0)))
        return float(np.mean(v))
    file=OUT/'fits'/f'{label}_six_target_final.json'
    if file.exists():fits=read(file)
    else:
        opt=differential_evolution(loss,bounds,seed=20261023,popsize=4,maxiter=16,polish=False)
        refine=minimize(loss,opt.x,method='Nelder-Mead',bounds=bounds,options=dict(maxfev=160,xatol=1.e-5,fatol=1.e-7))
        z=refine.x if refine.fun<opt.fun else opt.x;core=decode(z);fits=[dict(model='RWEC-0_six_targets',mode=0,core=core,alpha=0.,loss=loss(z))]
        for mode,name in [(1,'RWEC-L'),(2,'RWEC-V')]:
            opt=minimize_scalar(lambda a:loss(z,mode,a),bounds=(0,1),method='bounded',options=dict(maxiter=20,xatol=5.e-4))
            grid=np.r_[0.,.001,.005,.01,.025,.05,.1,opt.x];sc=np.array([loss(z,mode,float(a)) for a in grid]);alpha=float(grid[np.argmin(sc)])
            x0=np.r_[z,alpha];joint=minimize(lambda x:loss(x[:3],mode,x[3]),x0,method='Nelder-Mead',bounds=bounds+[(0.,1.)],options=dict(maxfev=70,xatol=1.e-5,fatol=1.e-7))
            choose=joint.fun<sc.min();fits.append(dict(model=f'{name}_six_targets',mode=mode,core=decode(joint.x[:3]) if choose else core,alpha=float(joint.x[3]) if choose else alpha,loss=min(float(joint.fun),float(sc.min())),joint_refinement_nfev=int(joint.nfev)))
        save(file,fits)
    rows=[];params=[]
    for f in fits:
        params.append(dict(model=f['model'],k=f['core'][0],n=f['core'][1],tau_per20km_h=f['core'][2],alpha=f['alpha'],train_only=True,posthoc_weight_diagnostic=True))
        for eid,g in groups.items():
            q=forecast(g['inflow'],1.,f['core'],f['alpha'],f['mode'],lookup)[0]
            for c in g['cases']:
                n=len(c['time']);mask=c['time']>=0
                for j,xc in enumerate(TARGETS):rows.append(dict(event=eid,subset=g['subset'],case=c['case'],condition=c['condition'],model=f['model'],xc_km=int(xc),**metrics(c['truth'][mask,j],q[:n][mask,(200-xc)//20],c['time'][mask])))
    frame(f'{label}_six_target_metrics.csv',rows);frame(f'{label}_six_target_parameters.csv',params)
    print(label,'six-section weight diagnostic complete',params,flush=True)
if __name__=='__main__':
    run('A');run('B')
