from common import *
from datetime import datetime

def main():
 plan=read(B/'config/pilot_frozen.json'); packets={}; records=[]; sources=[]
 protocol=dict(created=datetime.now().isoformat(),stage=9,train=plan['train'],validation=plan['validation'],order=plan['order'],scope='B joint state mechanism and forecast; C separate identifiability',future_HEC_inputs_allowed_only_for_oracle=True,sealed_test_opened=0,new_HEC_runs=0,window_h=[0,288],targets_outlet_km=TARGETS.tolist(),shape='initial profile + normal depth response + exponential downstream head response',shape_threshold_m=.1,dt_candidates_h=[.5,1.],max_objective_calls_per_model_per_dt=320,optimizer='48 seeded Latin-hypercube points; bounded Nelder-Mead restarts; training-only selection',new_model_parameters=['k','K20_h','lambda_multiplier'],feedback_off='same fitted parameters; upstream flux independent of current or historical backwater state',primary_objective='mean NSE over14trainingevents and7fixedsections',failures='retained, no clipping and no replacement',validation_selection=False)
 save(P/'provenance/frozen_protocol.json',protocol)
 for e in plan['order']:
  src=A/'05_routing_validation_slope075_20261022/data'/f'B_{e}.npz';d=load(src);use=(d['time']>=0)&(d['time']<=288);x=d['x'];q=d['q_control'][use];z=d['z_control'][use];zr=d['z_ref'][use];y=z-x[None,:]*.075;area=500*y+3*y*y
  ends=np.array([extent(h,x) for h in z-zr]);b=ends[:,0];cum=cumulative_trapezoid(area,x*1000,axis=1,initial=0)
  sb=np.array([np.interp(v,x,row) for v,row in zip(b,cum)]);qb=np.array([np.interp(v,x,row) for v,row in zip(b,q)])
  initial_volumes=np.array([np.interp(200-j*20,x,cum[0])-np.interp(180-j*20,x,cum[0]) for j in range(10)])
  initial_depth=np.array([y[0,(x>=180-j*20)&(x<=200-j*20)].mean() for j in range(10)])
  i=d['inflow'][use];hp=np.interp(np.arange(-48,1),d['time'],d['z_control'][:,0]);features=np.r_[i[::12]/5000,i.max()/5000,np.argmax(i)/288,np.trapz(i-1000)/1e6,hp[::6]/4]
  packets[e]=dict(time=GRID,x=x,inflow=i,past_H=hp,features=features,q=q,qref=d['q_ref'][use],z=z,zref=zr,area=area,b_true=b,b_connected=ends[:,1],qb_true=qb,sb_true=sb,sm_true=cum[:,-1]-sb,ab_true=np.array([np.interp(v,x,row) for v,row in zip(b,area)]),truth=q[:,(TARGETS/2).astype(int)].T,yinit=y[0],initial_volumes=initial_volumes,initial_depth=initial_depth)
  records.append(dict(event=e,subset='train' if e in plan['train'] else 'validation',minQ=float(q.min()),final_Q_error=float(np.max(abs(q[-1]-1000))),max_b_km=float(b.max()),detached_hours=int(np.sum((ends[:,0]>0)&(ends[:,1]==0))),island_hours=int(np.sum(ends[:,2]>1)),max_stage=float(z[:,0].max())))
  sources.append(dict(path=str(src.relative_to(W)),sha256=sha(src)))
 for e,d in packets.items():
  names=[v for v in plan['train'] if v!=e];dist=np.mean((np.array([packets[v]['features'] for v in names])-d['features'])**2,axis=1);ids=np.argsort(dist,kind='stable')[:8];w=1/(np.sqrt(dist[ids])+.02);w/=w.sum();bs=np.array([packets[names[j]]['b_true'] for j in ids]);d.update(b_pred=w@bs,b_low=np.quantile(bs,.1,axis=0),b_high=np.quantile(bs,.9,axis=0));np.savez_compressed(P/'data'/f'B_{e}.npz',**d)
 frame('event_physics.csv',records);save(P/'provenance/source_manifest.json',dict(sources=sources,source_plan_sha256=sha(B/'config/pilot_frozen.json'),sealed_test_opened=0,interpolation='hourly source, 2km fields; trapezoidal geometry area'))
 print(pd.DataFrame(records).to_string(index=False))
if __name__=='__main__':main()
