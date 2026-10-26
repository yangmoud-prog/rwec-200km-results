from common import *

def feature(i,hpast):
 assert len(i)==553 and len(hpast)==49
 return np.r_[i[::12]/5000.,i.max()/5000.,np.argmax(i)/552.,np.trapz(i-1000)/1.e6,hpast[::6]/4.]

def forecast(i,hpast,library_features,library_b):
 d=np.mean((library_features-feature(i,hpast))**2,axis=1)
 ids=np.argsort(d,kind='stable')[:min(8,len(d))];w=1/(np.sqrt(d[ids])+.02);w/=w.sum()
 return w@library_b[ids],*np.quantile(library_b[ids],[.1,.9],axis=0),ids,w

def main():
 plan=read(P/'provenance/C_frozen.json');rows=[];packets={};quality=[]
 for eid in plan['order']:
  dv=load(C/'results'/f'C_{eid}_var/full_fields.npz');dr=load(C/'results'/f'C_{eid}_ref/full_fields.npz')
  x=dv['x_m'][::-1]/1000;tt=dv['time'];use=tt>=-1.e-6;assert len(tt[use])==553
  z=dv['z'][use,::-1];zr=dr['z'][use,::-1];q=dv['q'][use,::-1];qr=dr['q'][use,::-1];area=dv['area'][use,::-1]
  endpoints=np.array([extent(h,x) for h in z-zr]);b=endpoints[:,0]
  qb=np.array([np.interp(bb,x,qq) for bb,qq in zip(b,q)]);ab=np.array([np.interp(bb,x,aa) for bb,aa in zip(b,area)])
  cumul=cumulative_trapezoid(area,x*1000,axis=1,initial=0);sm=np.array([v[-1]-np.interp(bb,x,v) for bb,v in zip(b,cumul)])
  init_area=area[0];V0=[];y0=[]
  for j in range(10):
   lo=(9-j)*20;hi=(10-j)*20;ix=(x>=lo)&(x<=hi)
   V0.append(np.trapz(init_area[ix],x[ix]*1000));y0.append(np.mean((z[0]-x*1000*plan['selected_slope'])[ix]))
  inflow=hec.boundary(plan['events'][eid],GRID,'var')[0]
  hpast=np.interp(np.arange(-48,1),tt,dv['z'][:,-1]);assert np.ptp(hpast)<1.e-5
  targetcols=np.array([np.argmin(abs(x-s)) for s in TARGETS]);targetq=q[:,targetcols].T
  packet=dict(time=GRID,x=x,q=q,qref=qr,z=z,zref=zr,area=area,b_true=b,b_connected=endpoints[:,1],islands=endpoints[:,2],censored=endpoints[:,3],qb_true=qb,ab_true=ab,sm_true=sm,
   inflow=inflow,past_H=hpast,truth=targetq,initial_volumes=np.array(V0),initial_depth=np.array(y0),features=feature(inflow,hpast))
  packets[eid]=packet
  subset='train' if eid in plan['train'] else 'validation'
  row=dict(event=eid,subset=subset,minQ=float(q.min()),negative_count=int((q<0).sum()),max_b_km=float(b.max()),detached_hours=int(np.sum((b>0)&(endpoints[:,1]==0))),multiple_component_hours=int(np.sum(endpoints[:,2]>1)),censored_hours=int(np.sum(endpoints[:,3])),past_H_range=float(np.ptp(hpast)))
  for kind in ['ref','var']:
   qual=read(C/'results'/f'C_{eid}_{kind}/quality.json');quality.append(dict(event=eid,kind=kind,**{k:qual[k] for k in ['basic_quality_pass','initialization_pass','tail_recovered','min_flow_m3s','initial_flow_error_m3s','last24h_flow_error_m3s','final_balance_residual_pct_total_inflow']}))
  rows.append(row)
 neighbors=[];predrows=[]
 for eid,packet in packets.items():
  names=[e for e in plan['train'] if e!=eid]
  bp,bl,bh,ids,w=forecast(packet['inflow'],packet['past_H'],np.array([packets[e]['features'] for e in names]),np.array([packets[e]['b_true'] for e in names]))
  packet.update(b_pred=bp,b_low=bl,b_high=bh);np.savez_compressed(P/'data'/f'{eid}.npz',**packet)
  predrows.append(dict(event=eid,subset='train' if eid in plan['train'] else 'validation',MAE_km=float(np.mean(abs(bp-packet['b_true']))),RMSE_km=float(np.sqrt(np.mean((bp-packet['b_true'])**2))),coverage=float(np.mean((packet['b_true']>=bl)&(packet['b_true']<=bh))),mean_width_km=float(np.mean(bh-bl))))
  for idx,weight in zip(ids,w):neighbors.append(dict(event=eid,neighbor=names[idx],weight=float(weight)))
 frame('event_physics.csv',rows);frame('HEC_quality.csv',quality);frame('endpoint_predictor_metrics.csv',predrows);frame('endpoint_neighbors.csv',neighbors)
 save(P/'provenance/predictor_frozen.json',dict(method='Stage7 fixed-origin8-neighbor predictor, fixed horizon extended to552h',features='upstream forecast samples and statistics; stage history[-48,0]h',future_stage_used=False,future_HEC_state_used=False,train='leave-one-event-out14 training events',validation='14 training events only',sealed_test_count=0))
 print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(predrows).groupby('subset').mean(numeric_only=True).to_string())

if __name__=='__main__':main()
