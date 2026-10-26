"""Post-hoc score alignment to an already-frozen historical predictor; no fitting or selection."""
from common10 import *
from evaluate_A import masks
def main():
 cfg=read(P/'provenance/frozen_protocol.json');rows=[];sources=[];old=P.parent/'05_routing_validation_slope075_20261022/predictions'
 for event in cfg['validation']:
  d=load(P/'data'/f'{event}.npz');f=old/f'A_{event}_RWEC-0.npz';p=load(f);ids=np.rint(d['time']-p['time'][0]).astype(int);assert np.all(ids>=0) and ids.max()<len(p['time']);assert np.max(abs(p['time'][ids]-d['time']))<1e-6;keep=d['time']>=0;qq=p['q'][ids][keep];ys=d['q'][keep];tm=np.arange(len(qq)-1)+.5
  for x in TARGETS:
   col=int(np.argmin(abs(d['x']-x)));y=.5*(ys[:-1,col]+ys[1:,col]);q=.5*(qq[:-1,int((200-x)//20)]+qq[1:,int((200-x)//20)])
   for window,mask in masks(ys[:,col]).items():rows.append(dict(event=event,model='historical_frozen_shared_original',window=window,x_from_outlet_km=int(x),distance_from_inlet_km=200-int(x),**metrics(y[mask],q[mask],tm[mask])))
  sources.append(dict(path=str(f),sha256=sha(f)))
 frame('historical_frozen_original_event_metrics.csv',rows);frame('historical_frozen_original_means.csv',pd.DataFrame(rows).groupby(['window','x_from_outlet_km','distance_from_inlet_km'])[['NSE','RMSE','KGE']].mean().reset_index());save(P/'HISTORICAL_REFERENCE.json',dict(posthoc_score_alignment=True,refitting=False,current_parameters_changed=False,selection_by_validation=False,reference='Stage5 original RWEC shared parameters trained across paired scenarios; separate context from new per-reach dynamic-only fit',sources=sources))
 print(pd.DataFrame(rows).query("window=='nonnegative_hour_endpoints' and x_from_outlet_km in [100,20,0]").groupby('x_from_outlet_km')[['NSE','RMSE']].mean().to_string())
if __name__=='__main__':main()
