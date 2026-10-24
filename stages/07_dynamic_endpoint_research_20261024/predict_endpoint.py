"""Fixed-origin forecaster. Only forecast I and H at/before origin enter query."""
from common import *

def features(inflow,past_H):
 assert len(inflow)==289 and len(past_H)==49
 # Fixed physical scales, fixed forecast horizon; no future outlet data.
 return np.r_[inflow[48::12]/5000., np.max(inflow[48:])/5000.,
              np.argmax(inflow[48:])/240.,np.trapz(inflow[48:]-1000)/1.e6,
              past_H[::6]/4.]

def predict(inflow,past_H,library):
 query=features(inflow,past_H)
 distances=np.mean((library['features']-query)**2,axis=1)
 ids=np.argsort(distances,kind='stable')[:min(8,len(distances))]
 w=1/(np.sqrt(distances[ids])+.02);w/=w.sum()
 trajectories=library['length'][ids]
 center=w@trajectories
 low,high=np.quantile(trajectories,[.1,.9],axis=0)
 return center,low,high,ids,w,distances[ids]

def main():
 dirs();meta=read(OUT/'data/metadata.json');rows=[];neighbors=[]
 for group in ['A','B']:
  tr=[m for m in meta if m['group']==group and m['subset']=='train']
  packets={m['event']:load(OUT/'data'/m['file']) for m in tr}
  for m in [x for x in meta if x['group']==group]:
   d=load(OUT/'data'/m['file']);available=[a for a in tr if a['event']!=m['event']]
   lib=dict(features=np.array([features(packets[a['event']]['inflow'],packets[a['event']]['past_H']) for a in available]),
            length=np.array([packets[a['event']]['true_lb'] for a in available]))
   central,low,high,ids,w,dist=predict(d['inflow'],d['past_H'],lib)
   np.savez_compressed(OUT/'predictions'/f"endpoint_{group}_{m['event']}.npz",time=GRID,central=central,low=low,high=high)
   for idx,weight,dd in zip(ids,w,dist):neighbors.append(dict(group=group,event=m['event'],subset=m['subset'],neighbor=available[idx]['event'],weight=weight,distance=dd))
   f=GRID>0;y=d['true_lb'][f]
   rows.append(dict(group=group,event=m['event'],subset=m['subset'],**metric(y,central[f]),
    MAE_km=float(np.mean(abs(y-central[f]))),coverage=float(np.mean((y>=low[f])&(y<=high[f]))),
    width_km=float(np.mean(high[f]-low[f])),history_H_range_m=float(np.ptp(d['past_H'])),
    peak_length_error_km=float(central[f].max()-y.max())))
 frame('endpoint_prediction_metrics.csv',rows);frame('endpoint_neighbors.csv',neighbors)
 save(OUT/'provenance/endpoint_prediction_inputs.json',dict(origin_h=0,history_h=[-48,0],forecast_h=[0,240],allowed_query=['upstream_flow_forecast','outlet_water_level_history_through_origin'],neighbors=8,training_predictions='event_leave_one_out',validation_predictions='training_library_only',future_truth_input=False))
 print(pd.DataFrame(rows).groupby(['group','subset'])[['RMSE','MAE_km','coverage','width_km','history_H_range_m']].mean().to_string())

if __name__=='__main__':main()
