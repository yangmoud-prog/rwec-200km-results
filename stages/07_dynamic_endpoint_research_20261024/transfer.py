"""Zero-shot A->B test. No B labels used in predictor or parameter fitting."""
from common import *
from predict_endpoint import features,predict
from evaluate import simulate,NODES

meta=read(OUT/'data/metadata.json');tr=[m for m in meta if m['group']=='A' and m['subset']=='train']
pack=[load(OUT/'data'/m['file']) for m in tr]
lib=dict(features=np.array([features(d['inflow'],d['past_H']) for d in pack]),length=np.array([d['true_lb'] for d in pack]))
rows=[];ends=[]
for m in [m for m in meta if m['group']=='B' and m['subset']=='validation']:
 d=load(OUT/'data'/m['file']);central,low,high,ids,w,dist=predict(d['inflow'],d['past_H'],lib);f=GRID>0
 ends.append(dict(event=m['event'],**metric(d['true_lb'][f],central[f]),MAE_km=float(np.mean(abs(d['true_lb'][f]-central[f]))),coverage=float(np.mean((d['true_lb'][f]>=low[f])&(d['true_lb'][f]<=high[f])))))
 for model in ['RWEC-0','RWEC-L-predicted']:
  fit=read(OUT/'fits'/f'A_{model}.json');q,quality=simulate(d,model,fit,central)
  for j,xc in enumerate(TARGETS):
   mask=f&(d['truth'][:,j]>=0)
   rows.append(dict(event=m['event'],xc_km=int(xc),model=model,**metric(d['truth'][mask,j],q[mask,NODES[j]])))
frame('zero_shot_A_to_B_flow.csv',rows);frame('zero_shot_A_to_B_endpoint.csv',ends)
save(OUT/'provenance/transfer_protocol.json',dict(training_domain='A49 only',target='B6 validation only',B_labels_used_in_fitting=False,geometry_correction='none; raw-transfer stress test',known_boundary_domain_shift=True))
print(pd.DataFrame(rows).groupby('model')[['NSE','RMSE']].mean().to_string());print('Endpoint MAE',np.mean([r['MAE_km'] for r in ends]))
