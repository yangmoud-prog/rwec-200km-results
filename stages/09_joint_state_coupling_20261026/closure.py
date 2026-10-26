from common import *
from joint import shape,atx,length_scale,rating
from datetime import datetime

def main():
 protocol=read(P/'provenance/frozen_protocol.json');rows=[]
 for e in protocol['order']:
  d=load(P/'data'/f'B_{e}.npz');x=d['x']*1000;sb=[];areas=[];bp=[];profiles=[];dr=[]
  for i,t in enumerate(GRID):
   b,a,c=shape(d['qb_true'][i],d['z'][i,0],1.,d['yinit']);sb.append(atx(c,d['b_true'][i]*1000));areas.append(a);bp.append(b/1000);profiles.append((-500+np.sqrt(500**2+12*a))/6+x*.000075)
   dr.append((d['z'][i,0]-5)*np.exp(-x/length_scale(d['qb_true'][i],1.)))
  sb=np.array(sb);areas=np.array(areas);profiles=np.array(profiles);bp=np.array(bp);dr=np.array(dr)
  sweep=[];sweep_true=[]
  for i in range(1,len(GRID)):
   b0=d['b_true'][i-1]*1000;b1=d['b_true'][i]*1000;xx=np.r_[min(b0,b1),x[(x>min(b0,b1))&(x<max(b0,b1))],max(b0,b1)]
   sweep.append(np.sign(b1-b0)*np.trapz(np.interp(xx,x,(areas[i]+areas[i-1])/2),xx));sweep_true.append(np.sign(b1-b0)*np.trapz(np.interp(xx,x,(d['area'][i]+d['area'][i-1])/2),xx))
  qnet=.5*(d['qb_true'][1:]+d['qb_true'][:-1]-d['q'][1:,0]-d['q'][:-1,0]);res=np.diff(sb)/3600-np.array(sweep)/3600-qnet;res_hec=np.diff(d['sb_true'])/3600-np.array(sweep_true)/3600-qnet
  mask=(x[None,:]<=d['b_true'][:,None]*1000)&(d['b_true'][:,None]>0)
  active=d['b_true']>.1
  row=dict(event=e,subset='train' if e in protocol['train'] else 'validation',shape_stage_RMSE_m=float(np.sqrt(np.mean((profiles[mask]-d['z'][mask])**2))),shape_deltaH_RMSE_m=float(np.sqrt(np.mean((dr-(d['z']-d['zref']))**2))),shape_endpoint_MAE_km=float(np.mean(abs(bp-d['b_true']))),shape_endpoint_active_MAE_km=float(np.mean(abs(bp[active]-d['b_true'][active]))),SB_RMSE_m3=float(np.sqrt(np.mean((sb-d['sb_true'])**2))),SB_NRMSE_range=float(np.sqrt(np.mean((sb-d['sb_true'])**2))/np.ptp(d['sb_true'])),physical_balance_RMSE_m3s=float(np.sqrt(np.mean(res**2))),HEC_hourly_balance_RMSE_m3s=float(np.sqrt(np.mean(res_hec**2))),rating_RMSE_m3s=float(np.sqrt(np.mean((np.array([rating(h) for h in d['z'][:,0]])-d['q'][:,0])**2))))
  rows.append(row);np.savez_compressed(P/'predictions'/f'closure_{e}.npz',time=GRID,b_shape=bp,z_shape=profiles,SB_shape=sb,balance=res,HEC_balance=res_hec)
 frame('known_endpoint_physical_closure.csv',rows)
 save(P/'provenance/closure_completed.json',dict(completed=datetime.now().isoformat(),lambda_multiplier=1,parameters_fitted=False,truth_used='diagnostic only: b,H0,Qb; not autonomous forecasting',gate='reported independent stage/storage/flux residuals; no claim that exact discrete accounting proves physics'))
 print(pd.DataFrame(rows).groupby('subset').mean(numeric_only=True).to_string())
if __name__=='__main__':main()
