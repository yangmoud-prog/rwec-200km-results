from pathlib import Path
import sys,json,hashlib,types,importlib.util,datetime
import numpy as np
import pandas as pd

OUT=Path(__file__).resolve().parent
A=OUT.parents[1];P5=A/'analysis/05_routing_validation_slope075_20261022';P3=A/'analysis/03_internal_discovery_20261021';P4=A/'analysis/04_autonomous_backwater_20261021'
for folder in ['tables','figures','provenance','data']:(OUT/folder).mkdir(exist_ok=True)
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,v):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def npz(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def frame(name,rows):pd.DataFrame(rows).to_csv(OUT/'tables'/name,index=False,encoding='utf-8-sig')

# Execute the unchanged kernel as ordinary Python when the local accelerator is absent.
# Neither package installation nor edits to previous stages are required.
stub=types.ModuleType('numba');stub.njit=lambda *a,**kw: (a[0] if a and callable(a[0]) else lambda fn:fn)
sys.modules['numba']=stub
spec=importlib.util.spec_from_file_location('frozen_routing_core',P5/'routing_core.py');kernel=importlib.util.module_from_spec(spec);spec.loader.exec_module(kernel)
split=read(P4/'provenance/event_split.json');events=split['validation'];assert len(events)==21
plan=read(A/'config/formal_plan.json');cases={}
for c in plan['formal_cases']:
 if c['event_id'] in events and c['role']=='main':cases.setdefault(c['event_id'],{})[c['kind']]=c
lookup=npz(P5/'data/A_lookup.npz')
fits={n:read(P5/'fits'/f) for n,f in [('linear_fixed','A_mode4_final.json'),('linear_joint','A_joint_mode4_final.json')]}
meta3=read(P3/'data/metadata.json');d3=npz(P3/'data/all.npz');j3={m['case']:j for j,m in enumerate(meta3)}
old={n:npz(P3/'predictions'/f'{file}.npz')['q'][9] for n,file in [('legacy0','RWEC0'),('legacy1','RWEC1'),('legacy0_balanced','RWEC0J')]}
names=['RWEC0','linear_fixed','linear_joint'];X=[20,40,60,80,100,120]
metric_rows=[];decomp=[];response=[];phases=[];source=[];checks=[];gallery={}

def metric(o,p,t):
 n=len(o)
 if n<3:return dict(n=n,NSE=np.nan,RMSE=np.nan,KGE=np.nan,bias=np.nan,peak_error_pct=np.nan,peak_time_h=np.nan)
 dif=p-o;den=np.sum((o-o.mean())**2);sd=o.std()
 corr=np.corrcoef(p,o)[0,1] if p.std()>1e-12 and sd>1e-12 else np.nan
 kge=1-np.sqrt((corr-1)**2+(p.std()/sd-1)**2+(p.mean()/o.mean()-1)**2) if sd>1e-12 and abs(o.mean())>1e-12 else np.nan
 return dict(n=n,NSE=1-np.sum(dif*dif)/den if den>1e-12 else np.nan,RMSE=np.sqrt(np.mean(dif*dif)),KGE=kge,bias=dif.mean(),peak_error_pct=(p.max()/o.max()-1)*100 if abs(o.max())>1e-12 else np.nan,peak_time_h=t[np.argmax(p)]-t[np.argmax(o)])
def windows(t,o):
 f=t>=0;neg=np.flatnonzero(o<0);last=float(t[neg[-1]]) if len(neg) else None
 return {'all':f,'nonnegative':f&(o>=0),'post_reverse':f if last is None else f&(t>last),'post_reverse12':f if last is None else f&(t>last+12)},last
def phase(t,e):
 return np.select([t<e['tb_h'],t<e['tb_h']+e['Tr_h'],t<e['tb_h']+e['Tr_h']+e['Th_h'],t<e['Hend_h']],['before','rise','hold','fall'],default='recovered')

for eid in events:
 ds={}
 for kind in ['var','ref']:
  c=cases[eid][kind];p=A/'results'/c['case']/'reach_timeseries.npz';d=npz(p);ds[kind]=d
  source.append(dict(path=str(p.relative_to(A)),sha256=sha(p)))
  qual=read(p.parent/'quality.json')
  checks.append(dict(event=eid,kind=kind,initialization_pass=qual['initialization_pass'],initial_flow_error_m3s=qual['initial_flow_error_m3s'],pre_boundary_flow_error_m3s=qual['pre_boundary_flow_error_m3s']))
 d=ds['var'];t=np.round(d['time_since_flood_h'],6);n=len(t);assert np.allclose(t,np.round(ds['ref']['time_since_flood_h'],6))
 saved=npz(P5/'predictions'/f'A_{eid}_RWEC-0.npz');assert np.allclose(saved['time'][:n],t)
 inf=saved['q'][:n,0];e=cases[eid]['var']['event']
 # Replay the originally imposed analytical boundary, matching stage 05 exactly.
 hh=np.full(n,5.);tb=e['tb_h'];tr=e['Tr_h'];th=e['Th_h'];tf=e['Tf_h'];dh=e['deltaH_m']
 rising=(t>=tb)&(t<tb+tr);holding=(t>=tb+tr)&(t<tb+tr+th);falling=(t>=tb+tr+th)&(t<tb+tr+th+tf)
 hh[rising]+=dh*(1-np.cos(np.pi*(t[rising]-tb)/tr))/2;hh[holding]+=dh;hh[falling]+=dh*(1+np.cos(np.pi*(t[falling]-tb-tr-th)/tf))/2
 pred={'RWEC0':saved['q'][:n]}
 pred_ref={'RWEC0':pred['RWEC0']}
 for name,f in fits.items():
  pred[name]=kernel.oracle(inf,1.,f['core'],f['alpha'],hh,lookup)[0]
  pred_ref[name]=kernel.oracle(inf,1.,f['core'],f['alpha'],np.full(n,5.),lookup)[0]
 for xc in X:
  o=d['flow_m3s'][:,xc//10];oref=ds['ref']['flow_m3s'][:,xc//10];j=(200-xc)//20
  ws,last=windows(t,o);f=ws['all'];nonneg=ws['nonnegative'];pp=pred['RWEC0'][:,j];err=pp-o;se=err**2
  decomp.append(dict(event=eid,xc_km=xc,total_hours=int(f.sum()),negative_hours=int((f&(o<0)).sum()),negative_fraction=float((f&(o<0)).sum()/f.sum()),min_q=float(o.min()),last_negative_h=last,
   SSE_all=float(se[f].sum()),SSE_negative=float(se[f&(o<0)].sum()),SSE_nonnegative=float(se[nonneg].sum()),negative_SSE_fraction=float(se[f&(o<0)].sum()/se[f].sum()),
   input_flood_peak=float(inf.max()),truth_flood_sd=float(o[f].std()),paired_delta_RMS=float(np.sqrt(np.mean((o[f]-oref[f])**2))),paired_delta_peak_abs=float(abs(o[f]-oref[f]).max())))
  for name in names:
   for wn,mask in ws.items():metric_rows.append(dict(event=eid,xc_km=xc,model=name,window=wn,**metric(o[mask],pred[name][mask,j],t[mask])))
   metric_rows.append(dict(event=eid,xc_km=xc,model=name,window='reference_all',**metric(oref[f],pred_ref[name][f,j],t[f])))
   common=nonneg&(oref>=0)
   metric_rows.append(dict(event=eid,xc_km=xc,model=name,window='reference_common_nonnegative',**metric(oref[common],pred_ref[name][common,j],t[common])))
   m=nonneg&(oref>=0);actual=(o-oref)[m];estimate=(pred[name][:,j]-pred_ref[name][:,j])[m]
   response.append(dict(event=eid,xc_km=xc,model=name,**metric(actual,estimate,t[m]),actual_peak_abs=float(abs(actual).max()),predicted_peak_abs=float(abs(estimate).max())))
   phase_names=phase(t,e)
   for ph in np.unique(phase_names):
    m=nonneg&(phase_names==ph)
    if m.sum():phases.append(dict(event=eid,xc_km=xc,model=name,phase=ph,n=int(m.sum()),RMSE=float(np.sqrt(np.mean((pred[name][m,j]-o[m])**2))),bias=float((pred[name][m,j]-o[m]).mean())))
  if xc==100:
   for name,v in old.items():
    k=j3[cases[eid]['var']['case']];p=v[k,:n]
    for wn,mask in ws.items():metric_rows.append(dict(event=eid,xc_km=xc,model=name,window=wn,**metric(o[mask],p[mask],t[mask])))
    m=nonneg&(oref>=0);actual=(o-oref)[m];estimate=(p-v[j3[cases[eid]['ref']['case']],:n])[m]
    response.append(dict(event=eid,xc_km=xc,model=name,**metric(actual,estimate,t[m]),actual_peak_abs=float(abs(actual).max()),predicted_peak_abs=float(abs(estimate).max())))
 # Full trajectories remain local and retain negative reference values.
 packet=dict(time=t,inflow=inf,H=hh,truth=d['flow_m3s'][:,np.array(X)//10],reference=ds['ref']['flow_m3s'][:,np.array(X)//10])
 for name in names:
  packet[name]=pred[name][:,(200-np.array(X))//20]
  packet[name+'_ref']=pred_ref[name][:,(200-np.array(X))//20]
 for name,v in old.items():packet[name]=v[j3[cases[eid]['var']['case']],:n]
 np.savez_compressed(OUT/'data'/f'{eid}.npz',**packet)
 print('Scored',eid,flush=True)

M=pd.DataFrame(metric_rows);D=pd.DataFrame(decomp);R=pd.DataFrame(response);F=pd.DataFrame(phases)
frame('event_metrics.csv',metric_rows);frame('negative_error_decomposition.csv',decomp);frame('paired_response_metrics.csv',response);frame('positive_phase_metrics.csv',phases);frame('initialization_checks.csv',checks)
cols=['NSE','RMSE','KGE','bias','peak_error_pct','peak_time_h','n']
summary=M.groupby(['xc_km','model','window'])[cols].mean().reset_index();summary['events']=M.groupby(['xc_km','model','window']).size().values
summary.to_csv(OUT/'tables/window_summary.csv',index=False,encoding='utf-8-sig')
R.groupby(['xc_km','model'])[['NSE','RMSE','actual_peak_abs','predicted_peak_abs']].mean().reset_index().to_csv(OUT/'tables/paired_response_summary.csv',index=False,encoding='utf-8-sig')
dd=D.groupby('xc_km').agg(total_hours=('total_hours','sum'),negative_hours=('negative_hours','sum'),negative_events=('negative_hours',lambda s:int((s>0).sum())),SSE_all=('SSE_all','sum'),SSE_negative=('SSE_negative','sum'),mean_negative_SSE_fraction=('negative_SSE_fraction','mean'),min_q=('min_q','min'),mean_paired_delta_RMS=('paired_delta_RMS','mean'),mean_paired_delta_peak_abs=('paired_delta_peak_abs','mean')).reset_index()
dd['negative_sample_pct']=100*dd.negative_hours/dd.total_hours;dd['pooled_negative_SSE_pct']=100*dd.SSE_negative/dd.SSE_all
dd.to_csv(OUT/'tables/negative_error_summary.csv',index=False,encoding='utf-8-sig')

rng=np.random.default_rng(20261023);improvements=[]
for (xc,wn),g in M.groupby(['xc_km','window']):
 pivot=g.pivot(index='event',columns='model',values='NSE')
 for nm,base in [('linear_fixed','RWEC0'),('linear_joint','RWEC0'),('legacy1','legacy0'),('legacy1','legacy0_balanced')]:
  if nm not in pivot or base not in pivot:continue
  delta=(pivot[nm]-pivot[base]).dropna().to_numpy();boots=delta[rng.integers(0,len(delta),(10000,len(delta)))].mean(1)
  improvements.append(dict(xc_km=int(xc),window=wn,model=nm,baseline=base,n=len(delta),mean_NSE_gain=float(delta.mean()),median_NSE_gain=float(np.median(delta)),improved_events=int((delta>0).sum()),CI95_low=float(np.quantile(boots,.025)),CI95_high=float(np.quantile(boots,.975))))
frame('paired_model_improvements.csv',improvements)

# Confirm this replay reproduces the previously frozen flood-window scores.
previous=pd.read_csv(P5/'tables/A_event_metrics.csv');match=[]
for name,oldname in [('RWEC0','RWEC-0'),('linear_fixed','RWEC-1_oracle'),('linear_joint','RWEC-1_oracle_joint')]:
 p=previous[(previous.subset=='validation')&(previous.condition=='independent_stage_main')&(previous.model==oldname)]
 q=M[(M.window=='all')&(M.model==name)].merge(p,left_on=['event','xc_km'],right_on=['event','xc_km'],suffixes=('_new','_old'))
 assert len(q)==126
 md=float(abs(q.NSE_new-q.NSE_old).max());assert md<1e-9
 match.append(dict(model=name,rows=len(q),maximum_NSE_difference=md))
scr=pd.read_csv(P3/'tables/section_case_screen.csv');sc=scr[scr.x_km==100]
verification=dict(completed=True,sealed_test_results_opened=0,validation_events=events,new_HEC_runs=0,new_parameter_fits=0,negative_truth_clipped=False,masked_before_state_recursion=False,
 causal_forecast_future_stage_used=False,linear_models_use_future_true_stage=True,old_sources_replayed=source,previous_flood_metrics_reproduced=match,
 x100_all_development_scenarios=len(sc),x100_all_development_negative_samples=int(sc.negative_samples.sum()),all_initialization_checks_passed=all(x['initialization_pass'] for x in checks),protocol_sha256=sha(OUT/'PROTOCOL.md'))
save(OUT/'provenance/verification.json',verification)
save(OUT/'checkpoint.json',dict(status='calculated',validation_events=21,sections=6,tests_opened=0))
print(summary[(summary.window=='nonnegative')&(summary.model.isin(names))].to_string(index=False),flush=True)
