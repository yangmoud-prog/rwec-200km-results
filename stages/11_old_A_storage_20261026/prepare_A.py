from common10 import *
from scipy.integrate import cumulative_trapezoid
import shutil,importlib.util,datetime

def main():
 for name in ['tables','data','predictions','figures','provenance','logs']:(P/name).mkdir(exist_ok=True)
 for name in ['formulas.py','legacy_models.py','test_formulas.py','EQUATIONS.md','EQUATIONS_EDITABLE.tex']:
  shutil.copy2(P10/name,P/name)
 eq=(P/'EQUATIONS.md').read_text(encoding='utf-8').split('\n河道 D 的库区采用')[0].replace('第十阶段方程推导与检查','旧河道 A：潜在滞蓄公式与初始化')
 eq+='\n\n本轮仅使用旧河道 A 的上游流量。下游为独立外部水位强迫，未来真实水位不进入模型。逐 20 km 子段独立率定，全部事件在 −48 h 恒定基流状态初始化，持续递推至完整洪水结束。训练统一采用 t≥0 完整小时区间，另以真实端点均非负的共同掩码作事后诊断，不在倒流处重置状态。\n'
 (P/'EQUATIONS.md').write_text(eq,encoding='utf-8')
 tex=(P/'EQUATIONS_EDITABLE.tex').read_text(encoding='utf-8').split('\\section*{Reservoir closure}')[0]+'\\end{document}\n';tex=tex.replace('Stage 10: Potential storage and moving-domain conservation','Old channel A: independent potential storage comparison');(P/'EQUATIONS_EDITABLE.tex').write_text(tex,encoding='utf-8')
 split=read(P.parent/'04_autonomous_backwater_20261021/provenance/event_split.json');dev=set(split['train']+split['validation']);assert len(dev)==70 and len(split['train'])==49 and len(split['validation'])==21
 old=read(P10/'provenance/frozen_protocol.json');protocol=dict(stage=11,created=datetime.datetime.now().astimezone().isoformat(),river='original A, 200 km, slope 0.05 per mille',train=split['train'],validation=split['validation'],old_sealed_tests_read=0,targets_outlet_km=TARGETS.tolist(),dt_candidates_h=[.25,.5,1.,2.],bounds=old['bounds'],seed=202610260,objective='event-equally weighted full post-origin hourly NSE at each 20 km reach outlet',calls_per_model_reach_dt=2000,total_calls_per_model=80000,independent_parameters_per_20km=True,fitter_reads_validation=False,training_boundary='main independent-stage case only; paired fixed5 evaluated with identical predictions',comparison_windows=['full','nonnegative_hour_endpoints','after_last_negative_plus12h'],positive_mask='Both hourly HEC endpoint flows >=0; same post-prediction mask for every model, no state reset',forecast_inputs=['upstream boundary inflow hydrograph','known pre-event constant-flow initialization','fixed geometry initial volumes','frozen trained parameters'],future_downstream_stage_input=False,initialization_start_h=-48,model_code_source=str(P10),HEC_rerun=False,validation_rule='one locked evaluation after all models frozen; this is further development on previously inspected events, not fresh independent testing')
 coeff=[float(np.log(.0001)),6.];shape=[float(np.log(.15)),float(np.log(150))];delay=[0.,8.]
 protocol['bounds']={'original':[coeff,shape,delay],'new':[coeff,shape,delay],'general':[coeff,coeff,shape,delay],'NLM':[[float(np.log(.05)),float(np.log(20))],[0.,.5],[.4,2.]]}
 protocol['bounds_note']='Set before any stage11 calibration: retain old A k/n/tau search range, extend positive coefficient maximum equally for all RWEC families; no validation-score-based adjustment'
 save(P/'provenance/frozen_protocol.json',protocol)
 plan=read(A/'config/formal_plan.json');cases={(c['event_id'],c['kind']):c for c in plan['formal_cases'] if c['event_id'] in dev and c['role']=='main' and c['kind'] in ['var','ref']}
 spec=importlib.util.spec_from_file_location('oldA_hec_boundary',A/'scripts/hec200.py');hec=importlib.util.module_from_spec(spec);spec.loader.exec_module(hec)
 manifest=[];quality=[]
 for e in split['train']+split['validation']:
  pair={}
  for kind in ['var','ref']:
   c=cases[e,kind];folder=A/'results'/c['case'];file=folder/'reach_timeseries.npz';raw=load(file);t=np.rint(raw['time_since_flood_h']);assert abs(raw['time_since_flood_h']-t).max()<1e-6 and t[0]==-48 and np.all(np.diff(t)==1)
   ql=read(folder/'quality.json');assert ql['basic_quality_pass'] and ql['initialization_pass'] and ql['native_balance_pass'] and ql['tail_recovered']
   x=raw['x_from_outlet_km'];assert np.array_equal(x,np.arange(0,201,10));inp=hec.boundary(c['event'],t,'ref')[0];cum=cumulative_trapezoid(raw['area_m2'][0],x*1000,initial=0);v0=np.array([np.interp(200-j*20,x,cum)-np.interp(180-j*20,x,cum) for j in range(10)])
   pair[kind]=dict(time=t,x=x,q=raw['flow_m3s'],z=raw['stage_m'],area=raw['area_m2'],inflow=inp,initial_volumes=v0)
   manifest.append(dict(event=e,subset='train' if e in split['train'] else 'validation',kind=kind,path=str(file),sha256=sha(file),bytes=file.stat().st_size))
   quality.append(dict(event=e,kind=kind,complete=True,initialization_pass=True,native_balance_pass=True,tail_recovered=True,min_flow_m3s=ql['min_flow_m3s'],start_h=float(t[0]),end_h=float(t[-1]),maximum_native_balance_fraction=max(ql['native_max_balance_residual_fraction_storage_by_length'][1:])))
  assert np.array_equal(pair['var']['time'],pair['ref']['time']) and np.array_equal(pair['var']['inflow'],pair['ref']['inflow'])
  np.savez_compressed(P/'data'/f'{e}.npz',**pair['var'],q_ref=pair['ref']['q'],z_ref=pair['ref']['z'],area_ref=pair['ref']['area'])
 save(P/'SOURCE_MANIFEST.json',dict(sources=manifest,HEC_cases=140,development_events=70,sealed_results_opened=0));frame('HEC_source_quality.csv',quality)
 current=read(W/'CURRENT_EXPERIMENT.json');current.update(analysis_previous_stage10=str(P10),analysis_stage='stage11_original_A_storage_formulas_running',analysis_latest=str(P),analysis_latest_status='in_progress_snapshot',stage11_data=70,stage11_old_sealed_tests_read=0);save(W/'CURRENT_EXPERIMENT.json',current)
 save(P/'checkpoint.json',dict(stage=11,status='data_prepared_training_not_started',train=49,validation=21,HEC_cases_read=140,sealed_test_read=0));print('OLD_A_READY_70_EVENTS_140_MAIN_CASES',flush=True)
if __name__=='__main__':main()
