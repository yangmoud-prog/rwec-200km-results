from common10 import *
import datetime,shutil
LABEL={'original':'原 RWEC','new':'新 SE RWEC','general':'广义 SE RWEC','NLM':'非线性马斯京根'}
def md(df,digits=4):
 def value(v):return f'{v:.{digits}f}' if isinstance(v,(float,np.floating)) and np.isfinite(v) else str(v)
 return '| '+' | '.join(df.columns)+' |\n|'+'|'.join(['---']*len(df.columns))+'|\n'+'\n'.join('| '+' | '.join(value(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
def main():
 cfg=read(P/'provenance/frozen_protocol.json');check=read(P/'provenance/evaluation_complete.json');df=pd.read_csv(P/'tables/validation_dynamic_event_metrics.csv');means=df.groupby(['window','model'])[['NSE','RMSE','KGE']].mean();frame('validation_overall_means.csv',means.reset_index())
 show=means.loc[['full','nonnegative_hour_endpoints']].reset_index();show['window']=show.window.map({'full':'完整洪水','nonnegative_hour_endpoints':'排除真实负值共同小时'});show['model']=show.model.map(LABEL);show.columns=['窗口','模型','平均 NSE','平均 RMSE m³/s','平均 KGE']
 loc=df[df.x_from_outlet_km.isin([100,20,0])];nse=loc.pivot_table(index=['x_from_outlet_km','distance_from_inlet_km','model'],columns='window',values='NSE',aggfunc='mean').reset_index();pos=loc[loc.window=='nonnegative_hour_endpoints'].groupby(['x_from_outlet_km','model'])[['RMSE','KGE','peak_pct','peak_time_h']].mean().reset_index();loc=nse.merge(pos,on=['x_from_outlet_km','model']);loc=loc[['x_from_outlet_km','distance_from_inlet_km','model','full','nonnegative_hour_endpoints','RMSE','KGE','peak_pct','peak_time_h']];loc['model']=loc.model.map(LABEL);loc.columns=['距出口 km','距入口 km','模型','完整 NSE','排负 NSE','排负 RMSE','排负 KGE','排负峰值误差 %','排负峰时误差 h']
 gains=pd.read_csv(P/'tables/event_mean_gains.csv').query("subset=='validation' and scenario=='dynamic'");counts=[]
 for window in ['full','nonnegative_hour_endpoints']:
  for m in ['new','general','NLM']:
   vals=gains[gains.window==window][m+'_minus_original'];counts.append(dict(window=window,model=m,events=len(vals),improved=int((vals>0).sum()),mean=float(vals.mean()),median=float(vals.median()),minimum=float(vals.min()),maximum=float(vals.max())))
 frame('event_gain_summary.csv',counts);count=pd.DataFrame(counts);count['window']=count.window.map({'full':'完整洪水','nonnegative_hour_endpoints':'排除真实负值共同小时'});count['model']=count.model.map(LABEL);count.columns=['窗口','比较模型','场数','改善场数','平均改善','改善中位数','最小改善','最大改善']
 negative=pd.read_csv(P/'tables/negative_flow_windows.csv').query("subset=='validation'");neg=negative.groupby(['x_from_outlet_km','distance_from_inlet_km']).agg(negative_events=('negative_nodes',lambda v:int((v>0).sum())),negative_nodes=('negative_nodes','sum'),excluded_hours=('excluded_hours','sum'),full_hours=('full_hours','sum'),minimum_Q=('minimum_Q_m3s','min')).reset_index();frame('validation_negative_summary.csv',neg)
 # Complete-event positive subgroup is a diagnostic, never a replacement of 21 validation events.
 allpos=negative[negative.whole_event_nonnegative][['event','x_from_outlet_km']];pos_events=df[df.window=='full'].merge(allpos,on=['event','x_from_outlet_km']);frame('whole_event_nonnegative_metrics.csv',pos_events);frame('whole_event_nonnegative_means.csv',pos_events.groupby(['model','x_from_outlet_km','distance_from_inlet_km'])[['NSE','RMSE','KGE']].mean().reset_index())
 for m in MODELS:
  trials=[read(f) for f in (P/'provenance').glob(f'fit_A_{m}_reach*_dt*.json')];assert len(trials)==40 and sum(r['evaluations'] for r in trials)==80000
 newgain=float(means.loc[('nonnegative_hour_endpoints','new'),'NSE']-means.loc[('nonnegative_hour_endpoints','original'),'NSE']);gengain=float(means.loc[('nonnegative_hour_endpoints','general'),'NSE']-means.loc[('nonnegative_hour_endpoints','original'),'NSE']);newcount=int(counts[3]['improved']);params=pd.read_csv(P/'tables/parameters_interpreted.csv');bound=params.groupby('model').boundary_parameter_count.sum().to_dict();n20=neg[neg.x_from_outlet_km==20].iloc[0];n100=neg[neg.x_from_outlet_km==100].iloc[0]
 historical=pd.read_csv(P/'tables/historical_frozen_original_means.csv').query("window=='nonnegative_hour_endpoints'").set_index('x_from_outlet_km');qmin=pd.read_csv(P/'tables/solver_quality.csv').groupby('model').min_Q.min();hist20=float(historical.loc[20,'NSE']);hist100=float(historical.loc[100,'NSE'])
 report=fr'''# 旧河道 A 的潜在滞蓄记忆结构再检验

本轮在最初的 200 km、0.05‰ 旧河道 A 上重新检验潜在滞蓄公式。原式、新式、广义式和非线性马斯京根均重新率定，使用原 49 场训练与 21 场开发内验证，全部预测不读取未来真实下游水位。排除真实负流量后的共同小时中，新式相对原式的七断面平均 NSE 变化为 {newgain:+.5f}，广义式为 {gengain:+.5f}；新式在 21 场中的 {newcount} 场取得七断面平均改善。均值、中位数及全部失败均保留，不能据某一场曲线决定公式优劣。

与前一轮 B、C、D 河道不同，A 的下游水位是独立外部强迫。模型只根据上游洪水和内部状态演进，因此本轮检验的是记忆结构对同一开发分布内的洪水演算是否有效，不能由高 NSE 推断模型已经自主恢复了未来独立水位过程。原 30 场封存测试未启用，原 HEC 工程未重跑。140 条主配对结果均通过既有初始化、原生水量平衡和恢复检查；32 场补充控制没有扩充本轮训练或验证样本。

严格比较的三个潜在滞蓄式为

\[
S_E=S_{{t-1}}+kV_\tau,\qquad S_E=kS_{{t-1}}+V_\tau,\qquad S_E=aS_{{t-1}}+bV_\tau.
\]

它们共享相同 \(S_W\)、互补关系和质量守恒求解。恒定入流下令 \(r=S^*/(I_0\Delta t)\)，各自独立求解

\[
\left(\frac r{{r+1}}\right)^n+\left(\frac r{{ar+b}}\right)^n=1.
\]

原式对应 \((a,b)=(1,k)\)，新式对应 \((k,1)\)。旧公式的数值 k 没有直接移植到新公式；原式的历史训练参数仅作为其搜索的一个起点。推导与局地记忆解释见 [可编辑方程](EQUATIONS.md)，TeX 文档已编译检查。

各模型都允许十个 20 km 子段分别率定参数及 0.25、0.5、1、2 h 时段。每模型、每段、每时段均为 2000 次目标函数调用，即每模型 80000 次，四个模型共 320000 次。为了保留原模型的既有自由度，任何本轮率定前就将正系数范围设为 10⁻⁴—e⁶，n 为 0.15—150、单段延迟 0—8 h；非线性马斯京根的蓄量尺度、入流权重与指数独立率定。下一子段始终使用模型自己的上段预测，不以未来真实内部流量替代。训练目标统一采用完整洪水过程，验证指标没有反过来决定参数范围或追加优化。

模型从 −48 h 恒定基流状态连续递推至各场原始完整终点，评分从 0 h 起。部分下游水位在洪水起点前已经上涨，它造成的真实抑流不能被解释为需要删除的数值启动扰动。按用户之前要求，同时计算两端 HEC 保存流量均非负的共同小时，以及最后一次负值后再缓冲 12 h 的连续窗口；三种窗口均在完整预测形成后评分，未在倒流处剪接或重置状态。后两种是真值条件下的事后诊断，不是可预先知道的业务预报条件，也不等于水力稳态。

全部模型和 HEC 都转换为相同的逐小时区间平均流量。与此前按保存节点评分、共享子段参数的历史数字相比，本轮的观察算子和参数结构有所区别，应以内轮同口径比较为准。完整与排负后平均结果如下，独立单位为 21 场事件，七个断面不扩充独立洪水数。

{md(show)}

![旧河道平均指标](figures/01_A_mean_skill_full_and_nonnegative.png)

距出口 100 km（距入口 100 km）处，21 场验证均没有负流量，最小保存流量为 {n100.minimum_Q:.2f} m³/s，因而完整和排负窗口相同。距出口 20 km（距入口 180 km）处有 {int(n20.negative_events)} 场出现负值，另有 {21-int(n20.negative_events)} 场完整洪水保持非负；这部分完整正流事件另列结果，但没有替换原 21 场验证集合。关键断面的结果为：

{md(loc)}

为与此前记录衔接，另外对第五阶段已经冻结的共享参数原 RWEC 按本轮小时均值和共同掩码重新评分，没有再训练或据此替换本轮模型。该历史参照在 100 km、20 km 处的排负 NSE 分别为 {hist100:.5f}、{hist20:.5f}。它原来在固定／变水位多场景上率定，当前四模型则使用主变水位、逐段局部率定，两者训练目标不同。更高的历史参照也说明，本轮扩大参数数量和逐段搜索并不保证整个级联目标已达到全局最优；不能将不同率定方案的差值归因于公式。历史冻结预测完整单列在 `historical_frozen_original_event_metrics.csv`。

本次四个新率定模型在这些洪水上的预测最小值均为约 {float(qmin.min()):.1f} m³/s，而 20 km 的完整正流事件中，真实流量也会明显降至基流以下。因而误差并不只发生在负流量阶段：在不倒流的区间，独立下游水位造成的抑流仍未由仅有上游输入的这四个模型还原。交换潜在滞蓄式的系数没有提供缺失的外部边界信息。

![代表洪水完整过程](figures/03_A_full_hydrographs_representative.png)

本图在形成验证结果前指定 E006、E037、E094，另图指定 E010、E051、E068；这些事件在 20 km 处属于全程正流背景，出口处仍可能出现负值，真值曲线完整保留。其余每一场验证洪水也各有 `event_场次_full_hydrographs.png`，同时展示 100 km、20 km 和出口。不能只看峰值接近而略去前期抑流和退水额外释放。

下表对每场先平均七个断面的 NSE 改善，再统计全部 21 场，保留最差场次和负改善。

{md(count,5)}

![全部验证洪水的改善](figures/02_A_all21_paired_gains.png)

广义式具有额外自由度，四种模型仍属于有限预算、逐段优化，而非全局最优证明。参数表列明全部边界命中；各模型边界参数累计数为 {bound}。固定参数将各子段时段统一减半或加倍的敏感性另列，不能在最有利的离散时段取得较高精度后就认定参数已具有稳定物理意义。最终 40 组子段参数的恒定流保持检查结果为 {check['steady_all_passed']}，预测失败子段数为 {check['failed_reach_runs']}。所有时段敏感性仍完整保留，没有据其验证精度重新率定。

![固定参数时段敏感性](figures/04_A_fixed_parameter_step_sensitivity.png)

有效滞蓄状态不自动等于几何河道蓄量。本轮同时记录储量异常、储量变化率和数值接口体积；独立时段之间的流量重采样导致的最大累计体积偏差占入流洪水增量体积约 {100*check['cascade_residual_max_fraction']:.4f}%。这项体积单独记账，不充当回水储量，也不用于修改已报告的流量预测。局部守恒及接口补记恒等式均不能单独证明水面闭合成立。

![有效储量与接口诊断](figures/05_A_storage_and_interfaces.png)

对同一洪水，固定与变动边界具有相同上游输入。本轮四个模型没有外部边界状态输入，因而给出同一条配对预测，其预测边界增量恒为零；真实配对出流差则保留在同一固定空间断面比较。固定 5 m 表中的结果明确标记为“相同预测的参照评分”，并非在固定边界上另外优化后的成绩。主体洪水拟合即使较好，也不能由此说明额外边界响应已被预测。

![独立边界响应的能力限制](figures/06_A_boundary_response_identifiability.png)

本轮结论应同时依据完整过程、排负共同窗口、逐场胜负和时段敏感性。若新式只改变平均传播或使个别事件更接近，尚不足以证明其比原式普遍更适合回水河道。由于这些开发内验证事件已用于此前研究，本次只能作为新公式的进一步开发证据，不能冒充封存测试或独立跨场景成功。全部代码、参数、逐场指标与完整预测见 [文件与复现入口](README.md)，前十阶段成果保持原样。
'''
 (P/'REPORT.md').write_text(report,encoding='utf-8');summary=dict(stage=11,river='original A',status='old_A_formula_comparison_completed',train=49,development_validation=21,old_sealed_tests_read=0,HEC_rerun=False,forecast_future_H_input=False,primary_validation_means=means.reset_index().to_dict('records'),new_minus_original_nonnegative_mean=newgain,general_minus_original_nonnegative_mean=gengain,checks=check,validation_not_unseen=True,model_revision_after_validation=False)
 summary['primary_validation_means']=means.reset_index().replace({np.nan:None,np.inf:None,-np.inf:None}).to_dict('records');summary['historical_frozen_reference_nonnegative_NSE']={'x100':hist100,'x20':hist20};save(P/'SUMMARY.json',summary);save(P/'delivery.json',dict(completed=True,stage=11,scientific_result='see report and all per-event metrics',old_sealed_tests_read=0,remote_verification_receipt='sharing/rwec-200km-results/.local/last_upload.json'))
 save(P/'checkpoint.json',dict(stage=11,status='local_comparison_complete',train=49,development_validation=21,models=4,objective_calls_per_model=80000,HEC_rerun=False,future_H_input=False,old_sealed_tests_read=0))
 current=read(W/'CURRENT_EXPERIMENT.json');current.update(analysis_stage='stage11_original_A_storage_formulas_delivered',analysis_latest=str(P),analysis_findings=str(P/'REPORT.md'),analysis_latest_status='local_complete_repository_sync_pending',stage11_summary=str(P/'SUMMARY.json'),forecast_window_h='individual original event full windows, common start 0 h',stage11_old_sealed_tests_read=0);save(W/'CURRENT_EXPERIMENT.json',current)
 print('REPORT_READY',newgain,gengain,flush=True)
if __name__=='__main__':main()
