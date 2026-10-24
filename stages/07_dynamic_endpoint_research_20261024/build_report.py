from common import *
from datetime import datetime

def table(df):
 d=df.copy()
 for c in d:
  if pd.api.types.is_float_dtype(d[c]):d[c]=d[c].map(lambda x:f'{x:.5g}')
 return '| '+' | '.join(d.columns)+' |\n|'+'|'.join(['---']*len(d.columns))+'|\n'+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in d.itertuples(index=False,name=None))

def main():
 flows=pd.read_csv(OUT/'tables/flow_event_metrics.csv');ev=flows[(flows.subset=='validation')&(flows.window=='positive_common')];end=pd.read_csv(OUT/'tables/endpoint_prediction_metrics.csv')
 summary=ev.groupby(['group','model'])[['NSE','RMSE','KGE']].mean().reset_index();frame('flow_overall_validation.csv',summary)
 ends=end[end.subset=='validation'].groupby('group')[['MAE_km','RMSE','coverage','width_km']].mean().reset_index();frame('endpoint_overall_validation.csv',ends)
 parameters=[];optim=[]
 for g in ['A','B']:
  for model in ['RWEC-0','RWEC-L-predicted','RWEC-L-oracle','Muskingum-linear']:
   f=read(OUT/'fits'/f'{g}_{model}.json');parameters.append(dict(group=g,model=model,**f['parameters'],training_loss=f['loss']))
   for r in f['runs']:optim.append(dict(group=g,model=model,seed=r['seed'],loss=r['loss'],de_success=r.get('de_success'),local_success=r.get('local_success',r.get('success')),nfev=r.get('de_nfev',0)+r.get('local_nfev',r.get('nfev',0))))
 frame('parameters.csv',parameters);frame('optimization_runs.csv',optim)
 save(OUT/'fits/parameters_final.json',dict(parameters=parameters,optimization_runs=optim))
 requested={'E006','E010','E014','E026','E032','E037','E042','E051','E066','E067','E068','E081','E094'}
 pair=pd.read_csv(OUT/'tables/indistinguishable_pair_bounds.csv');same=set(pair.event);strict=set(pair[pair.past_profile_difference_max_m==0].event)
 subgroup=[]
 for (g,model,xc),dd in ev.groupby(['group','model','xc_km']):
  groups={'all_validation_positive_hours':dd,'whole_forecast_positive_events':dd[dd.whole_forecast_positive]}
  if g=='A':groups.update(requested13=dd[dd.event.isin(requested)],same_past_profile_pairs=dd[dd.event.isin(strict)])
  for name,sub in groups.items():subgroup.append(dict(group=g,model=model,xc_km=int(xc),subgroup=name,events=len(sub),NSE=sub.NSE.mean(),RMSE=sub.RMSE.mean(),KGE=sub.KGE.mean()))
 frame('validation_subgroups.csv',subgroup)
 P=pd.read_csv(OUT/'tables/mechanism_process.csv')
 for g in ['A','B']:
  sub=P[(P.group==g)&(P.xc_km==20)&(P.subset=='validation')][['event','time_h','lb_km','delta_volume_m3','deltaQ_m3s']]
  frame(f'{g}_x20_validation_mechanism.csv',sub.to_dict('records'))
 mechanism=pd.read_csv(OUT/'tables/mechanism_events.csv');moving=pd.read_csv(OUT/'tables/moving_boundary_balance.csv');proxy=pd.read_csv(OUT/'tables/volume_proxy_metrics.csv')
 endpoint=pd.read_csv(OUT/'tables/endpoint_events.csv');ep=endpoint[endpoint.threshold_m==.1];sens=pd.read_csv(OUT/'tables/fixed_parameter_timestep.csv')
 lag=mechanism.copy();lag['flood_peak_lag_h']=lag.peak_time_h-lag.inflow_peak_time_h;lag['volume_peak_minus_length_peak_h']=lag.volume_peak_time_h-lag.length_peak_time_h
 frame('peak_timing_relations.csv',lag[['group','event','subset','xc_km','flood_peak_lag_h','volume_peak_minus_length_peak_h']].to_dict('records'))
 boundsummary=pair[pair.event.isin(strict)].groupby('xc_km')[['common_prediction_RMSE_lower_bound_m3s','max_flow_gap_m3s']].mean().reset_index()
 transfer=pd.read_csv(OUT/'tables/zero_shot_A_to_B_flow.csv').groupby('model')[['NSE','RMSE']].mean().reset_index();transferend=pd.read_csv(OUT/'tables/zero_shot_A_to_B_endpoint.csv')
 mm=moving.groupby('group')[['moving_mass_residual_RMSE_m3s','omitted_sweep_residual_RMSE_m3s']].mean().reset_index()
 p20=proxy[(proxy.subset=='validation')&(proxy.xc_km==20)].groupby(['group','diagnostic'])[['NSE','RMSE']].mean().reset_index()
 stat=ep.groupby('group')[['maximum_length_km','mean_length_km','censored_hours','zero_outlet_but_nonzero_length_hours']].mean().reset_index()
 water=pd.read_csv(OUT/'tables/routing_quality.csv');physical=water[water.model!='Muskingum-linear'];maxmass=physical.water_balance_max_m3.max()
 ss=sens.groupby(['group','model','dt_h'])[['NSE','RMSE_vs_hourly']].mean().reset_index()
 base20=ev[(ev.xc_km==20)&(ev.model=='RWEC-0')];sub20=pd.DataFrame(subgroup);sub20=sub20[(sub20.xc_km==20)&(sub20.group=='A')]
 bseed=sum(not bool(r['de_success']) for r in optim if r['de_success'] is not None)
 report=fr'''# 固定起点条件下的动态回水末端与洪水演进：首轮研究结果

2026年10月24日。本轮已经完成末端识别、蓄泄规律检验、固定起点末端范围预测、动态长度RWEC候选演算及开发内验证。结果尚不支持把“长度只调节传播延迟”作为优于原RWEC的新方法。较明确的物理证据是：配对附加蓄量的变化率决定配对流量差；随末端移动的控制体必须包含扫掠体积。末端可以预测到一定程度，并不意味着将该长度乘入传播时间就一定改善出流。

本轮使用旧河槽A的49场训练和21场开发内验证、新河槽B的14场训练和6场开发内验证。A坡度0.05‰，下游3—6 m涨幅独立施加；B沿用已冻结的0.075‰工程及固定控制关系 $H_0=3.5+1.5(Q_0/1000)^{{2/3}}$。两者分别建库、率定和检验，没有混合率定。未新增HEC-RAS计算，30场封存测试结果读取数为0，未制作动画。

预报起点固定为上游洪水起点0 h，出口水位只使用此前48 h至起点的历史，发出预报后不更新。允许完整上游洪水预报，统一比较0—240 h。这样不会把由未来水位决定的模拟结束时间带入特征。训练场次预测末端时按事件留一，验证场次只查询训练库。未来真实水位、内部水面线、流量、蓄量和实际末端均未进入可运行候选的预测输入。

## 末端识别与蓄泄联系

主末端仍以动态工况相对固定5 m参照的水面差大于0.10 m的最远位置定义，同时计算0.05与0.20 m、出口连通长度、影响区数量和200 km删失。逐场峰值及过程均值的等事件汇总如下；“删失小时”是每场0—240 h内的均值。

{table(stat)}

A的影响范围较长，部分时刻到达计算域上边界。A平均每场约{stat.loc[stat.group=='A','zero_outlet_but_nonzero_length_hours'].iloc[0]:.1f}个输出时刻，出口水面差已经低于阈值，上游仍有残余影响。因此“最远影响位置”与“从出口连续连接的回水段”在退水末期不能混用。没有出现同一时刻多个分离超阈值区，但单个残留区可以与出口脱离；两者并不矛盾。阈值末端突然消失会产生长度跳变，这也是将其直接视作平滑传播波速的限制。

为检验长度与传播时间的联系，另在各自训练集拟合“峰现时差—演算河长”的线性关系，再加入真实平均有效回水长度。A验证峰现时差RMSE由5.95 h变为5.98 h，B由1.83 h变为1.78 h；该简单关系没有显示稳定、明显的跨事件增益。峰现时差还受到洪峰变形及多峰响应影响，不能直接等同于特征线传播时间。相关系数、参数和逐场峰时均保留为机理证据，尚未用作可运行预报输入。

在固定预测断面上，两工况入流相同，连续方程给出 $\Delta Q=-d\Delta V/dt$。按六断面逐场检验，A的对应相关系数接近0.999998，20 km处差分残差RMSE均值约0.80 m³/s；B约0.055 m³/s。蓄入使出流受抑，附加蓄量回落时出流恢复，不能用长度峰值或洪峰时刻替代整个蓄泄过程。逐场峰时差见 `tables/peak_timing_relations.csv`。

只用有效长度表示体积，与加入真实高度乘长度的机理代理，在20 km处的未见验证结果为：

{table(p20)}

这些是使用真实几何与真实水位的诊断结果，不能作为自主预报精度。A中高度与长度共同解释体积通常优于仅长度；B不保证同样关系，水面形状与较小的附加蓄量仍有作用。几何附加蓄量也尚未被证明等于RWEC的潜在滞蓄。

移动控制体满足 $dS_M/dt=I-Q(b,t)-A(b,t)\dot b$。有、无扫掠项的逐场残差RMSE等事件均值为：

{table(mm)}

包含扫掠项后，残差为逐时插值与积分精度量级；省略后明显失衡。这里的相对界面通量与固定断面物理流量不同。退缩或阈值消失时相对通量可为负，不能据此认定河水发生真实倒流，更不能将图中的大扫掠尖峰直接加入固定断面的出流预报。

## 固定起点能预测多准确

首轮采用8邻居训练库类比，固定采样和尺度，不用复杂机器学习模型。下表的MAE和RMSE单位均为km，覆盖率为逐时真实长度落在邻居10%—90%经验范围中的比例；宽度单位km。各事件等权。

{table(ends)}

A的未知外部边界造成较宽的未来可能范围。B的控制关系固定，末端响应更可预测，但验证仅6场；其高覆盖率也包含范围偏宽及大量近零末端时段的贡献。经验范围没有经过概率校准，不能称为具有80%保证的置信带。

旧河槽中找到{len(strict)}组历史水面完全相同、同一上游洪水却施加不同未来边界的配对；另有{len(same)-len(strict)}组仅在0.001 m容差内相同。对完全相同配对，任意共同确定性预报都不能同时消除两个响应的差异。以下给出各配对最优共同预报RMSE下界的等事件均值；这是双场景平均意义的不可辨识界，不是单场NSE上限，也不是21场验证的误差下限。

{table(boundsummary)}

这一结果说明，“给定同一上游洪水是否存在唯一未来出流”必须由下游边界性质决定。增加回水长度变量本身不能创造未知外部强迫的信息。

## 动态长度接入RWEC后的实际结果

候选保留水能互补方程，将200 km划为10个20 km子段，以各段有效回水占比调节传播滞时 $\tau_j=\tau_{{20}}(1+\beta r_j)$。可变延迟用单调出发时钟和在途水量处理，不重复释放入流水量。此候选检验的是长度对传播时间的独立贡献；它尚未闭合整个移动双区的蓄泄模型。另列真实未来长度oracle，只判断现有结构在理想长度信息下的表现。

所有模型在各自训练集上联合率定全部参数，原RWEC的 $k,n,\tau_{{20}}$ 与新模型的 $k,n,\tau_{{20}},\beta$ 均一起调整。每场、每断面等权，主评分取真实非负时刻。完整洪水过程仍连续递推，不删改负值、不在转正时重新初始化，也没有用验证结果挑参数。采用三个固定种子的差分进化和局部细化；{bseed}次全局搜索达到预算或未达到全局停止标准，局部结果与各次训练损失已保存，不能据此证明全局最优。

各场六断面等权平均的开发内验证结果为：

{table(summary)}

旧河槽预测长度候选没有稳定增益，平均NSE从0.90486变为0.90423；真实长度oracle也只达到0.90667。说明本轮限制既包含末端预测误差，也包含“仅调节延迟”的表达不足。新河槽预测长度与真实长度oracle均将新增系数率定为0，恢复原RWEC，未证明长度因子的独立价值。线性马斯京根只是一个传统基线，不是此前的RWEC-1线性水位修正；未来水位修正的高精度不能作为本轮可运行模型的成绩。

20与100 km的主要结果为：

{table(pd.read_csv(OUT/'tables/flow_validation_summary.csv').query('xc_km in [20,100]')[['group','model','xc_km','NSE','RMSE','KGE']])}

本轮采用六断面共用参数、固定240 h窗口，目标与此前单独100 km率定不同，数值不能直接当作原先结果被推翻。主比较已经排除真实负流量时刻，近出口仍有明显误差；正流条件下的前期抑流、洪峰衰减与后期恢复仍需解释。既有13场正流子组和起点前无外部扰动子组均另存，不参与重新挑选或率定。

当前候选虽然使用历史水位来选择长度轨迹，却未将历史水位完整转换为起点真实几何储量与内部流量状态。起点前已经受边界扰动的事件可能有明显状态偏差；不能把这部分误差全部归因于未来不可辨识性。后续需要在允许的历史信息下作状态初始化，并与未来边界未知造成的误差分开检验。

## 迁移与数值限制

为避免A/B混合率定，另外冻结A的49场长度库及A的RWEC参数，直接预测B的6场验证洪水，未使用B标签修正。末端MAE约{transferend.MAE_km.mean():.2f} km，流量结果为：

{table(transfer)}

该未经几何和边界适配的直接迁移未通过。其结果不能证明任何物理长度模型都不能迁移，但足以说明当前A类比库不能直接用于B。分别在B训练得到较好成绩，不等同于从A向B成功迁移。

动态候选的离散水量账最大累计残差为{maxmass:.3g} m³，所有场次未触发在途量限制，未出现非收敛迭代，因为这个原型是显式递推而非水位—出流联合迭代。守恒通过不等于精度提高。相同参数从1 h改为0.5、0.25 h，RWEC结果明显变化，线性马斯京根相对稳定：

{table(ss)}

原互补式的储量尺度随时间步长、初始化和步平均约定变化，本轮不能宣称已获得时间步长独立的连续模型。需先建立可解释的时间尺度或明确参数随步长的转换关系，再讨论推广；不能为每个时间步重新率定后称为敏感性通过。

## 本轮结论与继续研究的位置

末端识别、蓄泄规律、固定起点范围预测和第一种动态长度候选的实际演算均已完成。完整动态河长RWEC的物理闭合仍未完成，本轮结果没有证明新公式优于原模型，也没有通过跨边界直接迁移及固定参数时间步长检查。

继续研究应从已验证的 $\Delta Q=-d\Delta V/dt$ 与移动控制体守恒出发，识别长度、剖面形状和初始储量分别提供的信息，再构造只由允许输入递推的蓄泄状态。需要同时保留原RWEC、仅改延迟、理想真实长度这三个对照，才能区分“长度预测不准”和“长度进入方程的方式不对”。对于旧河槽独立外部水位，范围预报和不可辨识界继续保留；固定控制的新河槽才适合检验内部水位—长度—出流闭合。当前不扩大HEC批量、不启用30场封存测试，先以现有数据完成这一闭合检验。

完整方程及单位见 [EQUATIONS.md](EQUATIONS.md)，率定参数与优化记录见 `tables/parameters.csv`、`tables/optimization_runs.csv`。全部逐场指标在 `tables/flow_event_metrics.csv`，末端指标在 `tables/endpoint_prediction_metrics.csv`；完整时序保存在本地 `data/`、`predictions/`。各脚本按 `prepare.py`、`predict_endpoint.py`、`calibrate.py --group A/B`、`evaluate.py`、`transfer.py`、`plot_results.py`、`build_report.py` 顺序可复核；两组率定可并行运行。

图1为E006水面差和蓄泄关系；图2为固定起点末端预测及全部开发内验证事件；图3为模型均值±事件标准差和逐场配对增益，误差条不代表置信区间；图4为预先选定的A-E037及B-E006过程，所有真实值均保留；图5为E006的24/72/120 h沿程水面线及移动末端相对通量，虚线水面为固定5 m参照。图件为Python生成的实际计算结果，英文数字Times New Roman、中文宋体，PDF/SVG文字可编辑，PNG为预览，TIFF为高分辨率输出。
'''
 (OUT/'REPORT.md').write_text(report,encoding='utf-8')
 save(OUT/'summary.json',dict(date='2026-10-24',origin_h=0,history_h=[-48,0],forecast_h=[0,240],A_train=49,A_validation=21,B_train=14,B_validation=6,test_opened=0,endpoint_validation=ends.to_dict('records'),flow_validation=summary.to_dict('records'),same_past_strict_pairs=len(strict),same_past_tolerance_pairs=len(same),mass_balance_max_m3=float(maxmass),dynamic_delay_stable_gain=False,zero_shot_transfer_pass=False,fixed_parameter_timestep_pass=False,full_moving_boundary_closure_complete=False,animation=False))
 print('Strict same-past pairs',len(strict),'tolerance',len(same));print('Report and summary written')

if __name__=='__main__':main()
