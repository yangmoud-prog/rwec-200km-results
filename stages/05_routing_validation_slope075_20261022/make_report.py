from common import *
import math
def table(df,digits=4):
    def fmt(v):return f'{v:.{digits}f}' if isinstance(v,(float,np.floating)) else str(v)
    return '| '+' | '.join(df.columns)+' |\n|'+'|'.join(['---']*len(df.columns))+'|\n'+'\n'.join('| '+' | '.join(fmt(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
def main():
    a=pd.read_csv(OUT/'tables/A_event_metrics.csv');b=pd.read_csv(OUT/'tables/B_event_metrics.csv');qa=pd.read_csv(OUT/'tables/B_HEC_quality.csv');sens=pd.read_csv(OUT/'tables/B_HEC_numerical_sensitivity.csv');response=pd.read_csv(OUT/'tables/B_20pilot_response.csv')
    av=a[(a.subset=='validation')&(a.condition=='independent_stage_main')];bv=b[b.subset=='validation'];models=['RWEC-0','RWEC-L','RWEC-V','RWEC-H','RWEC-1_oracle','RWEC-1_oracle_joint']
    summary=[];section=[];pars=[];numeric=[]
    for label,d,primary in [('A',av,a),('B',bv,b)]:
        for model in models:
            m=d[d.model==model]
            summary.append(dict(河槽=label,模型=model,独立洪水数=m.event.nunique(),NSE=float(m.NSE.mean()),RMSE_m3s=float(m.RMSE_m3s.mean()),KGE=float(m.KGE.mean()),洪峰误差_pct=float(m.peak_error_pct.mean()),峰时误差_h=float(m.peak_time_error_h.mean())))
        g=d[d.model=='RWEC-0'].groupby('xc_km')[['NSE','RMSE_m3s','KGE']].mean().reset_index();g.insert(0,'河槽',label);section.extend(g.to_dict('records'))
        p=pd.read_csv(OUT/'tables'/f'{label}_parameters.csv');p=p[~p.model.str.contains('_joint')];p.insert(0,'河槽',label);pars.extend(p.to_dict('records'))
        n=pd.read_csv(OUT/'tables'/f'{label}_numerical_checks.csv');numeric.append(dict(experiment=label,nonconverged_steps_max=int(n.nonconverged_steps.max()),iterations_max=int(n.max_iterations.max()),mass_error_m3=float(n.max_mass_residual_m3.max()),delay_balance_error_m3=float(n.max_full_balance_residual_m3.max()),domain_error_steps=int(n.out_of_lookup_steps.max())))
    summary=pd.DataFrame(summary);section=pd.DataFrame(section);parameters=pd.DataFrame(pars)
    frame('validation_primary_comparison.csv',summary.to_dict('records'));frame('validation_by_section.csv',section.to_dict('records'));frame('primary_parameters.csv',parameters.to_dict('records'))
    lt=pd.read_csv(OUT/'tables/B_length_metrics.csv');lv=lt[(lt.subset=='validation')&(lt.model=='RWEC-L')]
    tr=pd.read_csv(OUT/'tables/A_to_B_zero_refit_transfer.csv');tr=tr[(tr.subset=='validation')&(tr.model=='RWEC-0')]
    bounds=pd.read_csv(OUT/'tables/A_information_bounds_flood_window.csv');bg=bounds[(bounds.scenarios==10)&(bounds.xc_km==20)]
    six=[]
    for label in ['A','B']:
        d=pd.read_csv(OUT/'tables'/f'{label}_six_target_metrics.csv');d=d[d.subset=='validation']
        if label=='A':d=d[d.condition=='independent_stage_main']
        for model,g in d.groupby('model'):six.append(dict(河槽=label,模型=model,NSE=float(g.NSE.mean()),RMSE_m3s=float(g.RMSE_m3s.mean())))
    six=pd.DataFrame(six);tstep=[]
    for label,d in [('A',av),('B',bv)]:
        baseline=d[(d.model=='RWEC-0')&(d.xc_km==20)];tstep.append(dict(河槽=label,dt_h=1.,NSE=float(baseline.NSE.mean()),RMSE_m3s=float(baseline.RMSE_m3s.mean()),与1h过程差_RMSE_m3s=0.))
        s=pd.read_csv(OUT/'tables'/f'{label}_timestep_sensitivity.csv');s=s[(s.model=='RWEC-0')&(s.xc_km==20)]
        if label=='A':s=s[s.condition=='independent_stage_main']
        for dt,g in s.groupby('dt_h'):tstep.append(dict(河槽=label,dt_h=float(dt),NSE=float(g.NSE.mean()),RMSE_m3s=float(g.RMSE_m3s.mean()),与1h过程差_RMSE_m3s=float(g.difference_from_1h_RMSE_m3s.mean())))
    frame('timestep_comparison.csv',tstep)
    slope=read(B/'config/pilot_frozen.json');qc=read(B/'logs/pilot_checkpoint.json');nc=read(B/'logs/numerical_checks.json');Alook=loadnp(OUT/'data/A_lookup.npz');Blook=loadnp(OUT/'data/B_lookup.npz')
    equations=(P4/'EQUATIONS.md').read_text(encoding='utf-8').replace('实验B从出口传到200 km的初始近似时间为约12.50 h',f"本轮0.075‰实验B从出口传到200 km的初始近似时间为约{Blook['delay_h'][-1]:.2f} h")
    equations=equations.replace('四场的数值成功','已计算样本的数值成功')
    equations+='\n\n## 本轮验证口径\n\n旧A坡度0.00005，新B坡度0.000075；两套水面线查表、传播时间和率定相互独立。B为20场洪水、14训练／6开发内验证；A仍为49／21。主目标在6个指定内部断面及出口等权计算，另做只含6个内部断面的训练权重敏感性，后者为训练诊断触发的探索性分支，不替换冻结主结果。\n\nA未知独立外部边界下，5 m持续假设令附加项为零。移植C(H)只能作为明确的控制假设失配压力试验，不能宣称是A已知物理边界。参考真值和oracle只进入评价或标明的理想信息对照。主增益边界0≤α≤1；逐位置网格和固定非零α过程均保留，不用非物理负增益追求高NSE。\n'
    (OUT/'EQUATIONS.md').write_text(equations,encoding='utf-8')
    report=fr'''# 旧河道实际洪水演算与0.075‰河槽二十组HEC先导结果

2026年10月22日。已完成实际求解、洪水演算、逐场评价及静态图件；未启用30场封存测试。

本轮完成了旧200 km河道70场开发洪水的172条配对和补充控制场景评价，以及新0.075‰河槽20组配对先导、40次主求解和4次时间／空间检查。RWEC-0、L、V、整段长度H消融及真实未来水位oracle均实际执行，A另做5 m持续边界和固定非零增益压力试验。计算没有因既有时间步长问题中止。主要结果是：固定控制下可以自主预测动态回水范围，但在让原RWEC获得同等训练输出信息后，本轮简单正向潜在项修正没有保留非零共享增益。高流量NSE不能用来证明长度修正的独立价值。

## 实验对象与验证信息

新河槽为200 km、坡度0.075‰、底宽500 m、边坡3、糙率0.025、2 km断面、60 s计算、1 h保存。初始流量1000 m³/s、出口初始水位5 m，冻结的单值控制为

$$H_0(Q_0)=3.5+1.5(Q_0/1000)^{{2/3}}\quad[\mathrm m].$$

它是已知的理想堰流型评级关系，并未建立水库库容、调度或实体淹没堰模型。未来水位不作为时间序列输入HEC或自主RWEC；HEC的静态评级表在求解前已固定，自主模型由自己的出口流量求水位。涨幅1—4 m是设计目标，没有根据本20场未来结果反调参数。20组每组包含固定控制和固定5 m参照，共40次主计算。前四背景保留，余下依入流洪峰、历时、峰现比例的归一化最远点规则选取；14场来自原49训练集，6场来自原21开发内验证集，未访问原30场封存结果。

20组均已完成并通过本轮初始化、原生流量体积积分水量平衡、边界误差、断面覆盖和退水恢复检查。新河槽实际出口涨幅{response.outlet_stage_rise_m.min():.3f}—{response.outlet_stage_rise_m.max():.3f} m，参考回水最大长度{response.peak_reference_lb_km.min():.2f}—{response.peak_reference_lb_km.max():.2f} km，全河段最小流量{qa.min_flow_m3s.min():.3f} m³/s，无裁剪流量操作。20 km处最大抑流约{abs(response.suppression20_m3s.max()):.2f}—{abs(response.suppression20_m3s.min()):.2f} m³/s，退水恢复增流约{response.recovery20_m3s.min():.2f}—{response.recovery20_m3s.max():.2f} m³/s。

![20场水动力响应](figures/B_twenty_pilot_hydraulics.png)

E023另计算30 s时间步长和1 km空间间距，各含控制／参照两条，共4次附加HEC计算。相对于较细设置，原2 km／60 s在六个内部位置的最大RMSE为{sens.base_vs_finer_RMSE_m3s.max():.3f} m³/s，峰值相对差的最大绝对值{abs(sens.peak_difference_pct).max():.4f}%，最大水面差{sens.max_stage_difference_m.max():.6f} m；四次均通过基本质量检查。这是代表事件的数值敏感性证据，不能当作所有未采样水力状态均已验证。40场主计算原生平衡残差相对储量变化最大约{100*qa.native_balance_max_fraction.max():.5f}%。完整质量、分组及四个细化过程见[B_HEC_quality.csv](tables/B_HEC_quality.csv)、[B_HEC_numerical_sensitivity.csv](tables/B_HEC_numerical_sensitivity.csv)。

旧A没有重跑HEC。其已有固定5 m、独立3—6 m涨幅过程和32次控制输出全部保留为评价真值，仍用49训练／21开发内验证。A不存在唯一已知的出口控制，因此实际执行两种诚实的无未来水位闭合：5 m持续假设，以及将冻结堰流关系移植到A的压力试验。后者产生模型内部长度和体积，但控制假设与真实外部强迫不一致，不能冒称A真实固定边界。同一规定上游洪水在不同真实边界下始终使用同一条因果预测。负流量真值保留，详见[truth_reversal_by_case.csv](tables/truth_reversal_by_case.csv)。

![同一上游输入面对不同外部边界](figures/A_same_inflow_unknown_boundary.png)

## 公式、率定和实际演算

水能互补递推保持本工程原实现：

$$S_W=S_{{i-1}}+\overline I_{{\tau,i}}\Delta t,\qquad S_E=S_{{i-1}}+k\overline I_{{\tau,i}}\Delta t+B_i,$$
$$S_i=(S_W^{{-n}}+S_E^{{-n}})^{{-1/n}},\qquad \widehat Q_i=(S_W-S_i)/\Delta t.$$

L与V候选分别使用

$$B_{{L,i}}=\alpha B_0[\widehat H_0-5]_+L_{{\mathrm{{eff}},i}},\qquad B_{{V,i}}=\alpha\int\{{A(\widehat y)-A(y_r)\}}\,dx.$$

它们不叠加。几何面积差只是潜在项代理，不被重复加到实际储量，也不被假定等于RWEC状态。全河十个20 km子段由200 km向0 km递推，0 km预测流量闭合控制关系，再由渐变流和初始上行特征传播时间预测水面差；0.10 m阈值确定动态长度。预测函数没有真实未来水位、HEC内部出流或实际末端参数。完整单位、初始化、水量记账、迭代和信息边界推导见[EQUATIONS.md](EQUATIONS.md)。

为严格比较，本轮先让RWEC-0也在各自的训练场景输出上率定，避免只让新公式额外接触动态场景。A每个上游事件等权、同事件内两种或十种边界等权；B只用14个固定控制动态训练事件。主目标含六个指定内部位置及出口，各位置等权。随后冻结公共核心，分别率定每个正向潜在项的单一增益。还对每个模型执行同为70次目标评价预算的局部联合更新；它是同预算敏感性，不是全局最优证明。主核心先采用固定预算差分进化及局部细化，边界为$10^{{-4}}\le k\le20$、$0.15\le n\le150$、$0\le\tau_{{20km}}\le8$ h，$0\le\alpha\le1$。

A/B完全分开率定。主参数如下，所有内部位置共享相同子段参数，未逐位置另配一套最优参数。

{table(parameters[['河槽','model','k','n','tau_per20km_h','alpha']],6)}

## 原RWEC与新公式的验证结果

下表A为21场主动态边界事件、B为6场固定控制事件，在六个内部位置的事件均值。位置之间相关，不能把A的126条或B的36条过程当作126／36场独立洪水。RWEC-1使用真实未来水位，是oracle回算；它的有限公式结构不构成所有模型性能的数学上界。

{table(summary,5)}

本轮冻结主目标下，L、V和H的最优共享增益均取0；局部联合更新没有找到更好的非零自主解，因此它们的主流量过程与RWEC-0相同。这个结果不是没有执行新公式：A额外完成了$\alpha=0.01,0.05$的所有场景实际演算，两河槽均完成0至1的非零增益网格；详细过程误差、稳定性和训练损失见[A_nonzero_gain_stress_metrics.csv](tables/A_nonzero_gain_stress_metrics.csv)、[A_alpha_grid_by_event.csv](tables/A_alpha_grid_by_event.csv)、[B_alpha_grid_by_event.csv](tables/B_alpha_grid_by_event.csv)。单个事件或位置仍可局部改善，不能据此替换整个共享目标的选择。

![A实际演算与非零增益](figures/A_actual_routing_and_nonzero_gains.png)

对不同演算长度，RWEC-0及主L/V的验证均值如下。A越靠近出口误差越大，不能用上游位置较高的NSE掩盖出口附近的未知强迫；B则在已知控制下有较高小时精度。

{table(section,5)}

![B六场验证洪水过程](figures/B_six_validation_flows.png)

训练网格提示不同位置对增益的偏好不同，V在某些远端位置取较大的网格端点但总体改善很小。因此另做只含六个指定内部断面、排除出口误差权重的率定敏感性，重新率定原RWEC并分别搜索L/V，仍只用训练事件。该分支由训练诊断触发，属于探索性权重检查，不替换冻结主结果，也没有按验证精度挑选参数。其结果如下；L/V仍未保留非零共享增益，说明主结论不只是额外出口权重造成的。

{table(six,5)}

位置增益网格与边界命中见[A_position_alpha_stability.csv](tables/A_position_alpha_stability.csv)、[B_position_alpha_stability.csv](tables/B_position_alpha_stability.csv)。仅强制各位置共享参数不能证明参数稳定；本轮边界命中及零增益正是需要保留的可辨识性结果。

## 长度、水面线与迁移检验

B的六场验证，内部计算长度RMSE为{lv.RMSE_km.min():.3f}—{lv.RMSE_km.max():.3f} km，均值{lv.RMSE_km.mean():.3f} km；出口水位RMSE约{lv.H_RMSE_m.min():.4f}—{lv.H_RMSE_m.max():.4f} m。L/V增益虽然为零，内部水力状态和长度仍由各自预测出口过程计算，未使用真实末端，因而可以独立评价水力模块。它的长度精度不能自动转化为附加潜在项的流量增益。

![六场回水长度验证](figures/B_six_validation_lengths.png)

![沿程水面线](figures/B_longitudinal_profiles.png)

绝对水面线仍使用单一出口代表流量的准恒定基准，与非均匀非恒定HEC剖面存在偏差。几何体积代理逐位置的误差单独列于[A_profile_volume_errors.csv](tables/A_profile_volume_errors.csv)、[B_profile_volume_errors.csv](tables/B_profile_volume_errors.csv)，不将它们称为RWEC潜在储量真值。A阈值长度还有200 km右删失和可能不连通影响区，不能把最远超阈值点一概解释为连续传播波前。

另实际进行了A参数到新B的零重率定迁移，不混合两套目标。B六场验证中，迁移原参数的平均NSE={tr.NSE.mean():.5f}、RMSE={tr.RMSE_m3s.mean():.2f} m³/s，明显弱于在B独立训练后的结果。由于A学得的L/V增益为零，该检验不能证明非零长度作用可迁移，只说明共享参数未经重率定的表现。完整逐场和长度结果见[A_to_B_zero_refit_transfer.csv](tables/A_to_B_zero_refit_transfer.csv)、[A_to_B_zero_refit_length.csv](tables/A_to_B_zero_refit_length.csv)。

## 信息边界、守恒与失败范围

四个十边界背景的同一入流不可辨识性重新按洪水起点后共同保存时窗计算。20 km的等权经验RMSE下界为{bg.conditional_RMSE_floor_m3s.min():.2f}—{bg.conditional_RMSE_floor_m3s.max():.2f} m³/s。它们以洪水前48 h的相同物理初态为条件；误差从$t=0$以后统计，并不意味着提供了不同场景在$t=0$的完整内部状态。有限场景等权也不是现实概率。任何自主确定性模型都不能同时命中同一输入下的所有不同真实响应，详见[A_information_bounds_flood_window.csv](tables/A_information_bounds_flood_window.csv)。

两河槽主模型均无记录的不收敛或查表带外步，最大迭代23次。子段离散平衡和延迟在途体积平衡达到浮点误差量级；独立累计的在途体积平衡最大残差分别约{numeric[0]['delay_balance_error_m3']:.2e}、{numeric[1]['delay_balance_error_m3']:.2e} m³。它们验证算法记账，不证明算法状态等于HEC几何储量。级联接口的时段平均／梯形求积差也须按方程说明处理。改变24 h以后的上游输入，全部自主模型此前的结果不变；真实水位只进入标明的oracle。数值失败表与输入核验见[provenance/verification.json](provenance/verification.json)、[provenance/causality_check.json](provenance/causality_check.json)。

固定参数时间步长检验依然失败。下表取20 km各验证洪水均值；主L/V退化为原模型，因此具有相同变化。没有用各步长重新率定掩盖这种结构差异。

{table(pd.DataFrame(tstep),5)}

![B消融与固定参数步长](figures/B_ablation_and_timestep.png)

已知固定单值控制、初态稳态和温和亚临界正流使B的水力闭合可用；独立外部水位、未知调度、淹没控制、明显倒流、近临界或超临界状态、快速变化和沿程非均匀流量仍可能使当前模型失效。A已包含实际倒流反例，B20场无倒流不能代表那些状态已经通过验证。时间离散不一致是已经证实的另一项失效，而不是HEC求解失败。

本轮因此回答了四个研究问题：无未来水位时，A存在不可消除的信息误差，B在本固定控制样本的1 h演算可以取得较高精度；简单长度／体积附加项在同等训练信息下未显示稳定独立增益；零重率定迁移已计算但不支持非零长度因子跨河槽成功的结论；水力闭合、未知外部强迫和时间尺度均构成明确适用限制。它们是实际运行得到的结果，不把高NSE、参数取零或计算返回成功替代方法成立的证据。

## 完整结果入口

主指标表：[A逐模型文件索引](tables/A_metric_file_index.csv)、[B逐模型文件索引](tables/B_metric_file_index.csv)。本地A合表含12384行，阅读仓库按模型拆分以保留全部逐场内容；不存在因大文件限制丢掉低精度事件。另有洪峰与峰时、抑流恢复、长度、水位、几何体积、参数、数值及步长表。两河槽从未混合率定，原70／30及49／21分组保留。

[全部逐场静态图册](FIGURE_ATLAS.md)包含A全部21场主动态验证洪水的六位置过程及B全部20场的流量／长度过程。英文、数字和数学符号Times New Roman，中文SimSun，PDF/SVG保留文字，方程保留可编辑LaTeX；本轮无动画。原始HEC/HDF和完整预测过程保留在移动硬盘，阅读包保存正文、关键表格和PNG。工程control工况旧boundary.csv中的stage_m=5只是生成器初始占位，实际施加边界以冻结u01静态Rating Curve为准。
'''
    (OUT/'REPORT.md').write_text(report,encoding='utf-8')
    atlas=['# 逐场静态过程图册\n\n全部事件均保留，未按精度筛选。A中灰线代表RWEC-0与退化为零增益的主L/V；蓝／橙虚线是固定非零增益压力试验。B包含全20场，上游输入均从200 km开始，图册显示20 km流量和全程回水长度，其它位置逐场指标见表。\n\n## A：21场主动态开发内验证\n']
    for eid in read(P4/'provenance/event_split.json')['validation']:atlas.append(f'### {eid}\n\n![{eid}六个内部断面](figures/A_validation_{eid}.png)\n')
    atlas.append('\n## B：20场固定控制先导\n')
    for eid in slope['order']:atlas.append(f'### {eid}（{"开发内验证" if eid in slope["validation"] else "训练"}）\n\n![{eid}流量和长度](figures/B_flood_{eid}.png)\n')
    (OUT/'FIGURE_ATLAS.md').write_text('\n'.join(atlas),encoding='utf-8')
    save(OUT/'provenance/frozen_protocol.json',dict(protocol_sha256=sha(OUT/'PROTOCOL.md'),B_frozen=slope,A_split=read(P4/'provenance/event_split.json'),B_HEC_quality_pass_count=int(qa.basic_quality_pass.sum()),
       B_main_solver_cases=40,B_extra_numerical_cases=4,A_training_scenarios=122,A_validation_scenarios=50,original_control_equation_preserved=True,posthoc_six_target_diagnostic=True,animations=False))
    save(OUT/'checkpoint.json',dict(status='requested_A_routing_and_B20_pilots_completed',updated=now(),A_events=70,A_scenarios=172,B_floods=20,B_HEC_main_completed=40,B_HEC_extra_completed=4,
       B_HEC_main_basic_quality_pass=40,raw_A_sources_unchanged=True,original30_test_results_opened=0,main_alpha_L_A=0,main_alpha_V_A=0,main_alpha_L_B=0,main_alpha_V_B=0,timestep_consistency_passed=False,all_requested_computations_complete=True))
    save(OUT/'delivery.json',dict(completed=True,completed_at=now(),report='REPORT.md',equations='EQUATIONS.md',A_events=70,A_scenarios=172,B_pilot_floods=20,B_main_HEC_cases=40,B_numerical_HEC_cases=4,
       static_main_figures=7,A_gallery_events=21,B_gallery_events=20,sealed_test_results_opened=0,animation_requested=False,method_independent_gain_established=False))
    state=read(WORK/'CURRENT_EXPERIMENT.json');state.update(updated=now(),analysis_stage='actual_routing_and_20paired_pilots_complete',analysis_latest=str(OUT),analysis_findings=str(OUT/'REPORT.md'),
       autonomous_pilot_B_root=str(B),autonomous_pilot_B_flood_events=20,autonomous_pilot_B_HEC_cases=40,autonomous_B_numerical_check_cases=4,autonomous_B_slope=.000075,
       experiment_A_train_events=49,experiment_A_internal_validation_events=21,experiment_B_train_events=14,experiment_B_internal_validation_events=6,
       analysis_candidate_status='actual_comparison_no_nonzero_gain_retained_timestep_limit_preserved',requested_computations_complete=True,animation_requested=False)
    save(WORK/'CURRENT_EXPERIMENT.json',state)
    print('Actual-computation report and complete static atlas ready',flush=True)
if __name__=='__main__':main()
