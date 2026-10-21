"""Assemble the readable research record from the verified numerical tables."""
import math,json,datetime,importlib.metadata
from pathlib import Path
import numpy as np,pandas as pd
from prepare import OUT,ROOT,read,save,frame,sha
from evaluate import configs

def table(rows,headers):
    return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(str(v) for v in r)+' |' for r in rows)+'\n'

def main():
    model=configs();primary=read(OUT/'fits/neural_gate.json')['chosen_seed_by_training_loss'];nn=f'NN{primary}'
    model=[c for c in model if c['type']!='neural' or c['primary']];labels={c['id']:c['label'] for c in model}
    E=pd.read_csv(OUT/'tables/event_summary.csv');E=E[E.subset=='validation'].set_index('model')
    R=pd.read_csv(OUT/'tables/paired_response_metrics.csv');RR=R[(R.subset=='validation')&(R.Lr_km==100)].groupby(['model','category'])[['NSE','RMSE']].mean()
    S=pd.read_csv(OUT/'tables/section_screen.csv').set_index('x_km');V=pd.read_csv(OUT/'tables/response_summary.csv').set_index('x_km')
    report=r'''# 内部正向流断面上的动态回水增量关系

本轮在100 km处建立了内部预测目标，并从小型可微分模块中提取、回嵌检验了三条解析关系。开发内验证支持在原有传播与蓄泄基底上加入与外部水位、河长和输入流量共同相关的附加状态。然而，完整洪水过程的高NSE并不保证边界响应正确：两个只由外部水位驱动的简洁候选仍不能稳定解释配对差值；固定参数的时间细分也保留了明显的基础模型离散敏感性。因此，本轮得到的是一小时离散条件下的候选关系，还不能作为定稿的普适RWEC物理公式。

## 内部目标的确定

筛查覆盖70场变水位、70场固定5 m参照和32次控制计算。每个工况从前48 h到完整恢复尾段，逐个保存值检查 `Q<0`，没有负值容差、裁剪、取绝对值或删除时段。172份原生HEC HDF的全部21个研究位置流量已与提取文件逐值核对一致。求解时间步为60 s，保存时间步为1 h；本轮证明的是全部已保存时刻无负值，不能证明小时之间未保存的求解时刻绝无倒流。

'''
    rows=[]
    for x in [20,50,90,100,150]:
        s=S.loc[x];rows.append([x,f'{s.min_q_m3s:.3f}',int(s.negative_cases),int(s.negative_samples),f'{s.negative_duration_case_sum_h:.2f}', '通过' if s.passed else '不通过'])
    report+=table(rows,['距外部边界/km','最小流量/(m³/s)','负值工况数','负值保存点数','负值时长合计/h','门槛'])
    vr=V.loc[100];aux=V.loc[150]
    report+=rf'''
表中时长由相邻小时保存点作线性过零估算，并在各工况间累加，不是额外观测或单场最大持续时间。初选20、50、100 km中只有100 km合格；扩展筛查表显示90 km仍有负值，100—190 km合格。100 km保留上游100 km演算空间，同时70场配对的最大水位增量为{vr.max_stage_increment_m_min:.3f}—{vr.max_stage_increment_m_max:.3f} m，中位数{vr.max_stage_increment_m_median:.3f} m；配对最大流量差绝对值为{vr.max_abs_q_difference_m3s_min:.2f}—{vr.max_abs_q_difference_m3s_max:.2f} m³/s。150 km作为水动力辅助位置，最大水位增量中位数为{aux.max_stage_increment_m_median:.3f} m，剩余演算长度仅50 km。本轮不在辅助位置追加率定或推断空间衰减函数。

选点只依赖水动力方向、配对响应及剩余河长，不依赖RWEC拟合分数。主目标由此固定为 $x_c=100\ \mathrm{{km}}$。输入位置 $x_u=110,120,\ldots,200\ \mathrm{{km}}$，有效演算长度 $L_r=x_u-x_c$。每个工况使用自身的 $I(t)=Q_{{\mathrm{{HEC}}}}(x_u,t)$，目标为100 km处的流量。已知外部水位仍是0 km处的 $H_d(t)$，既不是100 km局部水位，也不被当成实际大坝位置。内部输入本身可能已经带有回水影响，短河段中基础模型已能继承部分响应；较长河段需要额外解释中间蓄泄。

## 共享基底与物理审查

70场开发洪水按预先冻结的事件随机排列划为49场训练、21场开发内验证，同一事件的所有河长、边界配对和补充控制始终同组。30场原独立测试完全封存。控制计算不参与训练，其四个背景中只有E006属于本轮验证组，因此该组的8次控制计算仍只有1场未见背景洪水，不能按8个独立洪水声称泛化。

主对比首先仅用49场训练固定水位参照率定每种河长的原三参数T-RWEC，然后所有增量候选共用并冻结这30个基础参数。原递推采用延迟后的步平均输入与平衡初始储量：

$$V_{{I,t}}=\bar I_{{\tau,t}}\Delta t,\quad S_{{W,t}}=S_{{t-1}}+V_{{I,t}},\quad S_{{E,t}}=S_{{t-1}}+kV_{{I,t}},$$
$$S_t=(S_{{W,t}}^{{-n}}+S_{{E,t}}^{{-n}})^{{-1/n}},\qquad Q_{{R,t}}=\frac{{S_{{t-1}}+V_{{I,t}}-S_t}}{{\Delta t}}.$$

参数 $k,n,\tau$ 的版本与上一轮已核对论文的递推一致，完整基础源码逐字复制后独立导入。没有采用未来目标出流、真实内部水位或HEC储量初始化。基础参数在10—100 km间的 $k$ 为0.0215—0.2465、$n$ 为3.645—14.295、$\tau$ 为0.904—12.812 h，均未命中边界；逐河长数值与近边界标记保留于表中。

附加状态 $B$ 定义为相对于基础解的储量异常，不能再把同一体积加到 $S_E$。出流由 $Q=Q_R-\Delta B/\Delta t$ 给出；再计入传播在途体积 $V_T$，有

$$\Delta(S+B+V_T)=(\bar I-Q)\Delta t,\qquad\Delta V_T=(\bar I-\bar I_\tau)\Delta t.$$

这检查代码的水量一致性，不等于验证该模型状态就是HEC的绝对几何储量。P1的即时体积交换与P2的一阶响应原式，都不能对任意强迫保证正向出流。本轮先书面审查，再采用明确标记的状态约束修正版：先限制附加蓄入体积不超过本步基础可供水量、释放不超过已有附加体积，再由体积差推得出流；原式未经约束的最小出流另存，不在输出端裁剪掩盖问题。P2对步内线性平衡容量采用解析积分。现有全部172个工况的所有拟合候选中，供水约束均未激活，因此该约束没有制造当前精度改善。

名义固定5 m边界在HEC保存为4.99993896484375 m。容量使用 $[H_d-5]_+$ 保持名义固定边界的零增量退化，差别已在PHYSICS_REVIEW.md中明确记录；没有改动任何原始水位或流量，也没有放宽负流量筛查。真正低于5 m的边界不属于本轮容量关系的适用范围。

P1得到 $\alpha={read(OUT/'fits/P1.json')['alpha']:.6f}$；P2得到 $\alpha={read(OUT/'fits/P2.json')['alpha']:.6f}$、$T_b={read(OUT/'fits/P2.json')['tb_h']:.3f}$ h。100 km位置固定，$\alpha$ 只表示该位置的有效响应，不能解释为已识别的空间衰减规律，$T_b$ 也不是波速。训练损失剖面提供条件参数检查，而不是正式置信区间。

共享固定基底的主对比是机制消融实验，并非所有模型分别完全联合优化的竞赛。看到主对比后另加了一次有限的RWEC-0均衡率定诊断，在同样49场训练洪水两类边界上优化其30个参数，采用相同范围和约2000次每河长的预算；不改变其他候选、事件划分和符号发现门槛。它的结果与原基准并列，限定“改善来自较弱基准”的解释，不能把这项事后诊断称为事前检验。

## 从可微分模块到解析关系

本轮借鉴[Hu等关于可微分混合模型用于水文知识发现的研究](https://doi.org/10.1016/j.jhydrol.2026.136147)，先诊断增量蓄泄关系，再提取可运行解析式。实现为3输入、4个tanh隐单元、1个sigmoid输出的小网络，21个网络参数和1个共享响应时间；基础传播与守恒状态递推保留。输入均无量纲且可在演算时获得：

$$l=\frac{{L_r}}{{100000\ \mathrm{{m}}}},\qquad h=\frac{{[H_d-5\ \mathrm{{m}}]_+}}{{6\ \mathrm{{m}}}},\qquad i=\frac{{\bar I_\tau}}{{6000\ \mathrm{{m^3/s}}}}.$$

内部目标位置固定，未把常数 $x_c$ 加为无从识别的网络变量；也没有输入HEC目标出流、内部真实水位或真实储量。网络学习容量系数 $c=2\sigma(g)$，并给 $B_{{eq}}=500L_r[H_d-5]_+c$，仍经附加状态方程和守恒递推产生出流。单位由河宽、河长、水位和时间显式恢复。

两个种子分别在21场未训练洪水中改善20场。相对固定基底的归一化平方误差降低约74.7%，事件平均NSE增量约0.00483，中位数约0.00143；按事件重采样的均值95%区间下限约0.00157。训练前已把不适合接近NSE上限的绝对0.01门槛改为正的事件改善区间、至少5%的相对误差降幅及至少60%的洪水改善，原协议和修改记录均保留。门槛修改发生在任何网络训练与验证模型评分之前。该区间只描述当前开发内分组，不是外部泛化保证。

依[PySR官方实现](https://github.com/astroautomata/PySR)，符号蒸馏只从49场训练事件抽取3920个活跃水位时刻的网络容量logit值，每事件最多80点；未使用验证目标。限定加、减、乘、除和复杂度上限12，按训练蒸馏损失在复杂度上限4、8、12各取一式，未按验证分数重新搜索。三式均回嵌独立递推计算完整过程。下式河宽为500 m，$L_r$与水位增量以m计，$B$以m³计，出流与体积差换算中的时间以秒计，列出的响应时间小时值在运行时换算：

$$g_1=h-3.5332525232272065,$$
$$g_2=h-i-l-2.6874291614539905,$$
$$g_3=h(2.1687527388137-i)-(i+l)-3.054277891787369,$$
$$c_j=\frac{{2}}{{1+\exp(-g_j)}},\qquad B_{{eq,j}}=500L_r[H_d-5]_+c_j,$$
$$\frac{{dB_j}}{{dt}}=\frac{{B_{{eq,j}}-B_j}}{{T_b}},\qquad T_b={math.exp(read(OUT/'fits'/f'neural_{primary}.json')['theta'][21]):.9f}\ \mathrm{{h}},\qquad Q_j=Q_R-\frac{{\Delta B_j}}{{\Delta t}}.$$

这些是共享物理外壳内的容量关系，PySR没有发现整套RWEC。sigmoid界限、守恒关系和状态约束是事先指定结构，响应时间继承自训练网络而未用验证重新率定。SR1仅随水位变化；SR2增加输入与河长；SR3进一步含水位与输入的交互。运行常数和结构复杂度不能代表整个发现过程的统计自由度，教师网络与结构搜索的选择成本仍存在。复杂度3、7、11指logit表达式的搜索复杂度，不是完整模型的参数数。

## 完整过程与边界差值

下表以洪水为单位，先在每场两类边界与10种河长中平均，再汇总21场验证洪水。它只反映100 km内部目标，不与上一轮0 km目标的NSE作直接高低比较。KGE和RMSE列为相同事件汇总的平均值；所有指标的平均、中位、最小、最大及逐场、逐工况、逐河长结果均在tables中。

'''
    rows=[]
    for c in model:
        r=E.loc[c['id']];rows.append([c['label'],f'{r.NSE_mean:.6f}',f'{r.NSE_median:.6f}',f'{r.NSE_min:.6f}',f'{r.NSE_max:.6f}',f'{r.KGE_mean:.6f}',f'{r.RMSE_mean:.3f}'])
    report+=table(rows,['模型','NSE均值','NSE中位','NSE最小','NSE最大','KGE均值','RMSE均值/(m³/s)'])
    report+='\n边界差值是更严格的机制检验。输入200 km时，各边界工况的上游输入相同，基础模型几乎给不出配对流量差。以下统计在有效演算长度100 km处，按21场主配对及验证背景E006的控制对比分别列出；后一组只有1场背景洪水，不能据其多组对比数增加独立样本量。\n\n'
    rows=[]
    for c in model:
        a=RR.loc[(c['id'],'main_pair')];b=RR.loc[(c['id'],'factorial_control')]
        rows.append([c['label'],f'{a.NSE:.4f}',f'{a.RMSE:.2f}',f'{b.NSE:.4f}',f'{b.RMSE:.2f}'])
    report+=table(rows,['模型','主配对差值NSE均值','差值RMSE/(m³/s)','E006控制差值NSE均值','差值RMSE/(m³/s)'])
    report+=r'''
P1、P2和经验项的完整过程分数有所改善，但主配对差值平均NSE仍为负；不同事件的差值方差不同，所以平均RMSE略降并不保证等事件平均NSE为正。小网络明显改善边界差值，SR2对主配对差值优于SR3，而SR3在全程精度与E006控制响应上更接近网络。离线蒸馏误差最小的式子未在全部机制指标上同时最优，说明回嵌递推和差值检验不可省略。SR1的简洁性不足以弥补其边界响应损失，SR2和SR3宜作为有取舍的开发候选并列保留，不能只按整体NSE宣布唯一优胜公式。

E094的小网络事件平均NSE略有退化，完整过程与逐时误差保留于图2、图5和图7；未删除该事件，也没有事后换断面。控制背景E007中的强幅度变化仍有明显抑流幅度和相位误差。部分零响应基准的差值KGE无定义，表内保留为空值，不以零或优良数值替代。

## 守恒、方向与离散限制

'''
    B=pd.read_csv(OUT/'tables/balance.csv');T=pd.read_csv(OUT/'tables/timestep_sensitivity.csv')
    report+=rf'''所有已回嵌模型在172个开发相关工况、10种演算长度中均无负预测值，供水约束激活次数为零，包含在途体积的最大逐步水量残差为{B.physical_control_balance_max_m3.max():.3g} m³。模型运行函数不接收真实目标作为输入，改变未来边界不会改变既有预测前缀；解析梯度已通过中央差分核对，供水约束另用强迫测试验证。以上结论只覆盖本轮正向流与名义不低于5 m的外部边界，不能外推至原0 km的倒流情况。

固定参数的时间细分使用同一分段线性输入，不补造更高频HEC真值。0.5和0.25 h输出在每个原小时内聚合后比较。4个背景主事件、3种河长的差异如下；相邻细分差异下降表明有数值细分趋势，但1 h率定结果仍不能直接跨步长移用。

'''
    rows=[]
    for c in model:
        u=T[T.model==c['id']]
        v1=u[(u.comparison=='vs_1h')&(u.dt_h==.5)].RMSE_vs_1h_m3s
        v2=u[(u.comparison=='vs_1h')&(u.dt_h==.25)].RMSE_vs_1h_m3s
        v3=u[u.comparison=='0.25h_vs_0.5h'].RMSE_vs_1h_m3s
        rows.append([c['label'],f'{v1.mean():.2f}',f'{v2.mean():.2f}',f'{v3.mean():.2f}',f'{v2.min():.2f}—{v2.max():.2f}'])
    report+=table(rows,['模型','0.5 h对1 h均值','0.25 h对1 h均值','0.25 h对0.5 h均值','0.25 h对1 h范围'])
    report+=r'''
单位均为出流差异RMSE（m³/s），完整的平均、中位、最小、最大见timestep_summary.csv。P1、P2的差异几乎完全继承基础解，小网络和解析容量的良好守恒也未消除基础递推的步长敏感性。本轮没有为了提高分数重新按步长率定k、n或τ；高一小时NSE与跨步长稳定性应分开判断。共享的集中响应也尚未分辨外部扰动的空间传播、波反射和输入流量耦合各自的贡献，不能把12.17 h解释为独立测得的水动力传播时间。

## 同步动画与研究判断

'''
    a=read(OUT/'animation/animation_summary.json')
    report+=rf'''E007动画使用同一事件、同一时刻的全部101个既有HEC水面线，配合200 km输入、0 km外部水位与100 km真实流量。视频为1920×1080、12帧每秒、{a['video_duration_s']:.2f} s，主体每1 s对应12 h模拟时间，末尾3 s为说明；GIF为同一过程的预览。100 km最大水位增量为{a['max_stage_increment_at_100km_m']:.3f} m，变动水位减固定水位流量在{a['min_difference_time_h']:.0f} h最小，为{a['min_var_minus_ref_flow_m3s']:.2f} m³/s，在{a['max_difference_time_h']:.0f} h最大，为{a['max_var_minus_ref_flow_m3s']:.2f} m³/s。抑流后的相对增流在外部水位刚开始消落时已经出现，而非只在边界完全回到5 m后才恢复。这提示关系不仅要响应水位变化，还需要表达洪水输入与附加状态之间的耦合。动画不提供额外验证样本，也没有把水面调整前缘当成水流速度。

本轮最有价值的结果是把“过程能拟合”与“边界响应能解释”分开，并得到可独立运行、可检查的SR2/SR3候选。它们在当前正向内部断面与一小时离散条件下值得继续审查；时间离散可迁移性、跨目标位置传播以及外部独立洪水泛化仍未解决。下一阶段应由研究者先审阅这些边界与失败，再决定是否冻结关系进入原30场独立测试。本轮在此停止，不增加HEC工况、不继续RWEC-2、不扩展网络或反复调整事件划分。

原始172份输入提取文件的哈希与修改时间均在交付检查中复核，原生HDF仅以只读方式打开。运行环境的Windows长路径与下载连接问题及修复日志均保留，PySR最终成功完成，未将依赖失败当成科学失败或用人工式子替代搜索。研究正文中的公式为可编辑LaTeX源码。完整来源、复算入口、冻结方案及调整记录见REPRODUCE.md、PROTOCOL.md和provenance目录。
'''
    (OUT/'REPORT.md').write_text(report,encoding='utf-8')
    # Complete distribution tables, preserving raw rows and event grouping.
    ts=T.groupby(['model','comparison','dt_h']).RMSE_vs_1h_m3s.agg(['mean','median','min','max']).reset_index()
    frame(OUT/'tables/timestep_summary.csv',ts)
    P=pd.read_csv(OUT/'tables/phase_metrics.csv');P=P[P.role=='main']
    g=P.groupby(['model','subset','kind','Lr_km','phase'])[['RMSE','bias_m3s']].agg(['mean','median','min','max']);g.columns=['_'.join(x) for x in g.columns]
    g.reset_index().to_csv(OUT/'tables/phase_summary.csv',index=False,encoding='utf-8-sig')
    gains=pd.read_csv(OUT/'tables/event_improvements.csv');gg=gains.groupby(['model','subset'])[['NSE_gain','RMSE_change']].agg(['mean','median','min','max']);gg.columns=['_'.join(c) for c in gg.columns]
    gg.reset_index().to_csv(OUT/'tables/event_improvement_summary.csv',index=False,encoding='utf-8-sig')
    meta=read(OUT/'data/metadata.json');C=pd.read_csv(OUT/'tables/case_metrics.csv');C=C[C.role!='main']
    cols=['NSE','KGE','RMSE'];cg=C.groupby(['model','event_id','subset'])[cols].mean().reset_index();frame(OUT/'tables/control_background_event_metrics.csv',cg)
    cs=cg.groupby(['model','subset'])[cols].agg(['mean','median','min','max']);cs.columns=['_'.join(x) for x in cs.columns]
    cs.reset_index().to_csv(OUT/'tables/control_background_summary.csv',index=False,encoding='utf-8-sig')
    # Preserve learned weights and costs as raw evidence, not physical parameters.
    nr=[];cost=[]
    for seed in [20261021,20261022]:
        rec=read(OUT/'fits'/f'neural_{seed}.json')
        nr.extend([dict(seed=seed,parameter=j,value=v,physical_interpretation='none: network weight' if j<21 else 'log shared response time in hours') for j,v in enumerate(rec['theta'])])
        cost.append(dict(model=f'NN{seed}',objective_evaluations=rec['nfev'],iterations=rec['nit'],elapsed_s=rec['elapsed_s'],includes_shared_core=False))
    frame(OUT/'tables/neural_parameters.csv',nr);frame(OUT/'tables/neural_training_costs.csv',cost)
    versions={x:importlib.metadata.version(x) for x in ['numpy','scipy','numba','pandas','matplotlib','h5py','Pillow','sympy']}
    save(OUT/'provenance/runtime_versions.json',dict(python=__import__('sys').version,packages=versions,
        optional=read(OUT/'fits/symbolic_candidates.json'),font_files=['C:/Windows/Fonts/times.ttf','C:/Windows/Fonts/simsun.ttc']))
    save(OUT/'summary.json',dict(main_xc_km=100,calibration_events=49,development_validation_events=21,sealed_test_events_used=0,
        cases_screened=172,pysr_completed=True,symbolic_candidates=read(OUT/'fits/symbolic_candidates.json')['candidates'],
        no_single_candidate_dominates_all_diagnostics=True,timestep_portability_passed=False,
        main_validation_mean_nse={c['id']:float(E.loc[c['id']].NSE_mean) for c in model},
        main_pair_response_nse_L100={c['id']:float(RR.loc[(c['id'],'main_pair')].NSE) for c in model},
        stopping_for_researcher_review=True))
    print('研究报告与完整汇总表已写入',flush=True)

if __name__=='__main__':main()
