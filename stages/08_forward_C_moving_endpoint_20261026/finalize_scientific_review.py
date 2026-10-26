from common import *
report=P/'REPORT.md';txt=report.read_text(encoding='utf-8')
lead='''
**本轮判断是“存在首轮预测增益，但还不能认定新模型已稳定优于传统方法”。** 自主候选在原率定时段下改善了出口附近的完整过程与抑流／释放误差，但固定参数改变计算时段后精度明显下降；理想末端还会在阈值消失时制造释放尖峰。这两项缺陷必须先解决，不能仅凭平均NSE接近0.99扩大成功结论。
'''
txt=txt.replace('## 河槽选择及数值证据',lead+'\n## 河槽选择及数值证据',1)
detail='''
六场验证、出口与20 km两断面合并统计，自主候选的抑流阶段RMSE为151.35 m³/s，原RWEC为205.17 m³/s，非线性马斯京根为157.68 m³/s；释放阶段分别为78.39、124.89和133.51 m³/s。这说明新增暂存机制对局部过程有作用，不能把结果完全归结为洪峰平移。释放阶段各模型NSE仍很低甚至显著为负，部分原因是该短窗口真实流量方差很小；因此同时保留RMSE和全过程图，不用阶段NSE与全程NSE直接比较。

失败同样具体：理想末端候选在E032出口退水期高估591.79 m³/s，恰逢0.10 m阈值末端由约32.48 km跳至0；E051对应高估552.47 m³/s。真实流量此时仍平滑，故这种峰值是阈值长度与暂存体积闭合产生的伪释放。未来改进应检验连续影响强度或水面形状状态，不能直接平滑真实评价数据来掩盖缺陷。

参数固定时，自主候选在出口与20 km处的平均NSE由0.9718降至时段减半的0.8900、加倍的0.7713。原RWEC对应0.9536、0.8769、0.7932；非线性马斯京根约0.9583且对这次时段变化几乎不变。RWEC系列多个最优参数接近k搜索下界，时间离散与经验参数之间仍存在强补偿。新增机制的初步增益不能视作已经具有时间尺度一致性；本轮保留原结果，不利用验证成绩继续修改公式或重新筛选河槽。

在共同真实末端位置上，非线性马斯京根的平均NSE为0.9908，自主候选为0.9896，原RWEC为0.9886。因此新候选没有在两个预测目标上全面胜出。下一轮若获授权，应先修正时段依赖和阈值伪释放，再使用新的、预先约定的开发检查；30场封存测试仍继续封存。
'''
txt=txt.replace('## 数据与复核入口',detail+'\n## 数据与复核入口',1);report.write_text(txt,encoding='utf-8')
r=pd.read_csv(P/'tables/event_raw_improvements.csv');r=r[r.subset=='validation'];s=pd.read_csv(P/'tables/endpoint_error_window_check.csv')
save(P/'SUMMARY.json',dict(conclusion='first-pass gain, not yet a timestep-robust superior model',HEC_runs=69,development_events=20,train_events=14,validation_events=6,sealed_tests_opened=0,
 original_mean_NSE=.9811157,autonomous_mean_NSE=.9868257,nonlinear_Muskingum_mean_NSE=.983081,
 improved_fixed_event_sections=int((r.MC_autonomous_minus_RWEC0>0).sum()),total_fixed_event_sections=len(r),
 endpoint_MAE_full_km=float(s.MAE_full_km.mean()),endpoint_MAE_active_km=float(s.MAE_0_360h_km.mean()),
 retained_failures=['0.05per mille CP01 reverse flow','initial overly strict per-event hydraulic response gate','threshold disappearance spurious release','RWEC timestep sensitivity','parameters at search bounds'],
 interpretation='Positive-flow C does not force original RWEC to perform poorly. Modest local gain, unresolved physical closure and timestep invariance.'))
save(P/'provenance/causality_check.json',dict(autonomous_query_inputs=['known upstream forecast','outlet history -48h through0','fixed geometry','initial steady state'],future_HEC_H_Q_storage_length_used_by_autonomous=False,ideal_true_length_separately_labeled=True,training_endpoint_forecast_event_leave_one_out=True,validation_library_training_events_only=True,all_models_same_evaluation_positions=True,no_sealed_test_read=True,models_parameters_frozen_before_validation=True,river_selected_before_routing_calibration=True))
print('Scientific limitations and failures added')
