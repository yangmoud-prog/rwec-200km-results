# 旧河道 A 新公式验证：文件入口

本目录独立于前十阶段，仅使用原旧河道 A 的 70 场开发洪水、140 条主固定／变水位结果。沿用原 49／21 划分，不读取原 30 场封存测试，不重建河槽、不重跑 HEC。`REPORT.md` 为本轮中文正文，`SUMMARY.json` 为结果摘要。

| 结果 | 文件 |
|---|---|
| 冻结方案、预算、范围及窗口 | `PROTOCOL.md`，`provenance/frozen_protocol.json` |
| 全部源数据散列及质量 | `SOURCE_MANIFEST.json`，`tables/HEC_source_quality.csv` |
| 全部逐场指标，本地完整表 | `tables/event_metrics.csv` |
| 阅读包内逐场指标 | `tables/event_metrics_full.csv`，`tables/event_metrics_nonnegative_hour_endpoints.csv`，`tables/event_metrics_after_last_negative_plus12h.csv` |
| 21 场主变水位验证指标 | `tables/validation_dynamic_event_metrics.csv` |
| 逐断面均值、中位数、最小值、最大值 | `tables/section_statistics.csv` |
| 同场同断面改善 | `tables/paired_improvements.csv` |
| 逐事件平均改善及完整胜负 | `tables/event_mean_gains.csv`，`tables/event_gain_summary.csv` |
| 负流量、删除时段及全程非负标记 | `tables/negative_flow_windows.csv`，`tables/validation_negative_summary.csv` |
| 整场无负值事件的补充评价 | `tables/whole_event_nonnegative_metrics.csv` |
| 原式、新式、广义式、非线性马斯京根参数 | `PARAMETERS_A_original.json`，`PARAMETERS_A_new.json`，`PARAMETERS_A_general.json`，`PARAMETERS_A_NLM.json` |
| 物理量和记忆参数解释 | `tables/parameters_interpreted.csv` |
| 历史冻结原 RWEC 的同口径参考 | `tables/historical_frozen_original_event_metrics.csv`，`tables/historical_frozen_original_means.csv` |
| 固定参数的时段敏感性 | `tables/fixed_parameter_dt_sensitivity.csv` |
| 有效储量及移动边界以外的接口数值体积 | `tables/storage_diagnostics.csv`，`tables/cascade_mass_balance.csv` |
| 恒定流、求解状态和因果性 | `tables/final_parameter_steady_checks.csv`，`tables/solver_quality.csv`，`provenance/causality_check.json` |

断面 x 从出口向上游计量，距入口为 200−x km。上游输入固定在 x=200 km。目标为 x=120、100、80、60、40、20、0 km。实际演算长度相应为 80、100、120、140、160、180、200 km。

评分采用相同逐小时区间平均流量。`nonnegative_hour_endpoints` 仅评价两个真实保存节点均非负的区间，三个模型和非线性马斯京根共用同一掩码。`after_last_negative_plus12h` 仅是倒流结束后再缓冲 12 h 的连续窗口，不称为水力稳态。模型始终完整递推，掩码只在评分时使用。参考工况标记 `fixed5_same_forecast`，使用相同入流模型预测再与固定边界真值比较，并未单独在固定边界优化参数。

全程无负值事件只作诊断分组，未取代全部 21 场验证，也未用于更改训练。开发内验证已经在前期研究中查看，本轮不是首次独立测试。

核心 `formulas.py` 与前一阶段公式内核完全一致。`prepare_A.py` 提取旧 A 数据并冻结方案，`fit_A.py` 进行逐段独立率定，`run_A.py` 以四个本地进程协调四种模型，`evaluate_A.py` 在全部参数冻结后统一评价。`figures_A.py` 生成中文静态图，每场验证均保留 100 km、20 km、出口三个断面的完整过程图。`report_A.py` 从原始表格生成正文，`final_checks.py` 检查源数据、预算、参数、评分窗口和导出图。

重新复现应先复制至独立目录并核对根路径；率定程序会复用本目录已有训练缓存，不应直接覆盖冻结结果。参数已在训练前允许原旧河道的宽范围，未根据验证表现修改。广义式的嵌套起点仅来自本轮原式／新式训练结果，所有目标函数调用均计入相同预算。

`EQUATIONS.md` 和 `EQUATIONS_EDITABLE.tex` 保留可编辑 LaTeX，`equations/main.tex` 是同源可编译文档。字体为英文／数字 Times New Roman、中文宋体。`figures/` 提供 PNG、SVG、PDF，本轮无动画。

完整逐步预测位于本地 `predictions/`，每文件包含四分之一小时节点流量、共同小时均值、有效状态、等效体积和接口重采样体积。原始 HEC 数据位于旧 A 的 `results/`，均保持只读。阅读仓库仅同步轻量正文、公式、代码、表格、JSON 和 PNG，超出表格体积上限的文件保留本地，清单会说明省略项。
