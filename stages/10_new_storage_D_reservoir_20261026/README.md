# 第十阶段：文件与复现入口

本目录为独立第十阶段。前九阶段原始成果未改写，原 30 场封存测试未启用。研究结论以 `REPORT.md`、`SUMMARY.json` 和逐场原始指标为准。自主联合模型尚未通过水面／储量／通量物理闭合，不存在已验证的耦合预测结果或关闭反馈消融分数。

`REPORT.md` 给出中文科学叙述、全部对照条件及失败解释。`EQUATIONS.md` 为可编辑 LaTeX 公式；本地 `equations/main.tex` 与根目录 `EQUATIONS_EDITABLE.tex` 为同源 TeX 文档，已用桌面内置编译器检查。`figures/` 保存中文 PNG、可编辑文字 SVG 和 PDF；英文、数字采用 Times New Roman，中文采用宋体。本轮无动画。

| 证据 | 文件 |
|---|---|
| 已有 B、C 公式逐场结果 | `tables/formula_event_metrics.csv` |
| B、C 逐断面均值、中位数、最小值、最大值 | `tables/formula_section_statistics.csv` |
| B、C 逐场相对原式改善 | `tables/formula_event_improvements.csv` |
| B、C 参数与局地记忆解释 | `tables/formula_parameters_interpreted.csv`，`PARAMETERS_B.json`，`PARAMETERS_C.json` |
| D 的冻结样本与设计 | `D_FROZEN_DESIGN.json`，`D_PHYSICAL_SELECTION.json`，`D_SELECTION_REPLAY.json` |
| D 的 80 次正式 HEC 质量检查 | `tables/D_formal_numerical_quality.csv`，`tables/D_independent_numerical_audit.csv` |
| HEC 时空敏感性与先导结果 | `tables/D_sensitivity_quality.csv`，`tables/D_pilot_pairs.csv` |
| D 水动力两因素配对差值 | `tables/D_factorial_physical_discharge.csv`，`tables/D_factorial_physical_storage.csv` |
| D 全模型逐场逐断面指标 | `tables/D_event_metrics.csv` |
| D 指标均值、中位数、最小值、最大值 | `tables/D_section_statistics.csv` |
| D 同场同断面原始改善 | `tables/D_paired_improvements.csv` |
| D 每场平均改善及全部胜负 | `tables/D_event_mean_gains.csv`，`tables/D_event_gain_summary.csv` |
| D 所有子段全部参数 | `tables/D_parameters_interpreted.csv`，`PARAMETERS_D??_*.json` |
| 蓄入与释放阶段误差 | `tables/D_storage_release_phase_errors.csv` |
| 固定参数时段敏感性 | `tables/D_fixed_parameter_dt_sensitivity.csv` |
| 有效储量异常与变化率 | `tables/D_storage_state_event_metrics.csv`，`tables/D_storage_state_statistics.csv` |
| 不同时段子段的整体水量衔接 | `tables/D_cascade_mass_balance.csv`，`STORAGE_AUDIT.json` |
| 相同真实移动评价位置的基准流量 | `tables/D_same_endpoint_baseline_metrics.csv` |
| 真实末端已知条件下的物理闭合 | `tables/D_oracle_physical_closure.csv` |
| 失败与纠错记录 | `FAILURES.md`，`provenance/invalid_time_alignment_attempt/` |
| 因果性、恒定流和非负性检查 | `provenance/causality_check.json`，`provenance/final_model_checks.json`，`tables/D_final_parameter_steady_checks.csv` |

DUF、DUR、DCF、DCR 依次表示均匀固定、均匀动态库区、复杂固定、复杂动态库区。断面坐标 x 从下游出口向上游增加；距上游入口为 200−x km。目标为距出口 120、100、80、60、40、20、0 km。流量统一以向下游为正。

移动评价位置表对四个基准统一使用每小时两端末端位置的平均值，再插值区间平均流量。理想物理闭合则使用各时刻自身末端处真实流量的梯形积分；两种时间离散诊断分别报告，不跨表混作模型排名。

`storage_audit.py` 记录全部冻结模型的有效储量，重算流量与原预测完全一致。原生子段局部守恒与整条串联河道守恒分开核验；子段时段不同时的流量重采样会产生小幅累计接口体积偏差，最大约占入流洪水增量体积 0.15%。该项单列为数值记账状态，不混入物理回水储量，也不用于回改已报告的预测精度。

本地完整工程为 `F:/YDX/200km/variants/compound_reservoir_D_20261026`。`models/` 中每个工况含 HEC `.prj/.g01/.p01/.u01`、原生 HDF、计算日志、输入散列、断面设计、提取时序、质量记录；动态工况另含与固定边界配对的 `paired.npz`。目录 `D001_complex_reservoir_v3` 等为正式工况，`DP*` 为先导与延长恢复试验，`DS*` 为数值敏感性试验，失败的最早 `v1` 工程也保留。20 场新洪水各有四个工况，共 80 个正式工程。文件存在与求解返回成功均不单独视为质量合格。

本地分析时序在本目录 `predictions/`，包括全部基准预测与每个理想条件水面候选的储量、扫掠和末端实际流量序列。真实末端只在理想物理诊断及统一移动评价位置时读取，未驱动任何宣称自主的预测结果。`legacy_models.py` 是保留的旧求解器副本，第十阶段仅调用其中 `muskingum` 与 `at`；该文件的其他历史函数未用于本轮新公式或水面候选。

复现环境为 Windows、HEC-RAS 6.3.1，Python 使用 NumPy、SciPy、pandas、Numba、h5py、matplotlib、pywin32、psutil。实际执行的 Python 为 `D:/04-Tools/Anaconda/python.exe`。HEC 原生工程不依赖另装新版 HEC；官方算例仅作为输入结构参考。

公式内核为 `formulas.py`，检查为 `test_formulas.py`。已有 B、C 的训练、等预算训练精化和一次锁定验证分别为 `fit_formulas.py`、`refine_training.py`、`evaluate_formulas.py`。D 的逐段公平率定及完整窗口评价为 `formal_baselines.py`，四种条件由 `baseline_coordinator.py` 并发协调；`PROTOCOL_D??.json` 记录完整预算和范围。已存在的率定文件会被用于续算，完整重新复现应复制到独立目录或先明确归档，避免覆盖当前结果。

河道 D 原生工程的构建、计算、提取脚本实际位于 D 目录的 `scripts/`。阅读包根目录 `D_source_*.py` 是对应源码副本，阅读用途保留此前缀；如需在独立副本运行，恢复原文件名并放回同一 `scripts/` 目录，核对根路径后执行。固定设计已记录在 `D_FROZEN_DESIGN.json`，不要再次随机生成样本或用验证精度修改设计。

`aggregate_formal.py` 汇总逐场指标、相同位置评价和因果性扰动检查；`audit_final_models.py` 复核全部已率定参数的恒定流与非负性；`factorial_analysis.py` 计算同空间断面的两因素水动力差。`oracle_closure.py` 仅进行真实末端已知条件下的几何与通量证伪，不能作为自主预报入口。`figures10.py`、`formal_figures.py` 和 `geometry_figure.py` 生成静态图件，`deliver_stage10.py` 从冻结结果生成正文。

阅读仓库只保存轻量正文、LaTeX 公式、源码、表格、JSON 和 PNG。原始 HEC/HDF、完整 NPZ、较大表格、PDF/SVG 以及无效试算原档留在本地，阅读包不是可脱离原始数据独立重算的完整水动力副本。实际远端状态以阅读仓库 `.local/last_upload.json` 的读回核验为准。
