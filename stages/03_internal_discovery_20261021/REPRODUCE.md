# 复算与来源

本目录是第二阶段的独立分析交付，活动工程由工作区根目录 CURRENT_EXPERIMENT.json 定位。原始结果只读；本轮既不调用HEC求解器，也不读取封存的30场测试结果。

## 从冻结参数复算

在本目录打开PowerShell，使用Python 3.12执行：

```powershell
$env:PYTHONIOENCODING = 'utf-8'
$env:OPENBLAS_NUM_THREADS = '2'
$env:OMP_NUM_THREADS = '2'
py -3.12 replay.py
```

默认入口读取已保存的拟合参数，重新运行完整递推、配对差值、步长检查、因果性检查、中文静态图、动画、报告及交付检查。不会重新训练或重新搜索表达式，也不会改动原始HEC文件。静态图由Python绘制，英文字母、数字、单位和符号使用Times New Roman，中文使用SimSun；系统必须有这两种字体。报告中的公式直接保留LaTeX源码。

精确软件版本见 provenance/runtime_versions.json；基础依赖在当前Python环境中，PySR和视频编码依赖独立保存在本目录 .deps。Julia运行时与包缓存位于工作区根目录 .runtime，以避免Windows长路径失败。optional_requirements.txt列出本次成功使用的独立依赖版本；不应将未约束的新numpy/scipy版本安装到现有环境后声称精确复现。

## 从头重复率定与发现

为保留原有结果，应将本目录的Python脚本、PROTOCOL.md、PHYSICS_REVIEW.md、ROBUSTNESS_ADDENDUM.md、FIGURE_LEGENDS.md、animation/make_animation.py及provenance/USER_SPEC.md复制到同一 analysis 目录下的新空阶段文件夹；保留相同层级，使 ROOT 指向同一活动工程。复制当前 .deps 或按下述独立依赖安装；使用工作区现有 .runtime，不能向旧研究档案写入。然后依次执行：

```powershell
py -3.12 prepare.py
py -3.12 verify.py native
py -3.12 fit.py base
py -3.12 fit.py increments
py -3.12 fit.py neural
# 只有 fits/neural_gate.json 的 passed 为 true 时执行下一行
py -3.12 symbolic.py
py -3.12 fit.py jointbase
py -3.12 replay.py
```

基础率定用49场训练固定水位参照，所有主候选冻结这组基础参数；增量拟合使用相同49场洪水两类边界，控制计算不训练。均衡率定基准是单次额外稳健性诊断，不改变门槛或重新选择候选。程序会复用已经存在的拟合检查点；从头率定须在新空目录中进行，不删除或覆盖失败记录来伪装独立重复。

神经模块是NumPy/Numba实现的3—4—1小网络，以解析反向链式求导和L-BFGS-B优化，中央差分核对见gradient_check.json。未依赖PyTorch。原始T-RWEC代码逐字复制为base_core.py，只改变导入方式以避免动态模块的Numba跨进程缓存问题，没有改变递推。两个固定种子分别最多400次迭代。若门槛失败，应报告该失败并停止符号发现，不能绕过门槛硬跑搜索。

PySR使用2.7.2、JuliaCall 0.9.36及Julia 1.11.9；只学习训练事件网络容量logit，不直接拟合验证HEC目标。搜索算子、样本索引、Pareto前沿与完整表达式分别保存在symbolic_protocol.json、symbolic_training_sample.csv、symbolic_pareto.csv和symbolic_candidates.json。固定seed与串行搜索提高复现性，但不同CPU/库版本可能产生不同的浮点搜索轨迹；冻结表达式的递推复算与从头重复搜索应区分。

若需新建独立依赖目录，可采用下列版本；不要把这些命令当作本轮已经执行了新的训练：

```powershell
py -3.12 -m pip install --target .deps --no-deps imageio-ffmpeg==0.6.0 pysr==2.7.2 juliacall==0.9.36 juliapkg==0.1.27 semver==3.1.0 tomli==2.5.0 tomlkit==0.15.2 sympy==1.14.0
```

Julia下载器使用现有Windows系统代理（不打印认证信息），包服务指定官方us-east服务器。首次安装的路径过长、包下载连接与Numba动态缓存失败日志均保留；这些是实现修复，不是修改水动力样本或隐藏科学失败。

## 数据与统计对应

- `provenance/input_manifest.json`：172份原始提取文件哈希、尺寸与修改时间；不包含测试结果。
- `tables/native_hdf_verification.csv`：172份原生HDF所有21个研究位置的完整时序逐值核对。
- `section_case_screen.csv`、`negative_episodes.csv`、`section_screen.csv`：严格负流量检查、负值时段和断面汇总。
- `paired_hydrodynamic_response.csv`、`response_summary.csv`：不使用模型精度的配对水动力选点依据。
- `provenance/event_split.json`：冻结49／21事件分组；`gate_amendment.json`记录网络训练前的门槛修订，原协议另存。
- `data/all.npz`与`metadata.json`：各工况自身的110—200 km输入、100 km目标、0 km外部边界及真实时序长度；填充尾部只用于矩阵存储，所有评价按ns截取。
- `predictions/*.npz`：每模型的全工况、全演算长度出流q、基础状态s、附加状态b、约束标记active；状态体积单位百万m³。
- `case_metrics.csv`、`event_metrics.csv`、`event_summary.csv`：逐工况逐河长、事件级与完整事件分布统计。
- `condition_summary.csv`、`length_summary.csv`、`phase_metrics.csv`、`phase_summary.csv`：两类边界、河长和各阶段比较。
- `paired_response_metrics.csv`、`paired_response_summary.csv`：70个主配对与48组控制因子对比，差值统一为明确标记的工况A减工况B；共有时间窗口写入每行。
- `control_background_event_metrics.csv`和`control_background_summary.csv`：先按背景洪水汇总控制条件，避免伪独立。逐控制工况分布另存。
- `timestep_sensitivity.csv`、`timestep_summary.csv`：固定参数与相同连续分段线性输入下的1、0.5、0.25 h比较。`0.25h_vs_0.5h`是真正相邻两细分输出之差，不代表另外计算了0.125 h。
- `balance.csv`：基础/附加状态和包括在途体积的物理控制体逐步残差。这是实现一致性，不等于HEC绝对储量验证。
- `P2_parameter_profile.csv`、`base_parameters.csv`、`neural_parameters.csv`、`model_complexity.csv`：条件剖面、边界命中、网络原始权重与最终运行复杂度。解析式的教师训练和结构搜索选择成本不能由少量运行常数消除。
- `animation/data_index.json`及`source_data_E007.npz`：实际HDF变量路径、全部101断面、逐小时物理时间与配对曲线。MP4、GIF、关键帧和make_animation.py一起交付。

所有主要表均提供平均、中位、最小、最大及原始逐例结果。总体模型比较先把同一洪水的河长和配对边界汇总，事件才是样本单位。对于基准零差值导致的KGE无定义，保留为空值而不伪造数值。未作跨河长伪独立显著性推断。

## 质量核查与本轮边界

verify.py核对原始文件不变、原生输出一致、完整递推守恒、预测方向、因果前缀、供水约束、实际PDF字体与字形尺寸、视频解码及GIF帧数。原有科研绘图工具的默认无衬线字体门槛被用户的Times New Roman／SimSun明确要求覆盖；其余问题必须修复。字体例外不被说成原默认门槛通过。

图件需检查完整页面和逐面板碰撞、差值方向及统计含义；自动源代码检查不能替代视觉检查。数值质量记录与科学适用性分开：本轮守恒和方向通过，跨步长参数可迁移性不标为通过；原封存测试保持未开启。阶段完成后等待研究者审阅，不自动进入独立测试或新一轮公式扩展。
