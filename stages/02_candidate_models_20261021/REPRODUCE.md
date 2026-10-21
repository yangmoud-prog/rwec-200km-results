# 复算与文件用途

本目录是独立的第一轮候选验证包。主要结论见REPORT.md，完整公式和限制见METHODS.md，冻结约定见PROTOCOL.md。现有HEC结果由程序以只读方式加载；30场测试目标不加载。

所需Python库为numpy、scipy、numba、pandas、matplotlib；使用本机已安装的Python3.12环境。字体要求Times New Roman与SimSun。绘图文字检查使用已安装科研绘图技能的检查器，该检查器不是率定依赖。

在当前目录运行以下步骤可使用已经保存的断点继续；fit中的既有主解与循环检查点会复用，不会主动删除任何结果。evaluate须在主率定完成时运行。如果要从零复算，应把代码和PROTOCOL.md复制到同一analysis目录下的另一新目录，保持原目录原状，并从prepare开始，而不是删除既有结果。

```powershell
py -3.12 campaign.py prepare
py -3.12 campaign.py preflight
py -3.12 campaign.py fit
py -3.12 evaluate.py
py -3.12 plot_results.py
py -3.12 verify_delivery.py
```

verify_delivery.py默认只完成可重复的数值、来源及导出检查，将视觉检查状态留为pending。逐面板检查实际图后可加--visual-review-confirmed；仅在需要同步项目入口时再加--update-route。本轮交付已完成33个面板的实际查看，不将自动检查冒充视觉检查。

campaign.py负责开发数据准备、与只读基准代码的一致性检查、固定预算主率定。models.py只包含输入、状态与参数到预测的演算，不接收目标出流。evaluate.py输出指标、控制差值、条件共享诊断与步长敏感性；plot_results.py生成完整图。主计算使用2个程序工作线程，不创建新的聊天或额外候选模型。

data/development.npz包含70场开发事件的20个入口流量、同一事件出口水位与出口目标。data/controls.npz为32项已有控制计算。metadata.json保存事件与边界时序。每个事件只使用ns指定的有效时段；矩形数组尾部填充值不参与目标或指标。

fits保存共同初始化、三模型各4轮改进、主解与条件诊断；每个主解记录目标值、评价次数、耗时及搜索轨迹。provenance保存源文件摘要、只读代码副本、论文提取、输入清单、冻结方案和求解检查。predictions包含开发与控制的完整流量和模型蓄量状态，以及全部条件差值和步长示例。模型状态不是HEC几何蓄量，不应将其误用于真实控制体蓄量检验。

tables中的event_length_metrics.csv是全部6120条模型—事件—河长主指标记录。metric_summary.csv与metric_by_length.csv给出主要指标的均值、中位数、最小、最大及有限数；phase_*为阶段明细。improvements.csv保存所有改善和退化。control_response_metrics.csv含48个条件对比×20河长×3模型=2880条差值指标。parameters.csv含原参数和主要共享参数，conditional_*仅为条件稳定性诊断，不能替换主要预测参数。saturation_product_profile.csv是固定乘积敏感性，不是另一候选模型。

本轮完成后按用户指定停止，不读取30场测试目标、不增加第三项修正、不扩大参数边界。复算若修改公式、输入、分组、参数边界或目标，应作为新的研究阶段另存；不能覆盖本轮结果后继续沿用原冻结标记。
