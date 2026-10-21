# 诊断量与数据使用边界

令 $x$ 从共同出口向上游增加，$L$ 为嵌套河段长度，流量 $Q$ 仍以下泄方向为正。变动水位与固定 5 m 参照分别记为 $v$ 和 $r$。本模型无侧向流入，各个短河段与完整 200 km 河段共享同一出口；每个短河段的入口位于自身的 $x=L$ 断面，不能用完整河道的规定入口流量替代。

\[
\Delta Z(x,t)=Z^v(x,t)-Z^r(x,t),\qquad
\Delta V_L(t)=V_L^v(t)-V_L^r(t).
\]

固定 5 m 是增量比较参照。使用时变水位过程作为下游边界，可表示河道汇入已知水位变化的湖库或潮汐环境，这是 [HEC-RAS 下游边界说明](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/latest/theoretical-basis-for-one-dimensional-and-two-dimensional-hydrodynamic-calculations/1d-unsteady-flow-hydrodynamics/implicit-finite-difference-scheme/downstream-boundary-conditions) 所描述的边界类型；本工程没有另行构建水库蓄水区。

梯形断面的底宽为 $B=500\,\mathrm{m}$，边坡系数为 $m=3$，水深为 $y$。原生断面计算结果满足以下几何关系；所有沿程积分都在原有 2 km 计算断面上用梯形积分完成，没有从 10 km 控制断面重新近似整个水面。

\[
A(y)=By+my^2,
\qquad
\Delta V_L(t)=\int_0^L
\left[(B+2my^r)\Delta Z+m(\Delta Z)^2\right],\mathrm dx.
\]

质量守恒要求蓄量变化等于净流入，依据见 [HEC-RAS 连续方程](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.5/theoretical-basis-for-one-dimensional-and-two-dimensional-hydrodynamic-calculations/1d-unsteady-flow-hydrodynamics/continuity-equation)。对同一河长的两个工况作差，有

\[
\Delta V_L(t)-\Delta V_L(t_0)
=\int_{t_0}^t\left[\Delta Q(L,u)-\Delta Q(0,u)\right],\mathrm du+\varepsilon_L(t).
\]

采用 HEC 原生累计通量计算该积分。小时流量梯形积分仅作为输出取样误差的辅助诊断，不替代原生通量守恒检查。200 km 入口的规定流量相同，短河段的入口流量一般不同；式中入口项不能删除。原始负流量按下泄方向的符号保留。

每增加 $\delta L=10\,\mathrm{km}$ 的同步附加蓄量定义为

\[
W_j(t)=\Delta V_{L_j}(t)-\Delta V_{L_{j-1}}(t),
\qquad \sum_jW_j(t)=\Delta V_{200}(t).
\]

分别保存空间段自身的最大值与相邻河长峰值的差：

\[
M_j=\max_t W_j(t),\qquad
\Delta P_j=\max_t\Delta V_{L_j}(t)-\max_t\Delta V_{L_{j-1}}(t).
\]

二者不是同一统计量。图中的 $L_{90}$、$L_{95}$ 取 10 km 网格上首次达到各自 200 km 峰值的 90%、95% 的河长，不做空间内插，不能视为相对于无限长河道的饱和长度。

影响范围以 $\Delta Z>\eta$ 分类，$\eta=0.05,0.10,0.20\,\mathrm m$。最远正增量断面位置定义为

\[
x_\eta(t)=\max\{x>0:\Delta Z(x,t)>\eta\}.
\]

无内部超阈值断面时记为 0。规定水位的出口 $x=0$ 不参与这个判据；另存从第一个内部断面连续超阈值的长度、区段数、负超阈值断面数以及 $x=200\,\mathrm{km}$ 是否仍超阈值。达到上游边界的影响范围作右删失记录，不从现有计算域外推终点。完整有符号水位响应没有按阈值删除。

控制实验的 A、B、C 分别取两个预定水平。以编码 $s_{iA},s_{iB},s_{iC}\in\{-1,+1\}$ 表示八个组合，某个因素或交互项的即时效应为高水平平均减低水平平均。例如

\[
E_A(L,t)=\frac14\sum_{i=1}^8s_{iA}\Delta V_{L,i}(t),\qquad
E_{AB}(L,t)=\frac14\sum_{i=1}^8s_{iA}s_{iB}\Delta V_{L,i}(t).
\]

另存每一固定其他因素条件下的对比，以及各工况先取自身峰值再作差的指标。对每个背景的四个对比平均后形成即时效应曲线；即时平均效应的最大值与各条件峰值差的平均值不可互换。时间积分单位为 m³·h，明确作为蓄量过程指标使用。

分组来自冻结清单，未按结果修改。开发数据为 70 场率定配对和 32 次预定控制计算；30 场验证洪水的目标序列不进入本阶段程序。完整性审计覆盖全部 232 次的质量元数据。统计图中的 10—90% 阴影为率定事件分布，不是置信区间；不为确定性控制计算赋予独立重复数或推断性显著性。插值只用于已注明的边界恢复时刻及等边界水位对比，先验证时间坐标严格递增。

方法输出依赖现有 2 km 网格、60 s 计算步长和 1 h 输出；前期代表工况的网格与时间比较已通过冻结的工作容差。它们是理想化单一几何下的机制证据，不能替代真实河道的独立检验。
