# 第九阶段：两区联合状态候选的方程与适用范围

本文件中的公式均为可编辑 LaTeX。以下为首轮待检验候选，不预设其优于原 RWEC。坐标 $x$ 从下游出口向上游增加，$L_0=200000\ \mathrm{m}$，$b$ 表示水位影响阈值边界，不是水体的物质边界。流量向下游为正。

## 连续水面形状与状态边界

梯形河槽满足

$$A(y)=B_0y+my^2,\quad T(y)=B_0+2my,\quad P(y)=B_0+2y\sqrt{1+m^2},\quad R=A/P.$$

本轮 $B_0=500\ \mathrm{m}$、$m=3$、$n_M=0.025$、$S_0=7.5\times10^{-5}$。正常水深 $y_n(q)$ 由曼宁式求得：

$$q=\frac{1}{n_M}A(y_n)R(y_n)^{2/3}\sqrt{S_0}.$$

以正常水深处的渐变流线性化长度为尺度，取

$$\lambda(q)=a_\lambda\frac{1-Fr_n^2}{-\left.\partial S_f/\partial y\right|_{y_n}},\qquad Fr_n^2=\frac{q^2T(y_n)}{gA(y_n)^3},$$

$$-\left.\frac{\partial S_f}{\partial y}\right|_{y_n}=2S_0\left[\frac{T}{A}+\frac{2}{3}\left(\frac{T}{A}-\frac{P'}{P}\right)\right].$$

这里 $\lambda$ 的单位为 m，$a_\lambda$ 无量纲。亚临界、缓变和近似均匀流是这一尺度的理论前提。实际洪水中的流量空间变化不能由该线性化自动恢复。

设预报起点已知水深为 $y_i(x)$，本步交界实际流量为 $q_b$、出口水位为 $h_0$。候选形状为

$$\widehat y(x)=y_i(x)+[y_n(q_b)-y_n(1000)]\left(1-e^{-x/\lambda(q_b)}\right)+(h_0-5)e^{-x/\lambda(q_b)}.$$

相同流量下的固定 5 m 参照候选由上式去掉最后一项得到。因此

$$\widehat{\Delta H}(x)=(h_0-5)e^{-x/\lambda(q_b)},$$

$$b(q_b,h_0)=\begin{cases}0,&h_0-5\leq\epsilon,\\\min[L_0,\lambda(q_b)\log((h_0-5)/\epsilon)],&h_0-5>\epsilon,\end{cases}\qquad\epsilon=0.10\ \mathrm{m}.$$

此参照是模型内部形状，不是未来 HEC 固定水位结果；配对 HEC 水面只用于评价。该候选仅表达正向水位抬升，负水位扰动不能用零末端解释为没有水力影响。

这里还有必须单列的结构近似：HEC参考是同一上游入流、固定5m边界的非恒定演算，而当前内部参照由同一本步交界流量构造，未独立递推参考洪水状态。两个参照并不等价。因此末端误差与指定末端处的阈值偏差，既包含水面形状误差，也包含参考状态定义的差异。本轮不将二者混同为某一参数的误差，更不以未来真实参考水面补齐自主模型。

下游区几何储量为

$$S_B(q_b,h_0)=\int_0^{b(q_b,h_0)}A(\widehat y(x))\,dx.$$

B 河道采用已知固定控制关系

$$Q_0=\mathcal C(h_0),\qquad h_0=3.5+1.5(Q_0/1000)^{2/3}.$$

计算使用原 HEC 工程中 241 对按 $50\ \mathrm{m^3/s}$ 间隔离散、按五位小数保存水位的线性插值表，避免把未来 HEC 水位作为边界输入。

## 上游区的互补储量递推

将全河分为十个 20 km 单元。单元 $j$ 位于 $[x_{j,l},x_{j,h}]$，末端上游的活动长度比例为

$$f_j(b)=\frac{[x_{j,h}-\max(x_{j,l},\min(x_{j,h},b))]_+}{20000}.$$

用 $U_j$ 表示活动部分的几何储量估计，用 $R_j$ 表示 RWEC 内部互补储量状态。二者通过初始状态偏移联系，而不直接混同：

$$U_j=V_{i,j}(b)+R_j-R_{0,j}(b),\qquad R_{0,j}=1000K_{20}f_j.$$

$K_{20}$ 在本式中以 s 计，是全长 20 km 单元的参考储量时间。$V_{i,j}$ 是已知初始水面在活动部分的体积，单位均为 $\mathrm{m^3}$。

为了在长度和时段改变时保留相同参考储量尺度，本候选将原互补指数重新参数化。令 $r_j=K_{20}f_j/\Delta t$，用下式确定 $n_j$：

$$\left(\frac{r_j}{r_j+1}\right)^{n_j}+\left(\frac{r_j}{r_j+k}\right)^{n_j}=1.$$

该处理保留互补函数，但并非原 RWEC 固定指数的原样复制。因此与原模型的差异不能全部归因于边界反馈；关闭反馈的同结构消融是必要对照。候选只有 $k,K_{20},a_\lambda$ 三个率定参数，不另加传播时滞参数。

末端移动后先重映射自身的上一步储量。记 $W_{j,t}$ 为从上游区转入下游区的有符号扫掠体积，正号对应末端向上游推进：

$$\widetilde U_{j,t-1}=U_{j,t-1}-W_{j,t},\qquad\widetilde R_{j,t-1}=\widetilde U_{j,t-1}-V_{i,j}(b_t)+R_{0,j}(b_t).$$

推进时按旧上游活动单元内的均匀面积密度移出体积；回撤时将旧下游水面形状在扫掠区间的积分转入上游。两侧几何是近似闭合，不能因重映射守恒就宣称两侧水面连续。末端穿过单元时出现的误差和求解失败均保留。

在重映射状态上执行原互补函数：

$$S_{W,j}=\widetilde R_{j,t-1}+\overline Q_{j-1,t}\Delta t,\qquad S_{E,j}=\widetilde R_{j,t-1}+k\overline Q_{j-1,t}\Delta t,$$

$$\left(\frac{R_{j,t}}{S_{W,j}}\right)^{n_j}+\left(\frac{R_{j,t}}{S_{E,j}}\right)^{n_j}=1,\qquad \overline Q_{j,t}=\frac{S_{W,j}-R_{j,t}}{\Delta t}.$$

上游各活动单元依次递推，最后一个活动单元得到 $q_b^{R}(b_t)$。不存在活动单元时令其等于上游入流。负储量、负出流或无效互补输入作为失败，不能通过裁剪改为成功。

## 同一步联合求解与流量含义

对于移动上游控制体，连续关系为

$$\frac{dS_U}{dt}=I-Q_b-A_b\dot b,$$

下游区为

$$\frac{dS_B}{dt}=Q_b-Q_0+A_b\dot b.$$

两式相加才得到全河固定控制体的连续方程。离散重映射用 $W_t=\sum_jW_{j,t}$ 近似扫掠体积。本步联合求解 $q_b,h_0$：

$$F_1=q_b-q_b^R[b(q_b,h_0)]=0,$$

$$F_2=\frac{S_B(q_b,h_0)-S_{B,t-1}-W_t}{\Delta t}-q_b+\mathcal C(h_0)=0.$$

这样，水位和末端改变上游活动长度、参考储量及重映射状态；上游出流又通过下游区的蓄泄方程影响水位和末端。上一步两区储量构成时间记忆，并未预先提供未来长度序列。

采用有限差分雅可比、带回溯的阻尼 Newton 法，最多 35 次迭代；两个流量量纲残差均小于 $10^{-5}\ \mathrm{m^3/s}$ 才接受本步。记录不收敛时间，不从未来 HEC 状态补初值。

实际断面流量为 $q_b$；移动界面相对通量为

$$q_b^{rel}=q_b+W_t/\Delta t.$$

二者不得互相替代。在固定空间断面 $x$，由模型自身的下游累计体积 $V_{0:x}$ 计算

$$\overline Q_x=\overline Q_0+\frac{V_{0:x,t}-V_{0:x,t-1}}{\Delta t}.$$

这是一种基于假定空间储量分布的流量重建，并不等于独立求解了完整动量方程。固定断面精度必须与 HEC 实际流量独立比较。体积平衡几乎为零是数值实现的必要条件，不是正确水面形状或洪水演进的充分证据。

## 理想长度、自治预测与消融

理想长度试验把 $b_t$ 替换为 HEC 参考末端，但 $q_b,h_0$ 仍由联合方程计算。此时还要检查形状在指定末端处是否满足 $0.10\ \mathrm m$ 阈值；已知真实末端并不保证形状闭合正确。

自治预测只读取上游流量、初始水面、固定河道和控制关系，末端由状态求得。新模型不调用八近邻预测。第八阶段八近邻串联在 B 训练集内重新建立邻居库，作为外部预处理对照；训练事件留一，验证事件只查训练库。

关闭反馈时，用相同参数的固定 200 km RWEC 计算上游通量；当前及历史回水状态不再改变该通量。下游蓄泄仍递推，上游储量按连续方程重映射，用于保持账面可比。此消融有意关闭上游储量对出流的本构反馈，不能作为完整物理模型。它与联合模型共享率定参数，衡量的是这一组参数下开启反馈的作用，不代表重新率定后单向模型的最优性能。

## C 河道的不可辨识性

C 的未来独立水位不是状态的已知固定函数。删除 $Q_0=\mathcal C(h_0)$ 后，水量方程不足以唯一确定 $Q_0,h_0,b$。给定起点前历史水位无法消除这一自由度。对同一入流和起点状态下两个合法未来边界产生的出流 $Q^{(1)},Q^{(2)}$，若仅作等概率两情景思想实验，任一无未来边界信息的确定性预测的最小期望均方误差满足

$$\mathrm{MSE}_{\min}=\frac{1}{4T}\sum_t[Q_t^{(1)}-Q_t^{(2)}]^2.$$

该下界不是实际边界概率已知的结论。两情景响应包络也不是置信区间。只有固定水力控制、已知未来边界规律或新的边界观测约束能进一步闭合此问题。

## 理论来源与候选性质

渐变流水面线与下游控制处理参照 [HEC-RAS 一维水面线理论](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/7.0/theoretical-basis-for-one-dimensional-and-two-dimensional-hydrodynamic-calculations/1d-steady-flow-water-surface-profiles) 和 [HEC-RAS 下游边界说明](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/latest/theoretical-basis-for-one-dimensional-and-two-dimensional-hydrodynamic-calculations/1d-unsteady-flow-hydrodynamics/implicit-finite-difference-scheme/downstream-boundary-conditions)。本文件的连续形状、两区近似和互补状态映射为本轮自行提出的候选，并非 HEC-RAS 官方模型或已经验证成立的推导结果。
