# 移动回水控制体与互补滞蓄候选：方程、单位及适用范围

本文件中的方程均为可编辑 LaTeX。它区分连续方程的严格结果、RWEC 原有经验关系和本轮需要检验的新增几何假设。严格守恒并不自动证明几何假设正确，也不能消除未知外部水位强迫造成的信息不足。

## 坐标与移动控制体

取距下游出口的距离 $x$ 为坐标，向上游为正，河道入口 $L_0=200000\,\mathrm m$，物理流量 $Q>0$ 表示向下游流动。固定断面连续方程为

$$\frac{\partial A}{\partial t}=\frac{\partial Q}{\partial x}.$$

以配对水面差 $\Delta H=H_{\mathrm{var}}-H_{\mathrm{ref}}\ge0.10\,\mathrm m$ 的最上游位置定义阈值末端 $b(t)$，同时保存出口连通末端、脱离出口的残留影响区、多个影响区和到达入口的删失。阈值边界不是物质边界，突跳不能直接解释为水体运动速度。上游移动控制体的长度和储量为

$$L(t)=L_0-b(t),\qquad S_M(t)=\int_{b(t)}^{L_0} A(x,t)\,dx.$$

应用变上限积分求导可得

$$\frac{dS_M}{dt}=\int_b^{L_0}\frac{\partial A}{\partial t}\,dx-A(b,t)\dot b
=I(t)-Q(b,t)-A(b,t)\dot b.$$

因此实际断面流量和移动界面相对通量分别是

$$\boxed{Q_b=I-\dot S_M-A_b\dot b},\qquad F_{\mathrm{rel}}=Q_b+A_b\dot b.$$

$S_M$ 的单位为 $\mathrm m^3$，$A_b\dot b$ 为 $\mathrm m^3\mathrm s^{-1}$。末端上移时 $\dot b>0$，控制体缩短；扫掠项不能再作为额外向下游的真实流量重复相加。

## 原 RWEC 与统一通量约定

每个20 km子段保留原始互补式，$R$ 表示递推储量变量，$n$ 和 $k$ 无量纲。滞后入流体积记为 $V_\tau$，计算时长 $\Delta t$ 以秒计：

$$S_W=R_{j,t-1}+V_{\tau,j,t},\qquad S_E=R_{j,t-1}+kV_{\tau,j,t},$$
$$\left(\frac{R_{j,t}}{S_W}\right)^n+\left(\frac{R_{j,t}}{S_E}\right)^n=1,
\qquad R_{j,t}=\left(S_W^{-n}+S_E^{-n}\right)^{-1/n},$$
$$Q_{j,t}=\frac{S_W-R_{j,t}}{\Delta t}.$$

滞后输入采用同一输入过程的梯形积分。为避免将滞后处理遗漏为不守恒的水量损失，显式保存输移储量 $T$：

$$T_{j,t}=T_{j,t-1}+\bar Q_{j-1,t}\Delta t-V_{\tau,j,t},
\qquad T_{j,0}=Q_0\tau.$$

初始 $R_0=rQ_0\Delta t$，$r$ 由下式求解，使恒定流量严格保持恒定：

$$\left(\frac r{r+1}\right)^n+\left(\frac r{r+k}\right)^n=1.$$

原 RWEC 出流是本步体积除以步长。为了在子段间和 HEC 时刻输出之间统一使用梯形通量，增加明确的记账转换项 $J_{j,t}=\tfrac12\Delta t\,Q_{j,t}$；它不是独立的物理水库。初始化几何体积为 $V_{j,0}$，模型总储量写为

$$\widehat V_{j,t}=V_{j,0}+(R_{j,t}-R_{j,0})+(T_{j,t}-T_{j,0})+(J_{j,t}-J_{j,0}).$$

于是 $\widehat V_{j,t}-\widehat V_{j,t-1}=\Delta t(\bar Q_{j-1,t}-\bar Q_{j,t})$。这说明递推变量、滞后储量和几何总储量不能不加说明地直接等同。原模型代码等价性另有数值核验。

## 本轮新增候选：空间分布暂存与释放

候选名称为 MC-RWEC。它在固定20 km网格保留内部状态，再对移动控制体进行部分单元积分；不在末端移动时清空或凭空重建状态。新增假设只有一个简化水深形状和一个响应时标：

$$\widehat h_b(x,t)=\alpha S_0[\widehat b(t)-x]_+,$$
$$B_j^{\mathrm{eq}}(t)=\int_{x_{j,d}}^{x_{j,u}}
\{[B_0+2my_{0,j}]\widehat h_b+m\widehat h_b^2\}\,dx,$$
$$B_{j,t}=B_{j,t-1}+\left(1-e^{-\Delta t/T_B}\right)
[B_j^{\mathrm{eq}}(t)-B_{j,t-1}].$$

$B_0=500\,\mathrm m$，边坡系数 $m=3$，$S_0=0.000075$，$y_{0,j}$ 仅来自起点已知稳定水深，$\alpha$ 无量纲，$T_B$ 为小时并在计算指数时与步长统一单位。$T_B=0$ 按瞬时响应处理。该形状并非渐变流方程的严格解；以阈值长度代替几何支撑范围、以初始水深代替随洪水变化的参照水深，均是本轮可被数据否定的近似。

暂存填充和释放通过水量收支进入原互补求解，不直接叠加一个必为正的潜在滞蓄项：

$$V_{j,t}^{\mathrm{net}}=V_{\tau,j,t}-(B_{j,t}-B_{j,t-1}),$$
$$S_W^*=R_{j,t-1}+V_{j,t}^{\mathrm{net}},\qquad
S_E^*=R_{j,t-1}+kV_{j,t}^{\mathrm{net}},$$
$$R_{j,t}=\left[(S_W^*)^{-n}+(S_E^*)^{-n}\right]^{-1/n},
\qquad Q_{j,t}=\frac{S_W^*-R_{j,t}}{\Delta t}.$$

增加的储量是 $B$ 本身，不能再同时把同一几何体积加进 $S_E$：

$$\widehat V_{j,t}=V_{j,0}+\Delta R_j+\Delta T_j+\Delta B_j+\Delta J_j,$$
$$\boxed{\widehat V_{j,t}-\widehat V_{j,t-1}
=\Delta t(\bar Q_{j-1,t}-\bar Q_{j,t})}.$$

若互补项非正、内部储量非正、计算出流为负或求解无效，则记录失败，不把流量截断为零，也不自动退回原模型。$\alpha=0$ 时严格退化为原 RWEC，因而另以固定候选参数、仅令 $\alpha=0$ 的消融试验检查新增暂存机制的独立作用。

## 末端移动时的状态更新和守恒核验

令 $f_j(b)$ 是第 $j$ 个20 km单元位于末端上游的长度比例，则

$$f_j(b)=\operatorname{clip}\left(\frac{x_{j,u}-b}{20000},0,1\right),
\qquad \widehat S_M(t)=\sum_j f_j[b(t)]\widehat V_j(t).$$

对于一个时间步，乘积差分恒等式为

$$\Delta\widehat S_M=\sum_j\bar f_j\Delta\widehat V_j+\sum_j\bar V_j\Delta f_j.$$

定义扫掠体积 $W_{\mathrm{sweep}}=-\sum_j\bar V_j\Delta f_j$。它等价于按单元内部面积估计对边界扫过的空间区间积分，可以处理跨越多个单元的末端突跳，不必把 $\dot b$ 数值尖峰当成水流速度。对应空间、时间双线性求积的断面流量满足

$$\boxed{\overline Q_b^{\mathrm{quad}}=\bar I-
\frac{\Delta\widehat S_M+W_{\mathrm{sweep}}}{\Delta t}}.$$

预报流量图采用固定网格流量在 $b(t)$ 处的空间插值 $\widehat Q[b(t),t]$。这是物理断面流量；$\overline Q_b^{\mathrm{quad}}$ 是守恒核验中的步平均求积量，两者不能偷换。报告同时提供采用时刻流量梯形平均后的离散残差。所有模型的移动位置评价采用相同位置，绝不只给新模型使用真实末端。

这一实现属于具有移动控制体守恒记账的分布暂存候选，并未求解完整动量方程。守恒恒等式本身不能构成优于传统模型的证据；必须看固定断面、共同移动位置、储量误差和消融结果是否一致支持改善。

## 非线性与线性马斯京根

非线性模型使用

$$S=K[XI+(1-X)Q]^p,\qquad
K=\frac{3600K_r}{pQ_r^{p-1}},\quad Q_r=1000\,\mathrm{m^3s^{-1}}.$$

$K_r$ 是在参考流量下的储量导数时间尺度（h），$X$ 为入流权重，$p$ 为独立率定的非线性指数。直接率定 $(K_r,X,p)$ 等价于给定参考流量后的 $(K,X,p)$ 一一变换，不固定指数。$K$ 的单位为 $\mathrm m^3/(\mathrm{m^3s^{-1}})^p$。$p=1$ 时为线性马斯京根。

隐式梯形连续方程为

$$S_t-S_{t-1}=\frac{\Delta t}{2}(I_{t-1}+I_t-Q_{t-1}-Q_t).$$

出流采用有界 Newton—二分法求根；无非负根或残差不满足要求时保留失败。线性极限已与标准显式系数公式逐点核验。

## 预测信息与不可辨识性

理想试验把 HEC 真实阈值末端作为给定条件，仅用于检验结构。自主试验只查询上游洪水预报与起点前48 h水位历史，沿用第七阶段的8近邻末端预测器。训练事件采用事件留一，开发验证只查询训练库；未来 HEC 水位、出流、储量和末端不进入自主模型。预测器的近邻范围不是经概率校准的置信区间。

相同入流与初始历史、不同未来外部水位能够产生 $Q_a\ne Q_b$。任何只有相同输入的确定性预测器都必须给出同一结果 $\widehat Q$。在两个情景等权时，

$$\frac{\|Q_a-\widehat Q\|^2+\|Q_b-\widehat Q\|^2}{2N}
\ge\frac{\|Q_a-Q_b\|^2}{4N}.$$

右端由独立水位反事实配对计算给出。增加移动末端机制能否改善条件平均过程，与是否能够唯一恢复每条未知未来边界，是两个不同问题。

## 技术来源

初始水面由稳定回水分析建立，另加24 h热启动与48 h恒定前史，并核验实际稳定程度；方法依据 [HEC-RAS 初始条件说明](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/latest/basic-data-requirements/unsteady-flow-data/initial-conditions)。线性马斯京根的系数及权重范围可与 [HEC Muskingum Routing](https://www.hec.usace.army.mil/confluence/cwmsdocs/rsum/muskingum-routing-58631546.html) 核对。移动控制体部分由本文所列连续方程直接推导，新增几何闭合属于本轮研究假设。
