# 动态末端、蓄泄与RWEC的可编辑方程

以下数学表达均保留LaTeX源码。空间坐标从下游出口向上游增加，$x=0$为出口，$L_0=200000\ \mathrm{m}$为上游输入断面；流量$Q$向下游为正。

## 水位差末端与有效长度

动态工况与相同入流、固定5 m参照工况的水面差为

$$\Delta H(x,t)=H_{\mathrm{var}}(x,t)-H_{5}(x,t).$$

本轮保持原定义，主阈值为$\varepsilon=0.10\ \mathrm{m}$，同时检查0.05与0.20 m。

$$l_b(t;\varepsilon)=\sup\{x\in[0,L_0]:\Delta H(x,t)>\varepsilon\}.$$

空集时取0；最上游仍超过阈值时记录200 km右删失，不能认为真实影响恰好止于200 km。网格之间按水位差线性插值。另计算从出口开始始终满足阈值的连通区长度$l_{b,c}$。$l_b$是配对水位差影响位置，不必是物质界面，也不必代表连续的水力波前。

$$L_{\mathrm{free}}(t)=L_0-l_b(t),\qquad L_r=L_0-x_c,$$

$$L_{\mathrm{eff}}(x_c,t)=\max\{0,\min[l_b(t),L_0]-x_c\}.$$

$L_{\mathrm{free}}$、$L_r$与$L_{\mathrm{eff}}$的含义不同。$x_c=20,40,60,80,100,120\ \mathrm{km}$；体积计算统一换算为m。

## 几何附加蓄量与配对流量

梯形河槽断面积为

$$A(y)=B_0y+my^2,\qquad B_0=500\ \mathrm{m},\quad m=3.$$

$$\Delta V(x_c,t)=\int_{x_c}^{L_0}\{A[y_{\mathrm{var}}(x,t)]-A[y_5(x,t)]\}\,dx.$$

连续方程及固定区间水量平衡为

$$\frac{\partial A}{\partial t}-\frac{\partial Q}{\partial x}=0,\qquad
\frac{dS}{dt}=I-Q(x_c,t).$$

两工况上游输入相同，故

$$\boxed{\Delta Q(x_c,t)=-\frac{d\Delta V(x_c,t)}{dt}}.$$

因此蓄入阶段$\dot{\Delta V}>0$对应出流抑制，释放阶段对应出流恢复。这个关系约束的是几何蓄量与真实流量，不能未经证明把$\Delta V$直接认作RWEC潜在滞蓄$S_E$。本轮用面积积分梯形法和逐时中心差分检验，误差包含时间采样与空间插值误差。

仅长度与高度乘长度代理分别为

$$\widehat{\Delta V}_L=a_L L_{\mathrm{eff}},\qquad
\widehat{\Delta V}_{HL}=a_{HL}B_0[H_0(t)-5]_+L_{\mathrm{eff}}.$$

$a_L$具有面积量纲，$a_{HL}$无量纲。只用各自训练集拟合，验证不重估；这里使用真实高度、真实长度是机理诊断，不能列为自主预报结果。

## 移动末端的质量守恒

取上游未覆盖控制体$[b(t),L_0]$，$b(t)=l_b(t)$。

$$S_M(t)=\int_{b(t)}^{L_0}A(x,t)\,dx.$$

Leibniz积分求导给出

$$\boxed{\frac{dS_M}{dt}=I(t)-Q[b(t),t]-A[b(t),t]\dot b(t)}.$$

定义穿过移动界面的相对通量

$$Q_{\mathrm{rel}}(t)=Q[b(t),t]+A[b(t),t]\dot b(t),$$

则$\dot S_M=I-Q_{\mathrm{rel}}$。相对通量为负可能来自边界向下游退缩，不能当作固定断面的真实倒流。阈值末端可突然消失，不能无条件把$\dot b$当成平滑波速。本轮逐时采用扫掠体积

$$V_{\mathrm{sweep},i}=\int_{b_{i-1}}^{b_i}\frac{A(x,t_{i-1})+A(x,t_i)}{2}\,dx,$$

并用$V_{\mathrm{sweep},i}/\Delta t$代替点导数项作诊断。积分上限小于下限时按有符号积分处理，不裁剪。

## 固定起点的末端范围预测

预报起点$t_0=0$，只接受

$$\mathcal I_0=\{I(t):-48\le t\le240\ \mathrm{h};\ H_0(t):-48\le t\le0\ \mathrm{h};\ \text{固定河道条件}\}.$$

上游完整洪水过程是已知预报输入，不允许输入$t>0$的真实出口水位、内部水面线、流量、蓄量或末端。未来水位是否可由水力控制关系确定，仍取决于边界性质；A中的独立外部水位不可能被上游洪水唯一恢复。

特征使用固定12 h入流采样、峰值、固定窗口内峰时、洪水体积及固定6 h历史水位采样，采用事先固定的物理尺度。距离为

$$d_e=\frac1p\sum_{j=1}^p[f_j(\mathcal I_0)-f_{e,j}]^2.$$

从训练库选择距离最小的8场，权重为

$$w_e=\frac{(\sqrt{d_e}+0.02)^{-1}}{\sum_{r\in\mathcal N}(\sqrt{d_r}+0.02)^{-1}},\quad
\widehat l_b(t)=\sum_{e\in\mathcal N}w_e l_{b,e}(t).$$

经验范围为这8场训练轨迹逐时的未加权10%—90%分位数。该范围不是经过校准的80%置信带，必须报告实测覆盖率及宽度。训练场次采用按事件留一预测，验证场次仅查询训练库。A/B分别建库。预报前48 h的类比轨迹只用于在起点以前初始化候选模型，不声称这些历史轨迹是在其各自时刻发出的实时预报。

## 保留互补方程的最小动态长度候选

将200 km分为10个20 km子段。第$j$段上下游坐标为$x_u,x_d$，预测回水占比为

$$r_{j,i}=\frac{\max\{0,\min(\widehat l_{b,i},x_u)-x_d\}}{x_u-x_d},\qquad
\tau_{j,i}=\tau_{20}(1+\beta r_{j,i}),\quad\beta\ge0.$$

这是“回水增加传播滞时”的待检验假设，不是从完整动量方程推导出的必然规律。$\beta=0$回到原有固定传播滞时实现。

直接取$I[t-\tau(t)]\Delta t$会在变滞时时重复或遗漏入流水量。本轮使用单调出发时钟，以时间步索引表示为

$$s_{j,i}=\min\left\{i,\max\left[s_{j,i-1},i-\frac{\tau_{j,i}}{\Delta t}\right]\right\}.$$

入流采用历史线性插值，子段本步进入互补计算的水量为

$$V_{j,i}=\frac{I_j(s_{j,i})+I_j(s_{j,i-1})}{2}(s_{j,i}-s_{j,i-1})\Delta t.$$

其中流量用$\mathrm{m^3/s}$、$\Delta t$用s，时钟差为无量纲索引。设在途量$T$，实际进入子段的梯形积分水量为$V^{\mathrm{raw}}$：

$$T_{j,i}=T_{j,i-1}+V^{\mathrm{raw}}_{j,i}-V_{j,i}.$$

若$V_{j,i}>T_{j,i-1}+V^{\mathrm{raw}}_{j,i}$，按可用在途水量限制并记录次数；不对出流作零值裁剪。本轮未触发这个限制。

互补关系保持

$$S_{W,j,i}=S_{j,i-1}+V_{j,i},\qquad S_{E,j,i}=S_{j,i-1}+kV_{j,i},$$

$$\left(\frac{S_{j,i}}{S_{W,j,i}}\right)^n+\left(\frac{S_{j,i}}{S_{E,j,i}}\right)^n=1,$$

$$S_{j,i}=\left(S_{W,j,i}^{-n}+S_{E,j,i}^{-n}\right)^{-1/n},\qquad
Q_{j,i}=\frac{S_{W,j,i}-S_{j,i}}{\Delta t}.$$

常流初始化由$r$满足$(r/(r+1))^n+(r/(r+k))^n=1$确定，$S_0=rI_0\Delta t$、$T_0=I_0\tau_{j,0}$。$S$与$T$单位均为$\mathrm{m^3}$，$k,n,\beta$无量纲，$\tau_{20}$存储时用h，进入体积计算时换为s。

级联中一个子段的输出是步平均量，下段入流采用梯形积分。为准确审计离散实现，内部界面修正为

$$V_{J,i}=\frac{\Delta t}{2}\sum_{j=1}^{9}Q_{j,i}.$$

于是总离散平衡为

$$\Delta\left[\sum_j(S_j+T_j)+V_J\right]=\frac{I_{i-1}+I_i}{2}\Delta t-Q_{10,i}\Delta t.$$

该等式说明算法水量账闭合，不证明算法状态等于真实河槽几何储量。当前候选在固定子段上运行，没有把$Q_{\mathrm{rel}}$作为预测闭合关系，不能称为已经建立完整的移动双区物理模型。

## 率定、对照与不可辨识下界

各模型联合调整全部参数。A用49场训练，B用14场训练；每场六个断面等权，真实非负时刻用于拟合和主评分，完整过程不截断、不重启。训练目标为

$$J=\frac1{6N}\sum_{e=1}^N\sum_{c=1}^6
\frac{\sum_{i\in\mathcal T^+_{e,c}}(\widehat Q_{e,c,i}-Q_{e,c,i})^2}
{\sum_{i\in\mathcal T^+_{e,c}}(Q_{e,c,i}-\overline Q^+_{e,c})^2}.$$

RWEC-0率定$k,n,\tau_{20}$；预测长度与真实长度oracle分别联合率定$k,n,\tau_{20},\beta$。Oracle仅诊断信息与结构，不是可运行预报。线性马斯京根另率定$K_{20},X$，不等同于此前称为RWEC-1的线性水位修正模型。

如果两个场景具有同一允许输入，却有不同未来真实流量$q_a,q_b$，任一共同确定性预报$p$满足

$$\frac{\|p-q_a\|_2^2+\|p-q_b\|_2^2}{2}\ge\frac{\|q_a-q_b\|_2^2}{4}.$$

故最优共同预报在双场景平均意义上的RMSE下界为$\operatorname{RMSE}(q_a,q_b)/2$。这是特定配对场景的不能同时准确预测界，不是对单一事件NSE的普适上限。本轮同时保存历史水面差，区分完全相同与数值容差内相同。

## 物理依据

连续方程与下游边界类型参考[HEC-RAS一维非恒定流连续方程](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.1/theoretical-basis-for-one-dimensional-and-two-dimensional-hydrodynamic-calculations/1d-unsteady-flow-hydrodynamics/implicit-finite-difference-scheme/continuity-equation0)及[HEC-RAS下游边界说明](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/latest/theoretical-basis-for-one-dimensional-and-two-dimensional-hydrodynamic-calculations/1d-unsteady-flow-hydrodynamics/implicit-finite-difference-scheme/downstream-boundary-conditions)。本轮类比预测与动态延迟候选是本地实现和检验的研究假设，不能援引手册作为其有效性的证明。
