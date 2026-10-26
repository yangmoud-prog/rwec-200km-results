from common import *
from joint import shape,length_scale
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont,FontProperties

plt.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':11,'axes.titlesize':11,'axes.labelsize':11,'legend.fontsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold','axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white'})
COL={'HEC':'#20252a','RWEC0':'#42739c','NLM':'#8a8c91','Serial8':'#c28b38','Joint':'#b94e58','Joint_oracle':'#397d76','Feedback_off':'#8b67a0'}
LABEL={'HEC':'HEC-RAS','RWEC0':'原 RWEC','NLM':'非线性马斯京根','Serial8':'第八阶段串联','Joint':'自主联合','Joint_oracle':'真实末端联合','Feedback_off':'关闭反馈'}
STYLE={'HEC':'-','RWEC0':'--','NLM':'-.','Serial8':':','Joint':'-','Joint_oracle':'--','Feedback_off':':'}
NAMES=['RWEC0','NLM','Serial8','Joint','Joint_oracle','Feedback_off']
exported=[]
def savefig(fig,name,caption):
 fig.savefig(P/'figures'/f'{name}.png',dpi=300,bbox_inches='tight')
 fig.savefig(P/'figures'/f'{name}.pdf',bbox_inches='tight')
 fig.savefig(P/'figures'/f'{name}.svg',bbox_inches='tight')
 exported.append(dict(file=name,caption=caption,panel_audit='pending',data_integrity='all six validation events; no validation-driven event exclusion'))
 plt.close(fig)
def line(ax,t,q,name):return ax.plot(t,q,color=COL[name],ls=STYLE[name],lw=1.7 if name in ['HEC','Joint'] else 1.2,label=LABEL[name])[0]
def decorate(ax,xlabel='预报起点后的时间（h）',ylabel='流量（m³/s）'):
 ax.set_xlabel(xlabel);ax.set_ylabel(ylabel);ax.grid(axis='y',color='#e3e7e9',lw=.55);ax.margins(x=.015)
def dual(ax):
 top=ax.secondary_xaxis('top',functions=(lambda v:200-v,lambda v:200-v));top.set_xlabel('距上游入口（km）',fontsize=10)
 if ax.get_xlim()[1]<150:ax.set_xticks(np.arange(0,121,20));top.set_xticks(np.arange(80,201,20))
 else:ax.set_xticks(np.arange(0,201,40));top.set_xticks(np.arange(0,201,40))
def legend(fig,ax,ncol=3):
 handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,1),ncol=ncol,frameon=False)

def main():
 cfg=read(P/'provenance/frozen_protocol.json');fits=read(P/'PARAMETERS.json');valid=cfg['validation'];table=pd.read_csv(P/'tables/fixed_section_means.csv');table=table[table.subset=='validation']
 assert np.all(np.diff(GRID)>0)
 for e in valid:assert np.all(np.diff(load(P/'data'/f'B_{e}.npz')['x'])>0)
 contract=dict(backend='python',archetype='quantitative grids',core_conclusion='The first coupled candidate does not outperform original RWEC/NLM and has endpoint and time-step stability failures; C retains boundary non-identifiability.',panels='mean skill; full hydrographs; endpoint location and same-location flux; storage/relative flux; water profiles; paired ablation and sensitivity; C two-scenario bounds',statistics='six fixed development-validation events, mean metrics without uncertainty bands; individual failures retained',font_files={f:findfont(FontProperties(family=f)) for f in ['Times New Roman','SimSun']},exports='PNG 300dpi, editable SVG/PDF; research review, not final journal submission',sealed_test_opened=0)
 save(P/'provenance/figure_contract.json',contract)
 fig,axs=plt.subplots(1,3,figsize=(13,4.7));fig.subplots_adjust(top=.75,bottom=.15,wspace=.34)
 for ax,metric,title in zip(axs,['NSE','RMSE','KGE'],['平均纳什效率','平均均方根误差','平均 KGE']):
  for name in NAMES:
   d=table[table.model==name].sort_values('x_from_outlet_km');line(ax,d.x_from_outlet_km,d[metric],name)
  ax.set_title(title,pad=42);decorate(ax,'距下游出口（km）',metric+('（m³/s）' if metric=='RMSE' else ''));dual(ax)
 legend(fig,axs[0],3);savefig(fig,'01_validation_mean_skill','6场验证洪水在7个固定断面的逐场指标算术平均。上下横轴分别表示距出口与距入口距离。均值不替代逐场失败记录。')
 fig,axs=plt.subplots(6,2,figsize=(12,17));fig.subplots_adjust(top=.95,bottom=.045,hspace=.57,wspace=.22)
 for row,e in enumerate(valid):
  d=load(P/'data'/f'B_{e}.npz');truth=.5*(d['truth'][:,1:]+d['truth'][:,:-1])
  for col,x in enumerate([20,0]):
   ax=axs[row,col];j=list(TARGETS).index(x);line(ax,GRID[1:]-.5,truth[j],'HEC')
   for name in ['RWEC0','NLM','Serial8','Joint']:
    p=load(P/'predictions'/f'B_{e}_{name}.npz');line(ax,GRID[1:]-.5,p['hourly_mean'][j+4],name)
   ax.set_title(f'验证洪水 {e}｜距出口 {x} km（距入口 {200-x} km）');decorate(ax);ax.set_xlim(0,288)
 legend(fig,axs[0,0],5);savefig(fig,'02_all_validation_hydrographs','全部6场验证洪水、距出口20km与0km的完整0—288h小时平均流量，无裁剪，无选择性省略。')
 fig,axs=plt.subplots(3,2,figsize=(12,10));fig.subplots_adjust(top=.92,hspace=.43,wspace=.25)
 for ax,e in zip(axs.flat,valid):
  d=load(P/'data'/f'B_{e}.npz');line(ax,GRID,d['b_true'],'HEC');line(ax,GRID,d['b_pred'],'Serial8')
  for name in ['Joint','Feedback_off']:
   p=load(P/'predictions'/f'B_{e}_{name}.npz');line(ax,p['time'],p['b'],name)
  ax.set_title(f'验证洪水 {e}｜动态末端位置');decorate(ax,ylabel='末端距下游出口（km）');ax.set_ylim(0,110)
 legend(fig,axs[0,0],4);savefig(fig,'03_endpoint_positions','末端坐标从出口向上游量取；对应入口至末端演算长度为200减去纵轴值。真实末端联合使用真实位置，故不作为位置预测曲线重复绘制。')
 fig,axs=plt.subplots(3,2,figsize=(12,10));fig.subplots_adjust(top=.92,hspace=.43,wspace=.25)
 for ax,e in zip(axs.flat,valid):
  d=load(P/'data'/f'B_{e}.npz');b=.5*(d['b_true'][1:]+d['b_true'][:-1]);q=.5*(d['q'][1:]+d['q'][:-1]);qt=np.array([np.interp(v,d['x'],row) for v,row in zip(b,q)]);line(ax,GRID[1:]-.5,qt,'HEC')
  for name in ['RWEC0','NLM','Serial8','Joint']:
   p=load(P/'predictions'/f'B_{e}_{name}.npz');pq=np.array([np.interp(v,np.arange(0,201,20),p['hourly_mean'][::-1,i]) for i,v in enumerate(b)]);line(ax,GRID[1:]-.5,pq,name)
  ax.set_title(f'验证洪水 {e}｜同一真实动态末端');decorate(ax)
 legend(fig,axs[0,0],5);savefig(fig,'04_same_dynamic_location_flux','所有模型在相同真实末端的小时中点位置插值评价。该位置只用于评价，不反馈自主预测；流量为实际断面流量，未加入扫掠相对通量。')
 e=valid[0];d=load(P/'data'/f'B_{e}.npz');p=load(P/'predictions'/f'B_{e}_Joint.npz');bp=np.interp(GRID,p['time'],p['b']);cum=cumulative_trapezoid(d['area'],d['x']*1000,axis=1,initial=0);hsb=np.array([np.interp(b,d['x'],v) for b,v in zip(bp,cum)]);hsu=cum[:,-1]-hsb
 fig,axs=plt.subplots(2,2,figsize=(12,9));fig.subplots_adjust(hspace=.4,wspace=.28,top=.92)
 for ax,key,hec,title in [(axs[0,0],'sb',hsb,'下游区储量：均按预测末端分区'),(axs[0,1],'su',hsu,'上游区储量：均按预测末端分区')]:
  line(ax,GRID,hec/1e8,'HEC');line(ax,p['time'],p[key]/1e8,'Joint');ax.set_title(title);decorate(ax,ylabel='储量（10⁸ m³）');ax.legend(frameon=False)
 ax=axs[1,0];line(ax,GRID,(cum[:,-1]-cum[0,-1])/1e8,'HEC');line(ax,p['time'],(p['sb']+p['su']-p['su'][0])/1e8,'Joint');ax.set_title('全河储量相对初始值的变化');decorate(ax,ylabel='储量变化（10⁸ m³）');ax.legend(frameon=False)
 ax=axs[1,1];ax.plot(p['time'],p['qb'],color=COL['Joint'],label='实际交界流量');ax.plot(p['time'],p['qb']+p['sweep_rate'],color=COL['Feedback_off'],ls='--',label='移动界面相对通量');ax.set_title('实际流量与相对通量分开');decorate(ax);ax.legend(frameon=False)
 fig.suptitle(f'验证洪水 {e}｜联合模型两区蓄泄与扫掠通量',y=.98);savefig(fig,'05_storage_and_relative_flux','预先固定展示首场验证洪水E006。前两图在同一预测边界位置积分HEC储量，避免把不同空间区间直接作差；图中接近不证明末端位置正确。')
 closure=load(P/'predictions'/f'closure_{e}.npz');fig,axs=plt.subplots(2,2,figsize=(12,9));fig.subplots_adjust(hspace=.65,wspace=.25,top=.86)
 for ax,t in zip(axs.flat,[24,48,96,168]):
  y=d['z'][t]-d['x']*.075;ax.plot(d['x'],y,color=COL['HEC'],label='HEC实际水深');ax.plot(d['x'],closure['z_shape'][t]-d['x']*.075,color=COL['Serial8'],ls='--',label='已知真实流量/水位的形状检验')
  index=np.argmin(abs(p['time']-t));bb,aa,cc=shape(p['qb'][index],p['h'][index],np.exp(fits['Joint']['theta'][2]),d['yinit']);xx=np.r_[d['x'][d['x']*1000<bb],bb/1000];yy=(-500+np.sqrt(500**2+12*aa))/6;ax.plot(xx,np.interp(xx,d['x'],yy),color=COL['Joint'],lw=2,label='自主联合的下游形状');ax.axvline(d['b_true'][t],color=COL['HEC'],ls=':',lw=.9)
  ax.set_title(f'{e}｜第 {t} h，虚线标真实末端',pad=42);decorate(ax,'距下游出口（km）','水深（m）');dual(ax);ax.set_xlim(0,200)
 legend(fig,axs[0,0],2);savefig(fig,'06_water_profile_closure','水位减河床高程得到水深。橙线使用真实Qb/H0作诊断，不是自主预报；红线仅在模型自身下游区绘制，未把未使用的形状延伸冒充上游RWEC状态。')
 paired=pd.read_csv(P/'tables/paired_NSE_changes.csv');paired=paired[paired.subset=='validation'];sens=pd.read_csv(P/'tables/fixed_parameter_dt_sensitivity.csv');fig,axs=plt.subplots(1,2,figsize=(12,5.8));fig.subplots_adjust(wspace=.28,bottom=.29,top=.83)
 ax=axs[0]
 for e in valid:
  v=paired[paired.event==e].sort_values('x_from_outlet_km');ax.plot(v.x_from_outlet_km,v.Joint_minus_Feedback_off,marker='o',ms=3,lw=.9,label=e)
 ax.axhline(0,color='#444444',lw=.8);ax.set_title('逐场开启反馈的 NSE 变化',pad=42);decorate(ax,'距下游出口（km）','联合减去关闭反馈');dual(ax);ax.legend(ncol=3,fontsize=8,frameon=False)
 ax=axs[1]
 for name in ['RWEC0','NLM','Joint','Feedback_off']:
  means=[]
  for factor in [.5,1,2]:
   v=table[table.model==name] if factor==1 else sens[(sens.model==name)&(sens.dt_factor==factor)]
   complete=bool(v.complete.all()) if factor!=1 else True
   means.append(v.NSE.mean() if complete else np.nan)
  ax.plot([.5,1,2],means,color=COL[name],ls=STYLE[name],marker='o',label=LABEL[name])
 fig.text(.72,.015,'半时段：联合 6/6 失败；关闭反馈 1/6 失败\n失败时不画幸存者均值',ha='center',fontsize=9);ax.set_xticks([.5,1,2]);decorate(ax,'相对于选定时段的倍数','平均 NSE');ax.set_title('固定参数的时段敏感性');ax.legend(fontsize=9,frameon=False,loc='upper center',bbox_to_anchor=(.5,-.23),ncol=2)
 savefig(fig,'07_feedback_ablation_and_dt','左图6场验证洪水逐场差值；右图仅在6场均完成时绘制均值，失败时留空并注明，未用幸存者均值代替总体结果。')
 c=load(P/'predictions/C_counterfactual.npz');fig,axs=plt.subplots(2,2,figsize=(12,9));fig.subplots_adjust(hspace=.4,wspace=.25,top=.93)
 for ax,x in zip(axs[0],[20,0]):
  j=int(x/2);ax.plot(c['time'],c['q1'][:,j],color=COL['RWEC0'],label='未来水位涨幅 1 m');ax.plot(c['time'],c['q4'][:,j],color=COL['Joint'],label='未来水位涨幅 4 m');ax.set_title(f'C 河道｜距出口 {x} km（距入口 {200-x} km）');decorate(ax);ax.legend(frameon=False)
 ax=axs[1,0];ax.plot(c['time'],c['b1'],color=COL['RWEC0'],label='涨幅 1 m');ax.plot(c['time'],c['b4'],color=COL['Joint'],label='涨幅 4 m');decorate(ax,ylabel='末端距下游出口（km）');ax.set_title('相同上游洪水与初始状态，不同动态末端');ax.legend(frameon=False)
 bounds=pd.read_csv(P/'tables/C_identifiability_two_scenario_bounds.csv').sort_values('x_from_outlet_km');ax=axs[1,1];ax.plot(bounds.x_from_outlet_km,bounds.equal_probability_RMSE_lower_bound_m3s,color=COL['Joint'],marker='o');ax.set_title('两情景等概率假设下的 RMSE 下界',pad=42);decorate(ax,'距下游出口（km）','不可消除的 RMSE（m³/s）');dual(ax)
 savefig(fig,'08_C_boundary_identifiability','沿用已完成的同入流、同初始状态反事实HEC配对。两种未来边界均物理允许；等概率仅用于计算信息不足的理论下界，不代表已知实际概率。')
 save(P/'provenance/figure_qa.json',dict(contract=contract,figures=exported,all_validation_events_shown=True,animation=False,render_review='pending'))
 print('FIGURES_READY',len(exported))
if __name__=='__main__':main()
