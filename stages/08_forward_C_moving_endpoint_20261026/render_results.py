from common import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont,FontProperties
import time

mpl.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':8,'axes.titlesize':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'svg.fonttype':'none','pdf.fonttype':42,'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'lines.linewidth':1.1,'axes.unicode_minus':False})
COL={'RWEC0':'#697a8c','NLM':'#006c91','LM':'#b2a690','MC_oracle':'#875aa3','MC_autonomous':'#d1493f'}
LAB={'RWEC0':'原 RWEC','NLM':'非线性马斯京根','LM':'线性马斯京根','MC_oracle':'新模型：理想末端','MC_autonomous':'新模型：自主末端'}
STYLE={'RWEC0':'--','NLM':'-.','LM':':','MC_oracle':'--','MC_autonomous':'-'}
qa=[]
def export(fig,name,claim):
 fig.canvas.draw()
 fig.savefig(P/'figures'/f'{name}.png',dpi=300)
 fig.savefig(P/'figures'/f'{name}.pdf')
 fig.savefig(P/'figures'/f'{name}.svg')
 qa.append(dict(figure=name,claim=claim,panels=len(fig.axes),final_width_mm=180,font_floor_pt=7,mean_or_raw='means explicitly labeled; event panels raw',source='tables and predictions',excluded_events=0))
 plt.close(fig)
def axes_grid(r,c,h):return plt.subplots(r,c,figsize=(180/25.4,h/25.4),squeeze=False,layout='constrained')
def legend(fig,ax):
 handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,loc='outside upper center',ncol=3,frameon=False)
def distance(ax):
 xx=np.arange(0,121,20);ax.set_xticks(xx,[f'{x}\n{200-x}' for x in xx]);ax.set_xlabel('距下游出口 / 距上游入口（km）')
def plotq(ax,eid,x):
 d=load(P/'data'/f'{eid}.npz');pr=load(P/'predictions'/f'{eid}_hydrographs.npz');j=np.argmin(abs(d['x']-x));node=int((200-x)/20)
 ax.plot(GRID/24,d['q'][:,j],color='#202020',lw=1.4,label='HEC-RAS')
 for name in COL:ax.plot(GRID/24,pr[name][node],color=COL[name],ls=STYLE[name],label=LAB[name])
 ax.set_title(f'洪水 {eid}：距出口 {x} km，距入口 {200-x} km');ax.set(xlabel='起涨后时间（d）',ylabel='流量（m³/s）');ax.grid(alpha=.15)

def main():
 plan=read(P/'provenance/C_frozen.json');fits=read(P/'provenance/frozen_parameters.json');v=plan['validation'];mean=pd.read_csv(P/'tables/fixed_section_mean_metrics.csv');ev=pd.read_csv(P/'tables/fixed_section_event_metrics.csv')
 save(P/'provenance/figure_contract.json',dict(backend='python',archetype='quantitative grid',claims=['physical screening uses HEC only','validation means and full hydrographs test accuracy','endpoint prediction is separated from oracle','storage and sweep test mechanism','raw improvements and dt changes expose failures'],n_validation=6,n_train=14,uncertainty='mean curves without bands as requested; all event differences separately shown; endpoint10-90% empirical neighbor envelope is not a calibrated interval',formats=['PNG300dpi','PDF','editableSVG'],font='SimSun Chinese; Times New Roman Latin/numerals/math',width_mm=180))
 # Hydraulic screening, all12 paired pilot scenarios retained.
 df=pd.read_csv(P/'tables/pilot_physical_screening.csv');fig,aa=axes_grid(1,2,80)
 for eid,g in df.groupby('event'):
  aa[0,0].plot(g.slope_per_mille,g.min_target_Q_m3s,'o-',label=eid);aa[0,1].plot(g.slope_per_mille,g.endpoint_excursion_km,'o-',label=eid)
 aa[0,0].axhline(0,color='black',lw=.7);aa[0,0].set(xlabel='河床坡度（‰）',ylabel='全目标断面最低流量（m³/s）',title='a  正向流动筛选')
 aa[0,1].set(xlabel='河床坡度（‰）',ylabel='末端最大移动范围（km）',title='b  配对水面差阈值 0.10 m');legend(fig,aa[0,0]);export(fig,'01_physical_screening','0.05per mille reverses in CP01; all pilot outcomes retained')
 fig,aa=axes_grid(3,1,175)
 for ax,met in zip(aa[:,0],['NSE','RMSE','KGE']):
  for name in COL:
   g=mean[(mean.subset=='validation')&(mean.model==name)].sort_values('x_from_outlet_km');ax.plot(g.x_from_outlet_km,g[met],marker='o',ms=3,color=COL[name],ls=STYLE[name],label=LAB[name])
  distance(ax);ax.set_ylabel('平均 '+met+('（m³/s）' if met=='RMSE' else ''));ax.grid(alpha=.15)
 aa[0,0].set_title('六场开发验证洪水：七个固定断面的等事件平均指标');legend(fig,aa[0,0]);export(fig,'02_validation_means','Mean validation accuracy at both distance coordinates')
 strength=[]
 for eid in v:
  d=load(P/'data'/f'{eid}.npz');j=np.argmin(abs(d['x']-20));strength.append((-np.min(d['q'][:,j]-d['qref'][:,j]),eid))
 ordered=[e for _,e in sorted(strength)];representatives=[ordered[-1],ordered[len(ordered)//2]]
 fig,aa=axes_grid(2,2,140)
 for row,eid in enumerate(representatives):
  for col,x in enumerate([20,0]):plotq(aa[row,col],eid,x)
 legend(fig,aa[0,0]);export(fig,'03_complete_hydrographs','Representatives selected by HEC suppression only; complete0-552h window')
 fig,aa=axes_grid(3,2,180)
 for ax,eid in zip(aa.flat,v):
  d=load(P/'data'/f'{eid}.npz');ax.fill_between(GRID/24,d['b_low'],d['b_high'],color='#e2e9ef',label='近邻 10—90% 范围')
  ax.plot(GRID/24,d['b_true'],color='black',label='HEC 真实末端');ax.plot(GRID/24,d['b_pred'],color=COL['MC_autonomous'],label='固定起点预测')
  ax.set(title=f'验证洪水 {eid}',xlabel='起涨后时间（d）',ylabel='距出口的末端位置（km）',ylim=(0,200))
 legend(fig,aa[0,0]);export(fig,'04_endpoint_positions','All6 validation endpoints; empirical range not confidence interval')
 fig,aa=axes_grid(3,2,180)
 from moving_balance import sample_flow
 for ax,eid in zip(aa.flat,v):
  d=load(P/'data'/f'{eid}.npz');pr=load(P/'predictions'/f'{eid}_hydrographs.npz');ax.plot(GRID/24,d['qb_true'],color='black',label='HEC 末端实际流量')
  for name in COL:ax.plot(GRID/24,sample_flow(pr[name],d['b_true']),color=COL[name],ls=STYLE[name],label=LAB[name])
  ax.set(title=f'验证洪水 {eid}：共同真实末端位置',xlabel='起涨后时间（d）',ylabel='断面流量（m³/s）')
 legend(fig,aa[0,0]);export(fig,'05_endpoint_physical_flow','Every model sampled at the same true moving position; physical not relative flux')
 eid=representatives[0];d=load(P/'data'/f'{eid}.npz');e=plan['events'][eid];times=[0,int(e['Hpeak_h']/2),int(e['Hpeak_h']),int(e['Hend_h'])]
 fig,aa=axes_grid(2,2,125)
 for ax,t in zip(aa.flat,times):
  ax.plot(d['x'],d['z'][t],color='#006c91',label='变动边界');ax.plot(d['x'],d['zref'][t],color='#a4866b',ls='--',label='固定 5 m 参照');ax.plot(d['x'],d['x']*1000*plan['selected_slope'],color='#555',lw=.6,label='河床')
  ax.axvline(d['b_true'][t],color='#d1493f',ls=':',lw=.8);xx=np.arange(0,201,40);ax.set_xticks(xx,[f'{x}\n{200-x}' for x in xx]);ax.set(title=f'洪水 {eid}，起涨后 {t/24:.1f} d',xlabel='距下游出口 / 距上游入口（km）',ylabel='水面高程（m）')
 legend(fig,aa[0,0]);export(fig,'06_longitudinal_profiles','Paired profiles and bed at4 actual times')
 fig,aa=axes_grid(2,2,140)
 for name in ['MC_oracle','MC_autonomous']:
  f=P/'predictions'/f'{eid}_{name}_states.npz'
  if not f.exists():continue
  s=load(f);ax=aa[0,0] if name=='MC_oracle' else aa[0,1];ax.plot(GRID/24,s['hec_storage_at_same_position']/1e6,color='black',label='HEC 同位置控制体');ax.plot(s['time']/24,s['sm']/1e6,color=COL[name],label=LAB[name]);ax.set(title=LAB[name],xlabel='起涨后时间（d）',ylabel='移动控制体储量（百万 m³）');ax.legend(frameon=False)
  ax=aa[1,0] if name=='MC_oracle' else aa[1,1];sec=float(s['time'][1]-s['time'][0])*3600
  ax.plot(s['time'][1:]/24,s['omitting_sweep_residual']/sec,color='#aaa',label='遗漏扫掠项');ax.plot(s['time'][1:]/24,s['balance_residual']/sec,color=COL[name],label='计入扫掠项');ax.set_yscale('symlog',linthresh=1.e-5);ax.set_yticks([-10000,-100,0,100,10000,100000],['−10000','−100','0','100','10000','100000']);ax.set(xlabel='起涨后时间（d）',ylabel='守恒残差（m³/s）');ax.legend(frameon=False)
 export(fig,'07_moving_storage_balance','Same-position geometry storage compared with own model states; algebraic closure distinguished from physical fit')
 raw=pd.read_csv(P/'tables/event_raw_improvements.csv');raw=raw[raw.subset=='validation'];fig,aa=axes_grid(1,2,85)
 for ax,name in zip(aa.flat,['MC_oracle','MC_autonomous']):
  for eid,g in raw.groupby('event'):ax.plot(g.x_from_outlet_km,g[f'{name}_minus_RWEC0'],'o-',ms=2.5,lw=.7,label=eid)
  ax.axhline(0,color='black',lw=.6);distance(ax);ax.set(title=LAB[name],ylabel='逐场 ΔNSE（新模型 − 原 RWEC）')
 legend(fig,aa[0,0]);export(fig,'08_raw_event_improvements','All event improvements and failures retained, no selected best cases')
 sens=pd.read_csv(P/'tables/model_dt_sensitivity.csv');fig,aa=axes_grid(1,2,85)
 for ax,x in zip(aa.flat,[20,0]):
  for name in COL:
   base=mean[(mean.subset=='validation')&(mean.model==name)&(mean.x_from_outlet_km==x)].NSE.iloc[0];yy=[]
   for factor in [.5,1.,2.]:
    gg=sens[(sens.model==name)&(sens.x_from_outlet_km==x)&(sens.dt_factor==factor)];yy.append(base if factor==1 else (gg.NSE.mean() if gg.valid.all() else np.nan))
   ax.plot([.5,1,2],yy,'o-',ms=3,color=COL[name],label=LAB[name])
  ax.set(title=f'距出口 {x} km，距入口 {200-x} km',xlabel='计算时段 / 已率定时段',ylabel='平均 NSE（参数固定）',xticks=[.5,1,2])
 legend(fig,aa[0,0]);export(fig,'09_timestep_sensitivity','Timestep changed with all physical and empirical parameters fixed')
 save(P/'provenance/figure_qa.json',dict(backend='python',figures=qa,font_paths=[findfont(FontProperties(family=f)) for f in ['Times New Roman','SimSun']],visual_inspection='pending',representative_selection='HEC20km suppression maximum and median among6validation; no model score used'))
 print('9 figures exported')

if __name__=='__main__':main()
