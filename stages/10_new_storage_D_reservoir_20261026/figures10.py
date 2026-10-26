from common10 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import MaxNLocator
matplotlib.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'axes.unicode_minus':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold'})
NAMES={'original':'原 RWEC','new':'新 SE RWEC','general':'广义 SE RWEC','NLM':'非线性马斯京根'}
COLORS={'original':'#4C566A','new':'#BC5A45','general':'#217F9A','NLM':'#8970A5'}
STYLES={'original':'--','new':'-.','general':'-','NLM':':'}
def export(fig,name):
 fig.savefig(P/'figures'/f'{name}.png',dpi=300,bbox_inches='tight')
 fig.savefig(P/'figures'/f'{name}.pdf',bbox_inches='tight')
 fig.savefig(P/'figures'/f'{name}.svg',bbox_inches='tight')
 plt.close(fig)

def formula():
 df=pd.read_csv(P/'tables/formula_event_metrics.csv');v=df.query("subset=='validation'");means=v.groupby(['domain','model','x_from_outlet_km'])[['NSE','RMSE','KGE']].mean()
 fig,axs=plt.subplots(2,3,figsize=(11.5,6.0),layout='constrained')
 for i,domain in enumerate(['B','C']):
  for j,metric in enumerate(['NSE','RMSE','KGE']):
   ax=axs[i,j]
   for name in NAMES:
    d=means.loc[(domain,name)];ax.plot(200-d.index,d[metric],label=NAMES[name],color=COLORS[name],ls=STYLES[name],marker='o',markersize=3)
   ax.set_xticks([80,100,120,140,160,180,200],['120/80','100/100','80/120','60/140','40/160','20/180','0/200'],rotation=35);ax.set_xlabel('距出口 / 距入口（km）');ax.set_ylabel(metric+('（m³/s）' if metric=='RMSE' else ''));ax.set_title(f'{chr(97+i*3+j)}  河道 {domain}：6 场验证洪水均值');ax.grid(alpha=.15)
 axs[0,0].legend(loc='lower left');export(fig,'01_existing_formula_comparison')
 imp=pd.read_csv(P/'tables/formula_event_improvements.csv').query("subset=='validation'");fig,axs=plt.subplots(1,2,figsize=(9.0,3.4),layout='constrained')
 for ax,domain in zip(axs,['B','C']):
  d=imp.query('domain==@domain');events=list(d.event.unique());yy=np.arange(len(events))
  for j,name in enumerate(['new','general','NLM']):
   g=d.groupby('event')[name+'_minus_original'].mean().reindex(events);ax.plot(g,yy+(j-1)*.16,'o',label=NAMES[name],color=COLORS[name],markersize=5)
  ax.axvline(0,c='.5',lw=.7);ax.set_yticks(yy,events);ax.set_title(f'河道 {domain}：每场七断面平均改善');ax.set_xlabel('相对原 RWEC 的 NSE 差值');ax.grid(axis='x',alpha=.15)
 axs[1].legend(loc='best');export(fig,'02_existing_paired_gains')
 for domain in ['B','C']:
  # Fixed a priori representative events; no selection based on validation gain.
  events=['E006','E051','E037'];fig,axs=plt.subplots(3,2,figsize=(10.8,8.2),layout='constrained')
  for i,e in enumerate(events):
   d=data(domain,e)
   for j,x in enumerate([20,0]):
    ax=axs[i,j];idx=int(np.argmin(abs(d['x']-x)));ax.plot(d['time']/24,d['q'][:,idx],c='#111111',lw=1.4,label='HEC-RAS')
    for name in NAMES:
     pp=load(P/'predictions'/f'{domain}_{e}_{name}.npz');node=int((200-x)/20);ax.plot(pp['time']/24,pp['q'][node],c=COLORS[name],ls=STYLES[name],lw=1.0,label=NAMES[name])
    ax.set_title(f'{e}：距出口 {x} km / 距入口 {200-x} km');ax.set_xlabel('预报起点后时间（天）');ax.set_ylabel('流量（m³/s）');ax.grid(alpha=.15)
  axs[0,0].legend(ncol=2);export(fig,f'03_{domain}_full_hydrographs')
 sens=pd.read_csv(P/'tables/formula_fixed_parameter_dt_sensitivity.csv');fig,axs=plt.subplots(1,2,figsize=(8.8,3.4),layout='constrained')
 for ax,domain in zip(axs,['B','C']):
  for name in NAMES:
   d=sens.query('domain==@domain and model==@name').groupby('dt_factor').NSE.mean();base=v.query('domain==@domain and model==@name').NSE.mean();ax.plot([.5,1,2],[d.loc[.5],base,d.loc[2]],color=COLORS[name],ls=STYLES[name],marker='o',label=NAMES[name])
  ax.set_xticks([.5,1,2]);ax.set_xlabel('固定参数的计算时段倍数');ax.set_ylabel('验证洪水平均 NSE');ax.set_title('河道 '+domain);ax.grid(alpha=.15)
 axs[1].legend();export(fig,'04_formula_time_step_sensitivity')

def pilots():
 rows=pd.read_csv(P/'tables/D_pilot_pairs.csv');v=rows.query("version=='v3'") if 'v3' in rows.version.values else rows.query("version=='v2'");fig,axs=plt.subplots(2,2,figsize=(10,6.4),layout='constrained')
 for ax,(metric,title) in zip(axs.flat,[('endpoint_max_km','回水末端距出口的最大距离'),('max_suppression_20_80_m3s','最大抑流：距出口 20—80 km / 距入口 180—120 km'),('max_recovery_20_80_m3s','最大恢复：距出口 20—80 km / 距入口 180—120 km'),('extra_storage_max_m3','河道附加储量最大值')]):
  for geom,color in [('uniform','#8193A7'),('complex','#217F9A')]:
   dd=v.query('geometry==@geom');ax.plot(dd.event,dd[metric]/(1e6 if metric=='extra_storage_max_m3' else 1),marker='o',color=color,label='均匀河槽' if geom=='uniform' else '复杂河槽')
  ax.set_title(title,fontsize=8);ax.grid(alpha=.15);ax.set_ylabel('百万 m³' if metric=='extra_storage_max_m3' else ('km' if metric=='endpoint_max_km' else 'm³/s'))
 axs[0,0].legend();export(fig,'05_D_physical_pilots')
 version=v.version.iloc[0];event='DP01';fig,axs=plt.subplots(2,2,figsize=(10.5,6.3),layout='constrained')
 for geom,color in [('uniform','#8193A7'),('complex','#217F9A')]:
  root=D/'models'/f'{event}_{geom}_reservoir_{version}';d=load(root/'extracted.npz');r=load(D/'models'/f'{event}_{geom}_fixed_v2'/'extracted.npz');pair=load(root/'paired.npz');t=d['time']/24;mask=t>=0;label='均匀河槽' if geom=='uniform' else '复杂河槽'
  axs[0,0].plot(t[mask],d['reservoir_H'][mask],c=color,label=label);axs[0,1].plot(t[mask],pair['b_true'][mask],c=color,label=label);axs[1,0].plot(t[mask],pair['extra_volume'][mask]/1e6,c=color,label=label)
  idx=int(np.argmin(abs(d['x']-20)));axs[1,1].plot(t[mask],d['q'][mask,idx],c=color,label=label+'：动态库区');axs[1,1].plot(t[mask],r['q'][mask,idx],c=color,ls='--',label=label+'：固定水位')
 for ax,title,unit in zip(axs.flat,['库水位由库容与泄量闭合','回水末端随时间移动','同一河道固定／动态边界的储量差','距出口 20 km / 距入口 180 km'],['库水位（m）','距出口（km）','附加储量（百万 m³）','流量（m³/s）']):ax.set_title(title);ax.set_ylabel(unit);ax.set_xlabel('起算后时间（天）');ax.grid(alpha=.15)
 axs[0,1].set_yticks([0,20,40,60,80,100],['0/200','20/180','40/160','60/140','80/120','100/100']);axs[0,1].set_ylabel('距出口 / 距入口（km）')
 axs[0,0].legend();axs[1,1].legend(fontsize=7);export(fig,'06_D_reservoir_backwater_process')

def closure():
 d=load(P/'predictions/oracle_DP01_complex_reservoir_v3_hydraulic_coordinate.npz');t=d['time']/24;x=d['x'];fig,axs=plt.subplots(2,2,figsize=(10.4,6.4),layout='constrained')
 im=axs[0,0].pcolormesh(t,x,d['delta_true'].T,shading='auto',cmap='Blues');fig.colorbar(im,ax=axs[0,0],label='水位差（m）');axs[0,0].plot(t,d['b'],c='#BC5A45',lw=1);axs[0,0].set_title('a  HEC 回水影响与真实末端');axs[0,0].set_ylabel('距出口（km）');axs[0,0].set_ylim(0,120)
 axs[0,0].set_yticks([0,20,40,60,80,100,120],['0/200','20/180','40/160','60/140','80/120','100/100','120/80']);axs[0,0].set_ylabel('距出口 / 距入口（km）')
 for h,style in [(48,'-'),(96,'--'),(168,':')]:
  i=int(np.argmin(abs(d['time']-h)));axs[0,1].plot(x,d['delta_true'][i],c='#111111',ls=style,label=f'HEC：{h} h');axs[0,1].plot(x,d['delta_shape'][i],c='#217F9A',ls=style,label=f'形状候选：{h} h')
 axs[0,1].set_xlim(0,120);axs[0,1].set_xticks([0,20,40,60,80,100,120],['0/200','20/180','40/160','60/140','80/120','100/100','120/80'],rotation=30);axs[0,1].set_xlabel('距出口 / 距入口（km）');axs[0,1].set_ylabel('水位差（m）');axs[0,1].set_title('b  相同非恒定参照下的水面形状');axs[0,1].legend(ncol=2,fontsize=7)
 axs[1,0].plot(t,d['extra_volume_true']/1e6,c='#111111',label='HEC 几何积分');axs[1,0].plot(t,d['extra_volume_shape']/1e6,c='#217F9A',ls='--',label='形状候选');axs[1,0].set_title('c  附加储量');axs[1,0].set_ylabel('附加储量（百万 m³）');axs[1,0].legend()
 axs[1,1].plot(t[1:],d['qb_shape']-d['qb_true'],c='#BC5A45',label='候选形状产生的通量误差');axs[1,1].plot(t[1:],d['HEC_balance'],c='#111111',label='HEC 小时离散残差');axs[1,1].set_title('d  储量变化与移动扫掠的独立检验');axs[1,1].set_ylabel('流量残差（m³/s）');axs[1,1].legend(fontsize=7)
 for ax in [axs[0,0],axs[1,0],axs[1,1]]:ax.set_xlabel('起算后时间（天）');ax.set_xlim(0,t.max())
 for ax in axs.flat:ax.grid(alpha=.12)
 export(fig,'07_oracle_shape_and_physical_flux')

if __name__=='__main__':
 if len(sys.argv)==1 or sys.argv[1]=='formula':formula()
 elif sys.argv[1]=='pilots':pilots()
 else:closure()
