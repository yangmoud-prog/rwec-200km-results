from common10 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':9,'axes.titlesize':10,'axes.labelsize':9,'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'axes.unicode_minus':False,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold'})
NAMES={'original':'原 RWEC','new':'新 SE RWEC','general':'广义 SE RWEC','NLM':'非线性马斯京根'}
COLORS={'original':'#4C566A','new':'#BC5A45','general':'#217F9A','NLM':'#8970A5'}
STYLES={'original':'--','new':'-.','general':'-','NLM':':'}
def export(fig,name):
 fig.savefig(P/'figures'/f'{name}.png',dpi=300,bbox_inches='tight')
 fig.savefig(P/'figures'/f'{name}.pdf',bbox_inches='tight')
 fig.savefig(P/'figures'/f'{name}.svg',bbox_inches='tight')
 plt.close(fig)
def spatial(ax):
 ax.set_xticks([80,100,120,140,160,180,200],['120/80','100/100','80/120','60/140','40/160','20/180','0/200'],rotation=30);ax.set_xlabel('距出口 / 距入口（km）');ax.grid(alpha=.15)
def main():
 df=pd.read_csv(P/'tables/validation_dynamic_event_metrics.csv');means=df.groupby(['window','model','x_from_outlet_km'])[['NSE','RMSE','KGE']].mean();fig,axs=plt.subplots(2,3,figsize=(11.5,6.0),layout='constrained')
 for i,(window,label) in enumerate([('full','完整洪水'),('nonnegative_hour_endpoints','真实端点均非负的共同小时')]):
  for j,metric in enumerate(['NSE','RMSE','KGE']):
   ax=axs[i,j]
   for name in MODELS:
    d=means.loc[(window,name)];ax.plot(200-d.index,d[metric],c=COLORS[name],ls=STYLES[name],marker='o',markersize=3,label=NAMES[name])
   spatial(ax);ax.set_ylabel(metric+('（m³/s）' if metric=='RMSE' else ''));ax.set_title(f'{chr(97+i*3+j)}  21 场验证：{label}',fontsize=9)
 axs[0,0].legend(fontsize=7);export(fig,'01_A_mean_skill_full_and_nonnegative')
 gains=pd.read_csv(P/'tables/event_mean_gains.csv').query("subset=='validation' and scenario=='dynamic'");fig,axs=plt.subplots(1,2,figsize=(10.3,6.4),layout='constrained')
 for ax,(window,label) in zip(axs,[('full','完整洪水'),('nonnegative_hour_endpoints','排除真实负流量后的共同小时')]):
  d=gains[gains.window==window].sort_values('event');yy=np.arange(len(d))
  for j,m in enumerate(['new','general','NLM']):ax.plot(d[m+'_minus_original'],yy+(j-1)*.18,'o',ms=3.5,c=COLORS[m],label=NAMES[m])
  ax.axvline(0,c='.5',lw=.7);ax.set_yticks(yy,d.event);ax.set_title(label);ax.set_xlabel('每场七断面平均 NSE 改善');ax.grid(axis='x',alpha=.15)
 axs[0].legend(fontsize=7);export(fig,'02_A_all21_paired_gains')
 cfg=read(P/'provenance/frozen_protocol.json')
 for group,events in [('representative',['E006','E037','E094']),('remaining_positive_examples',['E010','E051','E068'])]:
  fig,axs=plt.subplots(3,3,figsize=(11.5,8.3),layout='constrained')
  for i,event in enumerate(events):
   d=load(P/'data'/f'{event}.npz');keep=d['time']>=0
   for j,x in enumerate([100,20,0]):
    ax=axs[i,j];col=int(np.argmin(abs(d['x']-x)));q=d['q'][keep,col];ax.plot(d['time'][keep]/24,q,c='#111111',lw=1.4,label='HEC-RAS')
    for name in MODELS:
     pp=load(P/'predictions'/f'{event}_{name}.npz');mask=pp['time']>=0;ax.plot(pp['time'][mask]/24,pp['q'][int((200-x)/20),mask],c=COLORS[name],ls=STYLES[name],lw=1.,label=NAMES[name])
    if q.min()<0:ax.axhline(0,c='.6',lw=.6)
    ax.set_title(f'{event}：距出口 {x} km / 距入口 {200-x} km',fontsize=8.5);ax.set_xlabel('预报起点后时间（天）');ax.set_ylabel('流量（m³/s）');ax.grid(alpha=.15)
  axs[0,0].legend(fontsize=6.5);export(fig,'03_A_full_hydrographs_'+group)
 sens=pd.read_csv(P/'tables/fixed_parameter_dt_sensitivity.csv');fig,axs=plt.subplots(1,2,figsize=(10.3,3.8),layout='constrained')
 for ax,(window,label) in zip(axs,[('full','完整洪水'),('nonnegative_hour_endpoints','排除真实负流量后的共同小时')]):
  for m in MODELS:
   d=sens[(sens.model==m)&(sens.window==window)].groupby('dt_factor').NSE.mean();base=df[(df.model==m)&(df.window==window)].NSE.mean();ax.plot([.5,1,2],[d.loc[.5],base,d.loc[2]],c=COLORS[m],ls=STYLES[m],marker='o',label=NAMES[m])
  ax.set_xticks([.5,1,2]);ax.set_xlabel('各子段固定参数计算时段倍数');ax.set_ylabel('21 场七断面平均 NSE');ax.set_title(label);ax.grid(alpha=.15)
 axs[0].legend(fontsize=7);export(fig,'04_A_fixed_parameter_step_sensitivity')
 event='E037';d=load(P/'data'/f'{event}.npz');keep=d['time']>=0;from scipy.integrate import cumulative_trapezoid
 geom=cumulative_trapezoid(d['area'][keep],d['x']*1000,axis=1,initial=0);fig,axs=plt.subplots(1,3,figsize=(11.5,3.7),layout='constrained');true=geom[:,-1]-geom[:,2];true-=true[0];t=d['time'][keep]/24;axs[0].plot(t,true/1e6,c='#111111',lw=1.4,label='HEC 几何体积');axs[1].plot(t[1:],np.diff(true)/3600,c='#111111',lw=1.4)
 for name in MODELS:
  p=load(P/'predictions'/f'{event}_{name}.npz');mask=p['time']>=0;tm=p['time'][mask]/24;vv=np.cumsum(p['reach_bookkeeping_volume'],axis=0)[8,mask];vv-=vv[0];rr=np.cumsum(p['resampling_volume'],axis=0)[8,mask];rr-=rr[0];axs[0].plot(tm,vv/1e6,c=COLORS[name],ls=STYLES[name],label=NAMES[name]);axs[1].plot(tm[4::4],np.diff(vv[::4])/3600,c=COLORS[name],ls=STYLES[name]);axs[2].plot(tm,rr/1000,c=COLORS[name],ls=STYLES[name])
 for ax,title,unit in zip(axs,['有效储量异常的几何诊断','储量变化率诊断','独立时段衔接的数值体积'],['储量异常（百万 m³）','储量变化率（m³/s）','累计数值体积（千 m³）']):ax.set_title(title+'\n距出口 20 km / 距入口 180 km',fontsize=9);ax.set_ylabel(unit);ax.set_xlabel('预报起点后时间（天）');ax.grid(alpha=.15)
 axs[0].legend(fontsize=6.5);export(fig,'05_A_storage_and_interfaces')
 # Exact same predictor versus two different downstream boundary outcomes.
 fig,axs=plt.subplots(1,2,figsize=(10.3,3.8),layout='constrained');d=load(P/'data/E094.npz');keep=d['time']>=0
 for ax,x in zip(axs,[100,20]):
  j=int(np.argmin(abs(d['x']-x)));delta=d['q'][keep,j]-d['q_ref'][keep,j];ax.plot(d['time'][keep]/24,delta,c='#111111',label='HEC 动态减固定边界');ax.axhline(0,c=COLORS['original'],ls='--',label='四个入流模型的预测增量');ax.set_title(f'距出口 {x} km / 距入口 {200-x} km');ax.set_xlabel('预报起点后时间（天）');ax.set_ylabel('配对流量差（m³/s）');ax.grid(alpha=.15)
 axs[0].legend(fontsize=7);export(fig,'06_A_boundary_response_identifiability')
 for event in cfg['validation']:
  d=load(P/'data'/f'{event}.npz');keep=d['time']>=0;fig,axs=plt.subplots(1,3,figsize=(11.5,3.3),layout='constrained')
  for ax,x in zip(axs,[100,20,0]):
   col=int(np.argmin(abs(d['x']-x)));q=d['q'][keep,col];ax.plot(d['time'][keep]/24,q,c='#111111',lw=1.4,label='HEC-RAS')
   for name in MODELS:
    pp=load(P/'predictions'/f'{event}_{name}.npz');mask=pp['time']>=0;ax.plot(pp['time'][mask]/24,pp['q'][int((200-x)/20),mask],c=COLORS[name],ls=STYLES[name],lw=1.,label=NAMES[name])
   if q.min()<0:ax.axhline(0,c='.6',lw=.6)
   ax.set_title(f'{event}：距出口 {x} km / 距入口 {200-x} km',fontsize=9);ax.set_xlabel('预报起点后时间（天）');ax.set_ylabel('流量（m³/s）');ax.grid(alpha=.15)
  axs[0].legend(fontsize=6.5);export(fig,f'event_{event}_full_hydrographs')
if __name__=='__main__':main()
