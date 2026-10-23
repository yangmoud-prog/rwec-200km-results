from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager
OUT=Path(__file__).resolve().parent
fontManager.addfont(r'C:\Windows\Fonts\times.ttf');fontManager.addfont(r'C:\Windows\Fonts\simsun.ttc')
mpl.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':7.4,'axes.labelsize':7.4,'axes.titlesize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':6.8,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'lines.linewidth':1.1,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold'})
C={'RWEC0':'#595A7B','linear_fixed':'#C8823C','linear_joint':'#397D82'}
L={'RWEC0':'RWEC-0','linear_fixed':'RWEC-1 固定核心','linear_joint':'RWEC-1 联合率定'}
S=pd.read_csv(OUT/'tables/window_summary.csv');M=pd.read_csv(OUT/'tables/event_metrics.csv');D=pd.read_csv(OUT/'tables/negative_error_summary.csv')
def export(fig,name):
 fig.savefig(OUT/'figures'/f'{name}.pdf')
 fig.savefig(OUT/'figures'/f'{name}.svg')
 fig.savefig(OUT/'figures'/f'{name}.png',dpi=300)
 fig.savefig(OUT/'figures'/f'{name}.tiff',dpi=600)
 plt.close(fig)
def panel(ax,letter,title):ax.set_title(f'{letter}  {title}',loc='left',fontweight='bold',pad=7)

fig,axs=plt.subplots(2,2,figsize=(180/25.4,132/25.4));fig.subplots_adjust(left=.085,right=.97,bottom=.1,top=.92,wspace=.32,hspace=.52)
ax=axs[0,0]
for wn,lab,ls,col in [('reference_all','固定5 m参照','-.','#414141'),('all','变水位：全部时刻','--','#9A9AA5'),('nonnegative','变水位：仅正流','-','#595A7B')]:
 z=S[(S.model=='RWEC0')&(S.window==wn)].sort_values('xc_km');ax.plot(z.xc_km,z.NSE,ls,marker='o',ms=3,color=col,label=lab)
ax.set(xlabel='距下游出口距离 (km)',ylabel='平均 NSE',ylim=(.76,1.005));ax.legend(loc='lower right',frameon=False);panel(ax,'a','排除倒流后仍存在空间梯度')
ax=axs[0,1];rng=np.random.default_rng(20261023)
for name in C:
 z=M[(M.model==name)&(M.window=='nonnegative')].pivot(index='event',columns='xc_km',values='NSE');a=z.to_numpy();boots=a[rng.integers(0,len(a),(10000,len(a)))].mean(1)
 lo,hi=np.quantile(boots,[.025,.975],axis=0);x=z.columns.to_numpy()
 ax.fill_between(x,lo,hi,color=C[name],alpha=.08,lw=0);ax.plot(x,a.mean(0),marker='o',ms=2.5,color=C[name],label=L[name])
ax.set(xlabel='距下游出口距离 (km)',ylabel='正流时段平均 NSE',ylim=(.74,1.01));ax.legend(frameon=False,loc='lower right');panel(ax,'b','已知未来水位仍未稳定改善')
ax=axs[1,0];x=np.arange(len(D));ax.bar(x-.17,D.negative_sample_pct,width=.34,color='#B9BCC8',label='负流量小时占比');ax.bar(x+.17,D.pooled_negative_SSE_pct,width=.34,color='#C8823C',label='负流量时刻平方误差占比')
ax.set(xticks=x,xticklabels=D.xc_km.astype(int),xlabel='距下游出口距离 (km)',ylabel='占比 (%)',ylim=(0,52));ax.legend(frameon=False,loc='upper right');panel(ax,'c','倒流影响大但不能解释全部误差')
ax=axs[1,1];z=M[(M.xc_km==100)&(M.window=='nonnegative')].pivot(index='event',columns='model',values='NSE');z=z.assign(gain=z.linear_joint-z.RWEC0).sort_values('gain')
ax.axhline(0,color='#808080',lw=.65)
for name,dx in [('linear_fixed',-.1),('linear_joint',.1)]:ax.scatter(np.arange(len(z))+dx,z[name]-z.RWEC0,s=12,color=C[name],label=L[name],alpha=.85)
ax.set(xlabel='21场验证洪水（按联合修正增益排序）',ylabel='相对 RWEC-0 的 NSE 变化',xticks=[0,10,20],xticklabels=['1','11','21']);ax.legend(frameon=False,loc='upper left');panel(ax,'d','100 km平均改善由少数事件主导')
export(fig,'positive_flow_comparison')

best=z.index[-1];worst=z.index[0]
fig,axs=plt.subplots(2,2,figsize=(180/25.4,127/25.4));fig.subplots_adjust(left=.09,right=.97,bottom=.095,top=.88,wspace=.3,hspace=.5)
for row,eid in enumerate([best,worst]):
 with np.load(OUT/'data'/f'{eid}.npz') as zz:d={k:zz[k] for k in zz.files}
 t=d['time'];mask=t>=0;days=t[mask]/24;j=4
 ax=axs[row,0];ax.plot(days,d['truth'][mask,j],color='#232323',lw=1.4,label='HEC-RAS')
 for name in C:ax.plot(days,d[name][mask,j],color=C[name],label=L[name],ls={'RWEC0':'--','linear_fixed':':','linear_joint':'-'}[name])
 ax.set(xlabel='洪水起点后时间 (天)',ylabel='100 km断面流量 (m³/s)');panel(ax,chr(97+row*2),f'{eid} 整体洪水过程')
 ax=axs[row,1];ax.plot(days,(d['truth']-d['reference'])[mask,j],color='#232323',lw=1.4)
 for name in C:ax.plot(days,(d[name]-d[name+'_ref'])[mask,j],color=C[name],ls={'RWEC0':'--','linear_fixed':':','linear_joint':'-'}[name])
 ax.axhline(0,color='#A0A0A0',lw=.5,zorder=0);ax.set(xlabel='洪水起点后时间 (天)',ylabel='变动水位减固定水位流量 (m³/s)');panel(ax,chr(98+row*2),f'{eid} 回水增量响应')
handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,ncol=4,loc='upper center',bbox_to_anchor=(.52,.99),frameon=False)
export(fig,'x100_total_and_backwater_response')
(OUT/'provenance/figure_examples.json').write_text(__import__('json').dumps(dict(best_gain_event=best,worst_gain_event=worst,selection='Maximum/minimum validation NSE gain of frozen joint linear model at100 km; illustrations only'),indent=2),encoding='utf-8')
print('Two figures exported with all 21 events in the summary; example events',best,worst)
