"""Actual HEC/forecast evidence only; figure contract in FIGURE_CONTRACT.md."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib import font_manager

plt.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'lines.linewidth':1.2,'savefig.dpi':300,'legend.frameon':False})
COL={'RWEC-0':'#65757C','RWEC-L-predicted':'#0072B2','RWEC-L-oracle':'#CC79A7','Muskingum-linear':'#D55E00'}
LAB={'RWEC-0':'RWEC-0','RWEC-L-predicted':'RWEC-L predicted','RWEC-L-oracle':'RWEC-L oracle','Muskingum-linear':'Linear Muskingum'}
STYLE={'RWEC-0':'--','RWEC-L-predicted':'-','RWEC-L-oracle':':','Muskingum-linear':'-.'}
FIGS=[]

def axes2():
 fig,axs=plt.subplots(2,2,figsize=(180/25.4,142/25.4))
 fig.subplots_adjust(left=.10,right=.89,bottom=.15,top=.94,wspace=.48,hspace=.48)
 for ch,ax in zip('abcd',axs.flat):ax.text(-.15,1.06,ch,transform=ax.transAxes,fontweight='bold',fontsize=9)
 return fig,axs

def finish(fig,name):
 p=OUT/'figures'/name
 fig.savefig(str(p)+'.pdf')
 fig.savefig(str(p)+'.svg')
 fig.savefig(str(p)+'.png',dpi=300)
 fig.savefig(str(p)+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
 FIGS.append(name);plt.close(fig)

def packet(g,e='E006'):return load(OUT/'data'/f'{g}_{e}.npz')

def mechanisms():
 fig,axs=axes2()
 for c,g in enumerate(['A','B']):
  d=packet(g);f=GRID>=0;t=GRID[f]/24;ax=axs[0,c];v=6 if g=='A' else 2
  mesh=ax.pcolormesh(t,d['x'],d['delta_stage'][f].T,cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-v,vcenter=0,vmax=v),rasterized=True,shading='auto')
  ax.plot(t,d['true_lb'][f],color='black',label='Farthest > 0.10 m')
  ax.plot(t,d['true_connected'][f],color='#F0B442',ls='--',label='Outlet-connected')
  ax.set(xlabel='Time from forecast origin (d)',ylabel='Distance from outlet (km)',title=f'{g}: E006 水位差与回水末端',ylim=(0,200));ax.legend(loc='upper right',fontsize=7)
  cb=fig.colorbar(mesh,ax=ax,pad=.02,fraction=.046);cb.set_label('Water-level difference (m)')
  ax=axs[1,c];j=1;vol=d['dv'][:,-1]-d['dv'][:,10]
  ax.plot(t,vol[f]/1.e6,color='#0072B2',label='Additional volume')
  ay=ax.twinx();ay.plot(t,(d['truth'][:,0]-d['reference'][:,0])[f],color='#D55E00',ls='--',label='Flow difference')
  ay.plot(t,-np.gradient(vol,GRID*3600)[f],color='#111111',ls=':',label='Negative volume rate')
  ax.set(xlabel='Time from forecast origin (d)',ylabel='Additional volume (10⁶ m³)',title=f'{g}: 20 km断面以上的蓄入与释放');ay.set_ylabel('Flow difference (m³/s)')
  if c==0:
   handles,labs=ax.get_legend_handles_labels();hh,ll=ay.get_legend_handles_labels();fig.legend(handles+hh,labs+ll,loc='lower center',bbox_to_anchor=(.5,.01),ncol=3,fontsize=7)
 finish(fig,'01_endpoint_and_storage')

def endpoint():
 fig,axs=axes2();metrics=pd.read_csv(OUT/'tables/endpoint_prediction_metrics.csv')
 for c,g in enumerate(['A','B']):
  d=packet(g);p=load(OUT/'predictions'/f'endpoint_{g}_E006.npz');f=GRID>=0;ax=axs[0,c]
  ax.fill_between(GRID[f]/24,p['low'][f],p['high'][f],color='#56B4E9',alpha=.28,label='Empirical 10–90% range')
  ax.plot(GRID[f]/24,d['true_lb'][f],color='black',label='HEC-RAS reference')
  ax.plot(GRID[f]/24,p['central'][f],color='#0072B2',label='Predicted length')
  ax.set(xlabel='Time from forecast origin (d)',ylabel='Backwater extent (km)',title=f'{g}: E006 固定起点预测',ylim=(0,205 if g=='A' else 90))
  if c==0:fig.legend(*ax.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.01),ncol=3,fontsize=7)
  ev=metrics[(metrics.group==g)&(metrics.subset=='validation')].sort_values('event');ax=axs[1,c]
  x=np.arange(len(ev));ax.bar(x,ev.MAE_km,color='#0072B2',width=.65)
  ax.set_xticks(x);ax.set_xticklabels(ev.event,rotation=90,rotation_mode='anchor',ha='right');ax.set(ylabel='Length MAE (km)',title=f'{g}: 全部开发内验证洪水')
  ay=ax.twinx();ay.plot(x,ev.coverage*100,'o',color='#D55E00',ms=3,label='Range coverage');ay.set(ylim=(0,105),ylabel='Empirical coverage (%)')
 finish(fig,'02_fixed_origin_endpoint_prediction')

def comparison():
 fig,axs=axes2();d=pd.read_csv(OUT/'tables/flow_event_metrics.csv');d=d[(d.subset=='validation')&(d.window=='positive_common')]
 for c,g in enumerate(['A','B']):
  ev=d[d.group==g];ax=axs[0,c]
  for off,(m,color) in zip([-.9,-.3,.3,.9],COL.items()):
   mm=ev[ev.model==m].groupby('xc_km').NSE.agg(['mean','std']);ax.errorbar(mm.index+off,mm['mean'],yerr=mm['std'],color=color,ls=STYLE[m],lw=1.,marker='o',ms=2.5,capsize=1.5,label=LAB[m],elinewidth=.65)
  ax.set(xlabel='Distance from outlet (km)',ylabel='NSE (mean ± event SD)',title=f'{g}: 共同非负时刻的洪水演进精度');ax.legend(loc='lower right',fontsize=7)
  ax=axs[1,c]
  for k,(xc,marker,col) in enumerate([(20,'o','#0072B2'),(100,'s','#D55E00')]):
   sub=ev[ev.xc_km==xc].pivot(index='event',columns='model',values='NSE').sort_index();diff=sub['RWEC-L-predicted']-sub['RWEC-0'];xx=np.arange(len(diff))+(k-.5)*.18
   ax.scatter(xx,diff,marker=marker,s=12,color=col,label=f'{xc} km')
  ax.axhline(0,color='black',lw=.7);ax.set_xticks(np.arange(len(sub)));ax.set_xticklabels(sub.index,rotation=90,rotation_mode='anchor',ha='right');ax.set(ylabel='NSE change: predicted L − RWEC-0',title=f'{g}: 同场洪水的逐场增益');ax.legend(loc='best')
  if g=='B':ax.set_ylim(-.01,.006)
 finish(fig,'03_joint_calibrated_routing_comparison')

def hydrographs():
 fig,axs=axes2()
 for c,(g,eid) in enumerate([('A','E037'),('B','E006')]):
  d=packet(g,eid);f=GRID>=0
  for r,xc in enumerate([20,100]):
   ax=axs[r,c];j=int(np.flatnonzero(TARGETS==xc)[0]);ax.plot(GRID[f]/24,d['truth'][f,j],color='black',lw=1.6,label='HEC-RAS reference')
   for m in COL:
    p=load(OUT/'predictions'/f'flow_{g}_{eid}_{m}.npz');ax.plot(GRID[f]/24,p['q'][f,j],color=COL[m],ls=STYLE[m],lw=1.05,label=LAB[m])
   ax.axhline(0,color='#AAAAAA',lw=.6);ax.set(xlabel='Time from forecast origin (d)',ylabel='Flow (m³/s)',title=f'{g}: {eid}, x = {xc} km')
   if r==0 and c==0:ax.legend(loc='upper right',fontsize=7)
 finish(fig,'04_full_hydrograph_examples')

def waterlines():
 fig,axs=axes2()
 for c,g in enumerate(['A','B']):
  d=packet(g);ax=axs[0,c]
  for hour,col in [(24,'#0072B2'),(72,'#D55E00'),(120,'#009E73')]:
   i=int(np.flatnonzero(GRID==hour)[0]);ax.plot(d['x'],d['z'][i],color=col,label=f'{hour} h, dynamic');ax.plot(d['x'],d['zref'][i],color=col,ls='--',lw=.7)
  ax.set(xlabel='Distance from outlet (km)',ylabel='Water-surface elevation (m)',title=f'{g}: E006 沿程水面线');ax.legend(loc='upper left',fontsize=7)
  ax=axs[1,c];f=GRID>0;t=GRID[f]/24
  ax.plot(t,d['moving_Q'][f],color='#333333',label='Physical flow at endpoint')
  ax.plot(t,d['moving_Q'][f]+d['moving_sweep'][f],color='#0072B2',label='Relative moving-boundary flux')
  ax.axhline(0,color='#AAAAAA',lw=.7);ax.set(xlabel='Time from forecast origin (d)',ylabel='Flux (m³/s)',title=f'{g}: 移动末端的物理流量与相对通量')
  if c==0:fig.legend(*ax.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.01),ncol=2,fontsize=7)
 finish(fig,'05_waterlines_and_moving_boundary')

if __name__=='__main__':
 dirs();mechanisms();endpoint();comparison();hydrographs();waterlines()
 save(OUT/'provenance/figure_exports.json',dict(figures=FIGS,formats=['pdf','svg','png','tiff'],width_mm=180,font_english='Times New Roman',font_chinese='SimSun',min_text_pt=7,example_selection='A E037 retained prior case; A/B E006 first validation event, not error-ranked'))
 print('Exported',FIGS)
