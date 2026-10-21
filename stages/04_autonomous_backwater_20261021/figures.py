"""Quantitative evidence only; all event and section selection is declared."""
from common import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,FFMpegWriter
from prototype import forecast
import sys
mpl.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':8,'axes.titlesize':9,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'svg.fonttype':'none',
 'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':.7,'axes.unicode_minus':False,
 'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold'})
COL={'RWEC-0':'#787878','RWEC-L':'#247A9E','RWEC-V':'#C97B33','RWEC-H':'#9A83B3'}
IDS=['E023','E006','E007','E069']
def export(fig,name):
    fig.savefig(OUT/'figures'/f'{name}.png',dpi=300)
    fig.savefig(OUT/'figures'/f'{name}.pdf')
    fig.savefig(OUT/'figures'/f'{name}.svg')
    fig.savefig(OUT/'figures'/f'{name}.tiff',dpi=600)
    plt.close(fig)
def panel(ax,label,title):
    ax.set_title(f'{label}  {title}',loc='left',pad=7)
    ax.tick_params(direction='out',width=.6,length=3)
def main():
    dirs();lookup=loadnp(OUT/'data/hydraulic_lookup.npz');params=read(OUT/'fits/all_parameters.json')
    assert np.all(np.diff(lookup['q'])>0)
    data={e:loadnp(OUT/'data'/f'B_{e}.npz') for e in IDS}
    preds={(e,m):loadnp(OUT/'predictions'/f'{e}_mode{m}.npz') for e in IDS for m in [0,1,2,3]}
    fig,axs=plt.subplots(1,2,figsize=(180/25.4,78/25.4),layout='constrained')
    info=pd.read_csv(OUT/'tables/A_information_bounds.csv');i=info[info.group=='ten_boundary_scenarios']
    for eid,color in zip(IDS,['#247A9E','#C97B33','#5A9578','#9A83B3']):
        d=i[(i.event_id==eid)&(i.xc_km>0)]
        axs[0].plot(d.xc_km,d.optimal_pooled_RMSE_floor_m3s,'o-',ms=3,lw=1,color=color,label=eid)
    axs[0].set(xlabel='距出口位置 (km)',ylabel='不可消除的RMSE下界 (m³/s)');axs[0].legend(ncol=2)
    panel(axs[0],'a','相同洪水、十种下游边界')
    dg=pd.read_csv(OUT/'tables/A_geometric_diagnostic_metrics.csv');v=dg[(dg.subset=='validation')&(dg.xc_km==20)]
    dd=[v[v.diagnostic==s].NSE for s in ['length_only','boundary_height_times_length']]
    axs[1].boxplot(dd,labels=['仅有效长度','水位增量×有效长度'],widths=.5,showfliers=True,
       medianprops=dict(color='#247A9E'),flierprops=dict(marker='o',markersize=2))
    axs[1].axhline(0,color='.7',lw=.6);axs[1].set(ylabel='几何附加体积诊断NSE')
    panel(axs[1],'b','20 km：21场开发内验证')
    export(fig,'A_information_and_geometry')
    fig,axs=plt.subplots(2,2,figsize=(180/25.4,130/25.4),layout='constrained')
    for path in sorted((OUT/'data').glob('A_E*.npz')):
        ad=loadnp(path);axs[0,0].plot(ad['time_h'],ad['lb_reference_km'],color='#787878',alpha=.18,lw=.6)
    ad=loadnp(OUT/'data/A_E006.npz');at=ad['time_h'];al=ad['lb_reference_km'];av=ad['delta_volume_m3'][:,-1]-ad['delta_volume_m3'][:,2]
    axs[0,0].plot(at,al,color='#247A9E',lw=1.3,label='E006');axs[0,0].legend();axs[0,0].set(xlabel='洪水起点后时间 (h)',ylabel='配对参考长度 (km)');panel(axs[0,0],'a','A：70场开发洪水，15场触及200 km')
    ax=axs[0,1]
    ax.plot(at,al,color='#202020',label='总影响长度')
    for xc,color in [(20,'#247A9E'),(100,'#C97B33')]:ax.plot(at,np.maximum(al-xc,0),color=color,label=f'{xc} km的有效长度')
    ax.set(xlabel='洪水起点后时间 (h)',ylabel='长度 (km)');ax.legend(fontsize=6);panel(ax,'b','E006：有效长度随预测位置改变')
    ax=axs[1,0];ax.plot(at,av/1e6,color='#202020',label='几何附加体积')
    gain=pd.read_csv(OUT/'tables/A_geometric_diagnostic_parameters.csv');g=gain[(gain.xc_km==20)&(gain.diagnostic=='boundary_height_times_length')].gain.iloc[0]
    proxy=g*500*np.maximum(ad['boundary_increment_m'],0)*np.maximum(al-20,0)*1000
    ax.plot(at,proxy/1e6,color='#247A9E',label='水位×长度诊断估计')
    ax.set(xlabel='洪水起点后时间 (h)',ylabel='20—200 km附加体积 (百万m³)');ax.legend(fontsize=6);panel(ax,'c','同一长度不能唯一决定体积')
    ax=axs[1,1];eg=pd.read_csv(OUT/'tables/A_event_geometry.csv');eg=eg[eg.xc_km==20]
    for subset,color in [('train','#787878'),('validation','#C97B33')]:
        g=eg[eg.subset==subset];ax.scatter(g.max_lb_km,g.peak_delta_volume_m3/1e6,s=10,color=color,label='49场训练' if subset=='train' else '21场开发内验证')
    ax.set(xlabel='全过程最大参考长度 (km)',ylabel='20—200 km最大附加体积 (百万m³)');ax.legend(fontsize=6);panel(ax,'d','70场配对的峰值关系')
    export(fig,'A_length_storage_mechanism')
    fig,axs=plt.subplots(2,2,figsize=(180/25.4,130/25.4),layout='constrained')
    for ax,eid,label in zip(axs.flat,IDS,'abcd'):
        d=data[eid];t=d['time'];m=t>=0
        ax.plot(t[m],d['q_control'][m,10],color='#202020',lw=1.5,label='HEC固定控制')
        for name,mode in [('RWEC-0',0),('RWEC-L',1),('RWEC-V',2)]:
            ax.plot(t[m],preds[eid,mode]['q'][m,9],color=COL[name],lw=1,ls='--' if mode==0 else '-',label=name)
        ax.set(xlabel='洪水起点后时间 (h)',ylabel='20 km流量 (m³/s)')
        panel(ax,label,eid+(' 开发内验证' if eid=='E006' else ' 训练'))
    axs[0,0].legend(ncol=2,loc='upper right');export(fig,'B_four_flood_hydrographs')
    fig,axs=plt.subplots(2,2,figsize=(180/25.4,125/25.4),layout='constrained')
    for ax,eid,label in zip(axs.flat,IDS,'abcd'):
        d=data[eid];t=d['time'];m=t>=0
        ax.plot(t[m],d['lb'][m],color='#202020',lw=1.5,label='HEC配对参考')
        for mode,name in [(1,'RWEC-L'),(2,'RWEC-V')]:ax.plot(t[m],preds[eid,mode]['lb_km'][m],color=COL[name],lw=1,label=name)
        ax.set(xlabel='洪水起点后时间 (h)',ylabel='回水长度 (km)',ylim=(0,70))
        panel(ax,label,eid+(' 开发内验证' if eid=='E006' else ' 训练'))
    axs[0,0].legend(loc='upper right');export(fig,'B_dynamic_backwater_lengths')
    eid='E006';d=data[eid];p=preds[eid,2];t=d['time'];x=d['x']
    peak=int(np.argmax(d['q_control'][:,0]));times=[max(0,peak-20),peak,peak+20,peak+60]
    fig,axs=plt.subplots(2,2,figsize=(180/25.4,130/25.4),layout='constrained')
    for ax,k,label in zip(axs.flat,times,'abcd'):
        qp=p['q'][k,10];yr=np.array([np.interp(qp,lookup['q'],lookup['yref'][:,j]) for j in range(101)])
        model=x*.1+yr+p['delta_stage'][k]
        ax.plot(x,d['z_control'][k],color='#202020',lw=1.4,label='HEC固定控制')
        ax.plot(x,d['z_ref'][k],color='.6',ls='--',lw=1,label='HEC固定5 m')
        ax.plot(x,model,color=COL['RWEC-V'],lw=1,label='RWEC-V近似水面线')
        ax.plot(x,x*.1,color='.8',lw=.6,label='河床')
        ax.set(xlabel='距出口位置 (km)',ylabel='水面高程 (m)');panel(ax,label,f'E006：{t[k]:.0f} h')
    axs[0,0].legend(fontsize=6,loc='upper left');export(fig,'B_longitudinal_profiles')
    fig,axs=plt.subplots(1,2,figsize=(180/25.4,80/25.4),layout='constrained')
    f=pd.read_csv(OUT/'tables/B_event_metrics.csv');v=f[f.event=='E006']
    baseline=v[v.model=='RWEC-0'].set_index('xc_km').RMSE_m3s
    for name in ['RWEC-L','RWEC-V','RWEC-H']:
        s=v[v.model==name].set_index('xc_km').RMSE_m3s
        axs[0].plot(s.index,baseline-s,'o-',ms=3,lw=1,color=COL[name],label=name)
    axs[0].axhline(0,color='.7',lw=.6);axs[0].set(xlabel='距出口位置 (km)',ylabel='相对RWEC-0的RMSE减少 (m³/s)');axs[0].legend()
    panel(axs[0],'a','E006：长度因子的改善较小')
    ax=axs[1];ax.plot(t[t>=0],d['q_control'][t>=0,10],color='#202020',lw=1.4,label='HEC')
    for dt,ls in [(1.,'-'),(.5,'--'),(.25,':')]:
        e=read(B/'config/pilot_frozen.json')['events'][eid]
        sys.path.insert(0,str(B/'scripts'));import hec200
        fine=np.arange(t[0],t[-1]+dt/2,dt);inf=hec200.boundary(e,fine,'ref')[0]
        rr=forecast(inf,dt,params['core'],params['gains']['RWEC-V'],2,lookup)
        m=fine>=0;ax.plot(fine[m],rr[0][m,9],ls=ls,lw=1,label=f'RWEC-V {dt:g} h')
    ax.set(xlabel='洪水起点后时间 (h)',ylabel='20 km流量 (m³/s)');ax.legend(fontsize=6)
    panel(ax,'b','固定参数：时间步长尚不稳定');export(fig,'B_ablation_and_timestep')
    # One synchronized animation; frame interval3h, no scientific resampling.
    sys.path.insert(0,str(A/'analysis/03_internal_discovery_20261021/.deps'))
    import imageio_ffmpeg
    mpl.rcParams['animation.ffmpeg_path']=imageio_ffmpeg.get_ffmpeg_exe()
    fig,axs=plt.subplots(3,1,figsize=(180/25.4,170/25.4),layout='constrained')
    ax=axs[0];a,=ax.plot(x,d['z_control'][0],color='#202020',label='HEC固定控制');b,=ax.plot(x,d['z_ref'][0],color='.6',ls='--',label='固定5 m参照')
    c,=ax.plot(x,d['z_ref'][0],color=COL['RWEC-V'],label='RWEC-V近似')
    ax.plot(x,x*.1,color='.8',lw=.6);ax.set(xlabel='距出口位置 (km)',ylabel='水面高程 (m)',ylim=(0,27));ax.legend(loc='upper left',ncol=3)
    m=t>=0
    axs[1].plot(t[m],d['lb'][m],color='#202020',label='HEC参考');axs[1].plot(t[m],p['lb_km'][m],color=COL['RWEC-V'],label='RWEC-V')
    axs[1].set(xlabel='洪水起点后时间 (h)',ylabel='回水长度 (km)');axs[1].legend(loc='upper right')
    axs[2].plot(t[m],d['q_control'][m,10],color='#202020',label='HEC 20 km');axs[2].plot(t[m],p['q'][m,9],color=COL['RWEC-V'],label='RWEC-V 20 km')
    axs[2].set(xlabel='洪水起点后时间 (h)',ylabel='流量 (m³/s)');axs[2].legend(loc='upper right')
    v1=axs[1].axvline(0,color='#247A9E',lw=.8);v2=axs[2].axvline(0,color='#247A9E',lw=.8)
    title=axs[0].set_title('E006 开发内验证：0 h',loc='left')
    def update(k):
        a.set_ydata(d['z_control'][k]);b.set_ydata(d['z_ref'][k])
        qp=p['q'][k,10];yr=np.array([np.interp(qp,lookup['q'],lookup['yref'][:,j]) for j in range(101)])
        c.set_ydata(x*.1+yr+p['delta_stage'][k]);v1.set_xdata([t[k],t[k]]);v2.set_xdata([t[k],t[k]])
        title.set_text(f'E006 开发内验证：{t[k]:.0f} h');return a,b,c,v1,v2,title
    frames=np.flatnonzero((t>=0)&(t<=260))[::3]
    ani=FuncAnimation(fig,update,frames=frames,blit=False)
    ani.save(OUT/'animation/E006_synchronized.mp4',writer=FFMpegWriter(fps=8,bitrate=1600),dpi=300)
    update(int(frames[len(frames)//2]));fig.savefig(OUT/'figures/animation_preview.png',dpi=300);plt.close(fig)
    print('Figures and synchronized animation complete',flush=True)
if __name__=='__main__':main()
