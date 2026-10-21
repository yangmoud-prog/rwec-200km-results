from common import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
mpl.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,
 'legend.fontsize':6,'legend.frameon':False,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.right':False,'axes.spines.top':False,'axes.unicode_minus':False,
 'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold'})
C={'RWEC-0':'#737373','RWEC-L':'#247A9E','RWEC-V':'#C97B33','RWEC-H':'#9A83B3','RWEC-1_oracle_joint':'#8B5B83'}
def export(fig,name,gallery=False):
    folder=OUT/'figures';folder.mkdir(exist_ok=True)
    fig.savefig(folder/f'{name}.png',dpi=300);fig.savefig(folder/f'{name}.pdf')
    if not gallery:
        fig.savefig(folder/f'{name}.svg');fig.savefig(folder/f'{name}.tiff',dpi=600)
    plt.close(fig)
def setup(ax,label,title,xlabel='洪水起点后时间 (h)',ylabel='流量 (m³/s)'):
    ax.set_title(f'{label}  {title}',loc='left',pad=7);ax.set_xlabel(xlabel);ax.set_ylabel(ylabel)
    ax.tick_params(direction='out',width=.6,length=3)
def pred(label,eid,name):return loadnp(OUT/'predictions'/f'{label}_{eid}_{name}.npz')
def a_truth(eid):
    plan=read(A/'config/formal_plan.json');cases=[c for c in plan['formal_cases'] if c['event_id']==eid]
    packets=[]
    for c in cases:
        d=loadnp(A/'results'/c['case']/'reach_timeseries.npz');packets.append((c,d))
    return packets
def A_figures():
    fig,axs=plt.subplots(2,2,figsize=(180/25.4,130/25.4),layout='constrained')
    for ax,eid,letter in zip(axs.flat,['E023','E006','E007','E069'],'abcd'):
        raw=a_truth(eid);nt=min(len(d['time_since_flood_h']) for c,d in raw);t=raw[0][1]['time_since_flood_h'][:nt];m=t>=0
        q=np.array([d['flow_m3s'][:nt,2] for c,d in raw]);ax.fill_between(t[m],q[:,m].min(0),q[:,m].max(0),color='.85',label='10种真实边界响应范围')
        for name,mode in [('RWEC-L',1),('RWEC-V',2)]:
            p=loadnp(OUT/'predictions'/f'A_{eid}_nonzero_mode{mode}.npz');ax.plot(t[m],p['q'][:nt][m,9],color=C[name],lw=1,label=f'{name} α=0.01')
        p=pred('A',eid,'RWEC-0');ax.plot(t[m],p['q'][:nt][m,9],color='#202020',lw=1,label='RWEC-0／最优L、V')
        setup(ax,letter,eid+(' 开发内验证，20 km' if eid=='E006' else ' 训练背景，20 km'))
    axs[0,0].legend(loc='upper right');export(fig,'A_same_inflow_unknown_boundary')
    df=pd.read_csv(OUT/'tables/A_event_metrics.csv');stress=pd.read_csv(OUT/'tables/A_nonzero_gain_stress_metrics.csv')
    base=df[(df.subset=='validation')&(df.condition=='independent_stage_main')&(df.model=='RWEC-0')]
    fig,axs=plt.subplots(1,2,figsize=(180/25.4,82/25.4),layout='constrained')
    for name,d,color,ls in [('RWEC-0／最优L、V',base,'#202020','-'),('L，α=0.01',stress[(stress.subset=='validation')&(stress.condition=='independent_stage_main')&(stress.model=='RWEC-L_alpha0.01')],C['RWEC-L'],'--'),('V，α=0.01',stress[(stress.subset=='validation')&(stress.condition=='independent_stage_main')&(stress.model=='RWEC-V_alpha0.01')],C['RWEC-V'],':')]:
        g=d.groupby('xc_km').NSE;med=g.median();lo=g.quantile(.25);hi=g.quantile(.75)
        axs[0].plot(med.index,med,ls=ls,color=color,label=name,lw=1);axs[0].fill_between(med.index,lo,hi,color=color,alpha=.08)
    setup(axs[0],'a','21场验证：中位数与四分位距','距出口位置 (km)','NSE');axs[0].legend(loc='lower right')
    grid=pd.read_csv(OUT/'tables/A_alpha_grid_by_event.csv');grid=grid[grid.subset=='train']
    for mode,name in [(1,'RWEC-L'),(2,'RWEC-V')]:
        g=grid[grid['mode']==mode].groupby('alpha').loss_mean_without_constant.mean();axs[1].plot(g.index,g-g.iloc[0],'o-',ms=3,lw=1,color=C[name],label=name)
    setup(axs[1],'b','非零修正的训练损失变化','预先给定的α','相对α=0的损失增加');axs[1].legend();export(fig,'A_actual_routing_and_nonzero_gains')
    # All21 main-dynamic validation floods; no event omitted for low accuracy.
    for eid in read(P4/'provenance/event_split.json')['validation']:
        raw=next(d for c,d in a_truth(eid) if c['role']=='main' and c['kind']=='var');t=raw['time_since_flood_h'];m=t>=0;nt=len(t)
        base=pred('A',eid,'RWEC-0');lp=loadnp(OUT/'predictions'/f'A_{eid}_nonzero_mode1.npz');vp=loadnp(OUT/'predictions'/f'A_{eid}_nonzero_mode2.npz')
        fig,axs=plt.subplots(3,2,figsize=(180/25.4,175/25.4),layout='constrained')
        for ax,xc,letter in zip(axs.flat,TARGETS,'abcdef'):
            ax.plot(t[m],raw['flow_m3s'][m,xc//10],color='#202020',lw=1.2,label='HEC真实输出')
            ax.plot(t[m],base['q'][:nt][m,(200-xc)//20],color=C['RWEC-0'],lw=1,label='RWEC-0／最优L、V')
            ax.plot(t[m],lp['q'][:nt][m,(200-xc)//20],color=C['RWEC-L'],ls='--',lw=.9,label='L，α=0.01')
            ax.plot(t[m],vp['q'][:nt][m,(200-xc)//20],color=C['RWEC-V'],ls=':',lw=.9,label='V，α=0.01')
            setup(ax,letter,f'{eid}：{xc} km')
        axs[0,0].legend();export(fig,f'A_validation_{eid}',gallery=True)
    print('A static figures and all21 validation flood pages complete',flush=True)
def B_figures():
    frozen=read(B/'config/pilot_frozen.json');ids=frozen['order'];data={e:loadnp(OUT/'data'/f'B_{e}.npz') for e in ids};response=[]
    for eid,d in data.items():
        diff=d['q_control'][:,10]-d['q_ref'][:,10];response.append(dict(event=eid,subset='train' if eid in frozen['train'] else 'validation',Qp_m3s=frozen['events'][eid]['Qp_m3s'],D_h=frozen['events'][eid]['D_h'],r=frozen['events'][eid]['r'],
            outlet_stage_rise_m=float(d['z_control'][:,0].max()-5),peak_reference_lb_km=float(d['lb'].max()),min_flow_m3s=float(d['q_control'].min()),
            suppression20_m3s=float(diff.min()),recovery20_m3s=float(diff.max()),max_extra_volume20_m3=float((d['dv'][:,-1]-d['dv'][:,10]).max())))
    frame('B_20pilot_response.csv',response);r=pd.DataFrame(response)
    fig,axs=plt.subplots(1,2,figsize=(180/25.4,80/25.4),layout='constrained')
    axs[0].axhspan(1,4,color='.93',label='预定涨幅目标')
    for subset,color,marker in [('train','#737373','o'),('validation','#C97B33','s')]:
        s=r[r.subset==subset];axs[0].scatter(s.Qp_m3s,s.outlet_stage_rise_m,s=16,color=color,marker=marker,label='14场训练' if subset=='train' else '6场开发内验证')
        axs[1].scatter(s.Qp_m3s,s.peak_reference_lb_km,s=16,color=color,marker=marker)
    setup(axs[0],'a','全部20场动态水位响应','上游洪峰 (m³/s)','出口涨幅 (m)');axs[0].legend()
    setup(axs[1],'b','实际回水作用范围','上游洪峰 (m³/s)','全过程最大参考长度 (km)');export(fig,'B_twenty_pilot_hydraulics')
    for variable,name,ylab in [('q','B_six_validation_flows','20 km流量 (m³/s)'),('lb','B_six_validation_lengths','回水长度 (km)')]:
        fig,axs=plt.subplots(3,2,figsize=(180/25.4,165/25.4),layout='constrained')
        for ax,eid,letter in zip(axs.flat,frozen['validation'],'abcdef'):
            d=data[eid];t=d['time'];m=t>=0
            ax.plot(t[m],d['q_control'][m,10] if variable=='q' else d['lb'][m],color='#202020',lw=1.3,label='HEC参考')
            for model in ['RWEC-0','RWEC-L','RWEC-V']:
                if variable=='lb' and model=='RWEC-0':continue
                p=pred('B',eid,model);v=p['q'][:,9] if variable=='q' else p['lb_km']
                ax.plot(t[m],v[m],color=C[model],lw=.9,ls='--' if model=='RWEC-0' else '-',label=model)
            setup(ax,letter,eid,ylabel=ylab)
        axs[0,0].legend(ncol=2);export(fig,name)
    # All20 floods, each with flow and dynamic-length panels.
    for eid in ids:
        d=data[eid];t=d['time'];m=t>=0
        fig,axs=plt.subplots(1,2,figsize=(180/25.4,80/25.4),layout='constrained')
        for ax,var,letter in [(axs[0],'q','a'),(axs[1],'lb','b')]:
            ax.plot(t[m],d['q_control'][m,10] if var=='q' else d['lb'][m],color='#202020',lw=1.3,label='HEC参考')
            for model in ['RWEC-0','RWEC-L','RWEC-V']:
                if var=='lb' and model=='RWEC-0':continue
                p=pred('B',eid,model);ax.plot(t[m],p['q'][m,9] if var=='q' else p['lb_km'][m],color=C[model],lw=.9,ls='--' if model=='RWEC-0' else '-',label=model)
            setup(ax,letter,eid+(' 验证' if eid in frozen['validation'] else ' 训练'),ylabel='20 km流量 (m³/s)' if var=='q' else '回水长度 (km)')
        axs[0].legend(ncol=2);export(fig,f'B_flood_{eid}',gallery=True)
    # One declared validation event: absolute profiles and all model approximations.
    d=data['E006'];p=pred('B','E006','RWEC-V');lookup=loadnp(OUT/'data/B_lookup.npz');assert np.all(np.diff(lookup['q'])>0)
    t=d['time'];x=d['x'];peak=np.argmax(d['q_control'][:,0]);times=[peak-20,peak,peak+20,peak+60]
    fig,axs=plt.subplots(2,2,figsize=(180/25.4,130/25.4),layout='constrained')
    for ax,k,letter in zip(axs.flat,times,'abcd'):
        y=np.array([np.interp(p['q'][k,10],lookup['q'],lookup['yref'][:,j]) for j in range(101)])
        ax.plot(x,d['z_control'][k],color='#202020',lw=1.3,label='HEC固定控制');ax.plot(x,d['z_ref'][k],color='.65',ls='--',label='固定5 m参照')
        ax.plot(x,x*.075+y+p['delta_stage'][k],color=C['RWEC-V'],lw=1,label='RWEC-V近似');ax.plot(x,x*.075,color='.8',lw=.6,label='河床')
        setup(ax,letter,f'E006：{t[k]:.0f} h','距出口位置 (km)','水面高程 (m)')
    axs[0,0].legend();export(fig,'B_longitudinal_profiles')
    # Robustness: paired error improvement and fixed-parameter time-step failure.
    fig,axs=plt.subplots(1,2,figsize=(180/25.4,82/25.4),layout='constrained')
    df=pd.read_csv(OUT/'tables/B_event_metrics.csv');df=df[df.subset=='validation'];baseline=df[df.model=='RWEC-0'].set_index(['event','xc_km']).RMSE_m3s
    for model in ['RWEC-L','RWEC-V','RWEC-H']:
        z=df[df.model==model].set_index(['event','xc_km']).RMSE_m3s;paired=(baseline-z).reset_index(name='improvement');g=paired.groupby('xc_km').improvement
        center=g.mean();axs[0].plot(center.index,center,color=C[model],lw=1,label=model);axs[0].fill_between(center.index,g.quantile(.25),g.quantile(.75),color=C[model],alpha=.1)
    axs[0].axhline(0,color='.7',lw=.6);setup(axs[0],'a','6场验证：均值与四分位距','距出口位置 (km)','相对RWEC-0的RMSE减少 (m³/s)');axs[0].legend()
    s=pd.read_csv(OUT/'tables/B_timestep_sensitivity.csv');s=s[s.xc_km==20]
    for model in ['RWEC-0','RWEC-L','RWEC-V']:
        values=[df[(df.model==model)&(df.xc_km==20)].NSE.to_numpy()]+[s[(s.model==model)&(s.dt_h==dt)].NSE.to_numpy() for dt in [.5,.25]]
        center=[np.median(v) for v in values];axs[1].plot(range(3),center,'o-',ms=3,lw=1,color=C[model],label=model)
        axs[1].fill_between(range(3),[np.quantile(v,.25) for v in values],[np.quantile(v,.75) for v in values],color=C[model],alpha=.1)
    axs[1].set_xticks(range(3),['1','0.5','0.25']);setup(axs[1],'b','固定参数：中位数与四分位距','演算步长 (h)','20 km NSE');axs[1].legend();export(fig,'B_ablation_and_timestep')
    print('B static results plus all20 flood pages complete',flush=True)
if __name__=='__main__':A_figures() if sys.argv[1]=='A' else B_figures()
