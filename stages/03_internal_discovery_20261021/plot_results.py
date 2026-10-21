"""Chinese evidence figures. User-requested Times New Roman / SimSun override."""
from pathlib import Path
import json,math
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from prepare import OUT,read,save,frame
from fit import load
from evaluate import configs

plt.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':8,
    'axes.titlesize':9,'axes.labelsize':8,'legend.fontsize':7.3,'xtick.labelsize':7.5,'ytick.labelsize':7.5,
    'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65,'lines.linewidth':1.1,
    'legend.frameon':False,'svg.fonttype':'none','pdf.fonttype':42,'axes.unicode_minus':True,
    'savefig.facecolor':'white'})
WIDTH=180/25.4
COLORS={'RWEC0':'#8A8A8A','RWEC0J':'#4D4D4D','RWEC1':'#7992A6','P1':'#9A82A4','P2':'#3E8B83',
    'NN20261021':'#CF8D40','NN20261022':'#CF8D40','SR1':'#B74F4B','SR2':'#704F75','SR3':'#416D9B'}
STYLES={'RWEC0':'--','RWEC0J':'-.','RWEC1':':','P1':'-.','P2':'-','NN20261021':'-', 'NN20261022':'--','SR1':'--','SR2':'-.','SR3':':'}
audit=[]
def export(fig,name,claim,panels):
    fig.canvas.draw()
    fig.savefig(OUT/'figures'/f'{name}.pdf')
    fig.savefig(OUT/'figures'/f'{name}.svg')
    fig.savefig(OUT/'figures'/f'{name}.png',dpi=220)
    fig.savefig(OUT/'figures'/f'{name}.tiff',dpi=600)
    audit.append(dict(figure=name,claim=claim,archetype='quantitative grid',panels=panels,width_mm=180,
        font_override='Times New Roman Latin, SimSun Chinese',min_configured_font_pt=7.3,
        data_exclusion='none; model comparisons omit second seed only in shared legends; seed results retained in tables',
        formats=['PDF','SVG','600dpi TIFF','220dpi PNG']))
    plt.close(fig);print('中文图',name,flush=True)

def panel(ax,s):ax.text(-.16,1.06,s,transform=ax.transAxes,fontsize=10,fontweight='bold')
def grid(ax):ax.grid(axis='y',color='#EAEAEA',lw=.5,zorder=0)
def main():
    d,meta=load();M=pd.read_csv(OUT/'tables/case_metrics.csv');models=[c for c in configs() if c['type']!='neural' or c['primary']]
    names={c['id']:c['label'] for c in models}
    pred={c['id']:np.load(OUT/'predictions'/f'{c["id"]}.npz')['q'] for c in models}
    extra={c['id']:np.load(OUT/'predictions'/f'{c["id"]}.npz')['b'] for c in models}
    primary=read(OUT/'fits/neural_gate.json')['chosen_seed_by_training_loss'];nn=f'NN{primary}'
    sr=next((c['id'] for c in models if c['type']=='symbolic'),None)
    compact=['RWEC0','P2',nn]+([sr] if sr else [])
    def case(eid,kind='var'):
        return next(i for i,m in enumerate(meta) if m['role']=='main' and m['event_id']==eid and m['kind']==kind)
    def curves(ax,j,l,ids,reference=True):
        N=d['ns'][j];t=np.arange(N)-48
        ax.plot(t,d['obs'][j,:N],color='#222222',lw=1.6,label='HEC实际流量',zorder=4)
        if reference:
            r=case(meta[j]['event_id'],'ref');n=min(N,d['ns'][r]);ax.plot(t[:n],d['obs'][r,:n],color='#BBBBBB',ls=':',lw=1,label='固定5 m参照')
        for name in ids:ax.plot(t,pred[name][l,j,:N],color=COLORS[name],ls=STYLES[name],label=names[name])
        ax.set(xlim=(t[0],t[-1]),xlabel='距洪水起点时间（h）',ylabel='流量（m³/s）');grid(ax)
    def legend(fig,ax,ncol=3):
        handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=ncol,bbox_to_anchor=(.5,.005))

    S=pd.read_csv(OUT/'tables/section_screen.csv');R=pd.read_csv(OUT/'tables/paired_hydrodynamic_response.csv')
    fig,axs=plt.subplots(2,2,figsize=(WIDTH,140/25.4));fig.subplots_adjust(left=.11,right=.97,top=.92,bottom=.11,hspace=.54,wspace=.42)
    ax=axs[0,0];ax.plot(S.x_km,S.min_q_m3s,'o-',ms=3,color='#566F82');ax.axhline(0,color='#333333',ls='--',lw=.7)
    ax.scatter([100],[S[S.x_km==100].min_q_m3s.iloc[0]],s=30,color='#B74F4B',zorder=5)
    ax.set(title='严格负流量门槛：100 km首次合格',xlabel='距外部边界距离（km）',ylabel='全部172工况最小流量（m³/s）')
    ax=axs[0,1];ax.bar(S.x_km,S.negative_cases,width=6,color='#8DA3B1');ax.set(title='负流量工况数',xlabel='距外部边界距离（km）',ylabel='出现负值的工况数');grid(ax)
    for ax,col,title,ylabel in [(axs[1,0],'max_stage_increment_m','合格断面仍有水位响应','配对最大水位增量（m）'),
        (axs[1,1],'rms_q_difference_m3s','工况A流量减工况B流量','配对流量差均方根（m³/s）')]:
        g=R.groupby('x_km')[col];med=g.median();lo=g.min();hi=g.max()
        ax.fill_between(med.index,lo,hi,color='#E6ECEF',label='70场最小至最大');ax.plot(med.index,med,color='#3C7C9C',label='70场中位数')
        ax.axvline(100,color='#B74F4B',ls=':',lw=.8);ax.set(title=title,xlabel='距外部边界距离（km）',ylabel=ylabel);grid(ax)
    axs[1,0].legend(loc='upper right',fontsize=6.8)
    for ax,s in zip(axs.flat,'abcd'):panel(ax,s)
    export(fig,'01_section_selection','100 km兼具已保存时序正向流与可辨识回水响应',4)

    fig,axs=plt.subplots(2,2,figsize=(WIDTH,145/25.4));fig.subplots_adjust(left=.10,right=.98,top=.94,bottom=.17,hspace=.5,wspace=.30)
    for ax,eid,s in zip(axs.flat,['E007','E006','E025','E094'],'abcd'):
        j=case(eid);curves(ax,j,9,compact);ax.set_title(f'{eid}｜'+('训练' if meta[j]['subset']=='train' else '开发内验证')+'｜演算100 km');panel(ax,s)
    legend(fig,axs[0,0],3)
    export(fig,'02_full_hydrographs','代表过程与退化案例均保留完整前后时段',4)

    fig,axs=plt.subplots(1,2,figsize=(WIDTH,85/25.4));fig.subplots_adjust(left=.11,right=.98,top=.88,bottom=.28,wspace=.35)
    v=M[(M.subset=='validation')&(M.role=='main')&(M.kind=='var')]
    for c in models:
        g=v[v.model==c['id']].groupby('Lr_km')
        for ax,col in zip(axs,['NSE','RMSE']):
            med=g[col].median();low=g[col].quantile(.25);high=g[col].quantile(.75)
            ax.fill_between(med.index,low,high,color=COLORS[c['id']],alpha=.10,lw=0)
            ax.plot(med.index,med,color=COLORS[c['id']],ls=STYLES[c['id']],label=c['label'])
    for ax,title,y,s in zip(axs,['变动水位下的过程拟合','变动水位下的流量误差'],['NSE','RMSE（m³/s）'],'ab'):
        ax.set(title=title,xlabel='有效演算河长（km）',ylabel=y);grid(ax);panel(ax,s)
    legend(fig,axs[0],4)
    export(fig,'03_accuracy_by_length','按21场洪水比较有效演算长度效应；线为中位数，带为四分位范围',2)

    fig,axs=plt.subplots(2,2,figsize=(WIDTH,145/25.4));fig.subplots_adjust(left=.12,right=.98,top=.93,bottom=.19,hspace=.58,wspace=.35)
    scenarios=[('main',None,None,'变动水位流量减固定5 m参照流量'),
        ('幅度',5,1,'工况A（高幅）流量减工况B（低幅）流量'),
        ('时机',3,1,'工况A（较晚）流量减工况B（较早）流量'),
        ('维持',2,1,'工况A（维持）流量减工况B（不维持）流量')]
    for ax,(name,a,b,title),s in zip(axs.flat,scenarios,'abcd'):
        if name=='main':a=case('E007');b=case('E007','ref')
        else:
            a=next(j for j,m in enumerate(meta) if m['event_id']=='E007' and f'_C{a:02d}_' in m['case'])
            b=next(j for j,m in enumerate(meta) if m['event_id']=='E007' and f'_C{b:02d}_' in m['case'])
        N=min(d['ns'][a],d['ns'][b]);t=np.arange(N)-48
        ax.plot(t,d['obs'][a,:N]-d['obs'][b,:N],color='#222222',lw=1.6,label='HEC实际差值')
        for key in compact:ax.plot(t,pred[key][9,a,:N]-pred[key][9,b,:N],color=COLORS[key],ls=STYLES[key],label=names[key])
        ax.axhline(0,color='#777777',lw=.5);ax.set(title=title.replace('减工况B','减\n工况B'),xlim=(t[0],t[-1]),xlabel='距洪水起点时间（h）',ylabel='流量差（m³/s）');panel(ax,s);grid(ax)
    legend(fig,axs[0,0],3);fig.text(.5,.112,'差值为工况A流量减工况B流量；正值表示A增流，负值表示A抑流。背景洪水均为E007。',ha='center',fontsize=7)
    export(fig,'04_paired_boundary_responses','候选必须解释同一洪水的边界变化差值，而非只拟合洪峰',4)

    fig,axs=plt.subplots(2,2,figsize=(WIDTH,145/25.4));fig.subplots_adjust(left=.12,right=.98,top=.93,bottom=.19,hspace=.52,wspace=.34)
    ids=['P1','P2',nn]+([sr] if sr else [])
    for ax,eid,s in zip(axs[0],['E025','E094'],'ab'):
        curves(ax,case(eid),9,ids,False);ax.set_title(f'{eid}｜物理候选与学习模块');panel(ax,s)
    j=case('E025');N=d['ns'][j];t=np.arange(N)-48
    for key in ids:axs[1,0].plot(t,extra[key][9,j,:N],color=COLORS[key],ls=STYLES[key],label=names[key])
    axs[1,0].set(title='候选附加储量是模型状态异常',xlabel='距洪水起点时间（h）',ylabel='附加储量（百万m³）');grid(axs[1,0]);panel(axs[1,0],'c')
    phase_data=pd.read_csv(OUT/'tables/phase_metrics.csv');labels=['水位上涨','高位维持','水位消落','边界恢复后']
    for key in ids:
        g=phase_data[(phase_data.model==key)&(phase_data.subset=='validation')&(phase_data.role=='main')&(phase_data.kind=='var')&(phase_data.Lr_km==100)].groupby('phase').RMSE
        med=g.median().reindex(labels);axs[1,1].plot(range(4),med,marker='o',ms=2.5,color=COLORS[key],ls=STYLES[key],label=names[key])
    axs[1,1].set(title='开发内验证的分阶段误差',xticks=range(4),xticklabels=['上涨','维持','消落','恢复'],ylabel='RMSE中位数（m³/s）');grid(axs[1,1]);panel(axs[1,1],'d')
    legend(fig,axs[0,0],3)
    export(fig,'05_physical_and_learned_states','学习的是附加蓄泄状态而非直接替换全部流量预测',4)

    fig,axs=plt.subplots(2,2,figsize=(WIDTH,148/25.4));fig.subplots_adjust(left=.12,right=.96,top=.93,bottom=.15,hspace=.65,wspace=.80)
    p=pd.read_csv(OUT/'tables/P2_parameter_profile.csv');pp=p.pivot(index='tb_h',columns='alpha',values='train_loss')
    image=axs[0,0].pcolormesh(pp.columns,pp.index,pp.values,shading='auto',cmap='cividis')
    axs[0,0].set(yscale='log',xlabel='附加储量系数',ylabel='响应时间（h）',title='P2训练损失剖面');fig.colorbar(image,ax=axs[0,0],label='平均1−NSE',fraction=.05,pad=.04)
    ts=pd.read_csv(OUT/'tables/timestep_sensitivity.csv');ts=ts[(ts.comparison=='vs_1h')&(ts.dt_h==.25)]
    vals=[];errs=[]
    for c in models:
        u=ts[ts.model==c['id']].RMSE_vs_1h_m3s;vals.append(u.mean());errs.append([u.mean()-u.min(),u.max()-u.mean()])
    axs[0,1].barh(range(len(models)),vals,color=[COLORS[c['id']] for c in models],height=.6)
    axs[0,1].errorbar(vals,range(len(models)),xerr=np.array(errs).T,fmt='none',color='#444444',lw=.6,capsize=2)
    axs[0,1].set(title='固定参数的0.25 h与1 h差异',xlabel='输出差异RMSE（m³/s）',yticks=range(len(models)),yticklabels=[c['label'] for c in models])
    complexities=[]
    for c in models:
        total=30+(0 if c['type'] in ['base','jointbase'] else 1 if c['type']=='empirical' else c['mode'] if c['type']=='physical' else 22 if c['type']=='neural' else 1+len(__import__('sympy').sympify(c['equation_logit']).atoms(__import__('sympy').Float)))
        complexities.append(dict(model=c['id'],label=c['label'],core_parameters=30,increment_parameters=total-30,total_parameters=total,
            expression_complexity=c.get('complexity',0)))
    frame(OUT/'tables/model_complexity.csv',complexities)
    axs[1,0].barh(range(len(models)),[c['total_parameters'] for c in complexities],color=[COLORS[c['id']] for c in models]);axs[1,0].set(yticks=range(len(models)),yticklabels=[c['label'] for c in models],xlabel='最终运行参数数目',title='额外复杂度与运行参数')
    pars=pd.read_csv(OUT/'tables/base_parameters.csv');axs[1,1].plot(pars.Lr_km,pars.tau_h,'o-',ms=3,color='#566F82',label='有效传播时间')
    axs[1,1].set(xlabel='有效演算河长（km）',ylabel='基础传播时间（h）',title='基础参数跨河长变化');grid(axs[1,1])
    for ax,s in zip(axs.flat,'abcd'):panel(ax,s)
    export(fig,'06_complexity_and_robustness','参数、步长与复杂度检查限制精度结论',4)

    fig,axs=plt.subplots(2,1,figsize=(WIDTH,125/25.4));fig.subplots_adjust(left=.11,right=.98,top=.94,bottom=.18,hspace=.55)
    gains=pd.read_csv(OUT/'tables/event_improvements.csv');g=gains[(gains.model==nn)&(gains.subset=='validation')].sort_values('event_id')
    axs[0].bar(range(len(g)),g.NSE_gain,color=['#668A9C' if v>=0 else '#B96D69' for v in g.NSE_gain]);axs[0].axhline(0,color='#555555',lw=.6)
    axs[0].set(xticks=range(len(g)),xticklabels=g.event_id,ylabel='相对RWEC-0的NSE增量',title='21场开发内验证：改善分布与保留的退化');axs[0].tick_params(axis='x',labelrotation=50);panel(axs[0],'a')
    j=case('E094');N=d['ns'][j];t=np.arange(N)-48
    for key in compact:axs[1].plot(t,pred[key][9,j,:N]-d['obs'][j,:N],color=COLORS[key],ls=STYLES[key],label=names[key])
    axs[1].axhline(0,color='#555555',lw=.5);axs[1].set(title='E094：附加项可能扩大局部误差',xlabel='距洪水起点时间（h）',ylabel='预测减实际流量（m³/s）');grid(axs[1]);panel(axs[1],'b')
    legend(fig,axs[1],4)
    export(fig,'07_preserved_failure','总体改善不代表每场洪水、每段过程都改善',2)
    save(OUT/'provenance/figure_contracts.json',audit)

if __name__=='__main__':main()
