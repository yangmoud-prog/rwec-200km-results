"""Development comparisons with paired metrics and complete flow processes."""
from pathlib import Path
import json,subprocess,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from campaign import OUT,load_np,read,save

FIG=OUT/'figures';W=180/25.4
SKILL=Path.home()/'.codex/skills/nature-figure/scripts'
COLORS=['#777777','#4477AA','#CC6677'];STYLES=['--','-','-.']
NAMES=['RWEC-0','RWEC-1','RWEC-2'];AUDITS=[];PANELS=[]
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','SimSun'],
 'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'legend.fontsize':7.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
 'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':600,
 'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold'})

def label(ax,c,claim):
    ax.text(-.16,1.06,c,transform=ax.transAxes,fontweight='bold',fontsize=9)
    PANELS.append(dict(panel=c,claim=claim,visual_review='pending'))

def legend(fig,y=.99,truth=False):
    hh=[Line2D([],[],color=c,ls=s,lw=1.4,label=n) for n,c,s in zip(NAMES,COLORS,STYLES)]
    if truth:hh.insert(0,Line2D([],[],color='#111111',lw=1.4,label='HEC-RAS'))
    fig.legend(handles=hh,ncol=len(hh),loc='upper center',bbox_to_anchor=(.55,y),frameon=False)

def export(fig,name):
    fig.canvas.draw()
    fig.savefig(FIG/f'{name}.pdf');fig.savefig(FIG/f'{name}.svg')
    fig.savefig(FIG/f'{name}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    fig.savefig(FIG/f'{name}.png',dpi=300)
    for p in PANELS:
        if 'figure' not in p:p['figure']=name
    plt.close(fig)
    p=subprocess.run([sys.executable,str(SKILL/'audit_pdf_text.py'),str(FIG/f'{name}.pdf'),'--min-pt','5','--json'],
        capture_output=True,text=True,encoding='utf-8')
    AUDITS.append(dict(figure=name,returncode=p.returncode,report=json.loads(p.stdout)))
    assert p.returncode==0

def model_lines(ax):
    for m,(c,s) in enumerate(zip(COLORS,STYLES)):yield m,c,s

def main():
    fonts={s:font_manager.findfont(s,fallback_to_default=False) for s in ['Times New Roman','SimSun']}
    df=pd.read_csv(OUT/'tables/event_length_metrics.csv');dev=df[df.dataset=='development']
    imp=pd.read_csv(OUT/'tables/improvements.csv');imp=imp[imp.dataset=='development']
    phase=pd.read_csv(OUT/'tables/phase_event_metrics.csv');phase=phase[(phase.dataset=='development')&(phase.L_km==200)]
    fig,axs=plt.subplots(2,3,figsize=(W,150/25.4))
    fig.subplots_adjust(left=.085,right=.985,bottom=.10,top=.88,wspace=.44,hspace=.43)
    for ax,metric in zip(axs[0],['NSE','KGE','RMSE']):
        for m,c,s in model_lines(ax):
            z=dev[dev.model==NAMES[m]].pivot(index='event_id',columns='L_km',values=metric)
            x=z.columns.to_numpy();a=z.to_numpy();lo,hi=np.quantile(a,[.1,.9],axis=0)
            ax.fill_between(x,lo,hi,color=c,alpha=.12)
            ax.plot(x,np.median(a,axis=0),color=c,ls=s,lw=1.2)
        ax.set(xlabel='Reach length (km)',ylabel=metric+(' (m³/s)' if metric=='RMSE' else ''))
    for m in [1,2]:
        z=imp[imp.model==NAMES[m]].pivot(index='event_id',columns='L_km',values='delta_NSE')
        x=z.columns.to_numpy();a=z.to_numpy();lo,hi=np.quantile(a,[.1,.9],axis=0)
        axs[1,0].fill_between(x,lo,hi,color=COLORS[m],alpha=.15)
        axs[1,0].plot(x,np.median(a,axis=0),color=COLORS[m],ls=STYLES[m],lw=1.2)
        worse=np.mean(a< -1e-6,axis=0)
        axs[1,1].plot(x,worse,color=COLORS[m],ls=STYLES[m],lw=1.2)
    axs[1,0].axhline(0,color='#555555',lw=.6)
    axs[1,0].set(xlabel='Reach length (km)',ylabel='Paired NSE improvement')
    axs[1,1].set(xlabel='Reach length (km)',ylabel='Fraction worse than RWEC-0',ylim=(-.02,1.02))
    for m,c,s in model_lines(axs[1,2]):
        vals=[];lows=[];highs=[]
        for ph in ['stage_rise','stage_hold','stage_drawdown','after_stage_return']:
            a=phase[(phase.model==NAMES[m])&(phase.phase==ph)].RMSE.to_numpy()
            vals.append(np.median(a));lows.append(np.quantile(a,.25));highs.append(np.quantile(a,.75))
        vals=np.array(vals)
        axs[1,2].errorbar(np.arange(4)+(m-1)*.12,vals,yerr=[vals-np.array(lows),np.array(highs)-vals],
            color=c,ls=s,lw=.9,marker='o',ms=2,capsize=2)
    axs[1,2].set(xticks=range(4),xticklabels=['Rise','Hold','Fall','After'],ylabel='200 km phase RMSE (m³/s)')
    for i,ax in enumerate(axs.flat):label(ax,chr(97+i),['NSE distribution','KGE distribution','RMSE distribution','Within-event improvement',
        'Retained degradation cases','Boundary-phase errors at fixed200km'][i])
    legend(fig);export(fig,'01_accuracy_by_length')

    meta=read(OUT/'data/metadata.json');d=load_np(OUT/'data/development.npz');pred=load_np(OUT/'predictions/development.npz')['flow']
    chosen=read(OUT/'provenance/representative_events.json');eids=chosen['preset']+[chosen['worst_additional']]
    fig,axs=plt.subplots(2,2,figsize=(W,145/25.4))
    fig.subplots_adjust(left=.105,right=.98,bottom=.10,top=.88,wspace=.30,hspace=.42)
    for i,(ax,eid) in enumerate(zip(axs.flat,eids)):
        j=next(j for j,r in enumerate(meta['development']) if r['event_id']==eid);e=meta['development'][j]['event'];N=int(d['ns'][j]);t=np.arange(N)-48
        display_L=200 if i<3 else int(imp[(imp.event_id==eid)&(imp.model=='RWEC-2')].sort_values('delta_NSE').iloc[0].L_km)
        li=display_L//10-1
        top=e['tb_h']+e['Tr_h'];fall=top+e['Th_h']
        for start,end,color in [(e['tb_h'],top,'#EAF1F6'),(top,fall,'#F6EEE7'),(fall,e['Hend_h'],'#EEF2EB')]:
            ax.axvspan(start,end,color=color,zorder=0)
        ax.plot(t,d['obs'][j,:N],color='#111111',lw=1.25)
        for m,c,s in model_lines(ax):ax.plot(t,pred[m,li,j,:N],color=c,ls=s,lw=1.1)
        ax.axhline(0,color='#777777',lw=.5)
        ax.set(xlabel='Time since flood onset (h)',ylabel='Outlet flow (m³/s)',title=f'{eid}: {display_L} km'+(' (degradation)' if i==3 else ''))
        label(ax,chr(97+i),'Complete outlet hydrograph, including reversals and poor performance')
    legend(fig,truth=True);export(fig,'02_complete_hydrographs')

    for eid in ['E023','E006','E007','E069']:
        z=load_np(OUT/'predictions'/f'control_response_{eid}.npz');t=z['time_h'];a=z['differences']
        fig,axs=plt.subplots(1,3,figsize=(W,85/25.4))
        fig.subplots_adjust(left=.09,right=.985,bottom=.22,top=.71,wspace=.42)
        fig.suptitle(eid,y=.99,fontsize=10);legend(fig,y=.93,truth=True)
        for fi,ax in enumerate(axs):
            ax.plot(t,a[fi,3,0,19],color='#111111',lw=1.2)
            for m,c,s in model_lines(ax):ax.plot(t,a[fi,3,m+1,19],color=c,ls=s,lw=1.1)
            ax.axhline(0,color='#888888',lw=.5)
            ax.set(xlabel='Time since flood onset (h)',ylabel='Outlet-flow difference (m³/s)',title=['A: C08 − C04','B: C08 − C06','C: C08 − C07'][fi])
            label(ax,chr(97+fi),'Fixed-parameter within-background boundary-response difference')
        export(fig,'control_response_'+eid)

    p=pd.read_csv(OUT/'tables/parameters.csv');cond=pd.read_csv(OUT/'tables/conditional_alpha_by_length.csv')
    gr=pd.read_csv(OUT/'tables/conditional_shared_parameters.csv')
    fig,axs=plt.subplots(2,3,figsize=(W,145/25.4));fig.subplots_adjust(left=.09,right=.985,bottom=.10,top=.88,wspace=.48,hspace=.45)
    for ax,field in zip(axs[0],['k','n','tau_h']):
        for m,c,s in model_lines(ax):
            z=p[p.model==NAMES[m]];ax.plot(z.L_km,z[field],color=c,ls=s,lw=1.2)
        ax.set(xlabel='Reach length (km)',ylabel={'k':'k','n':'n','tau_h':'Propagation time (h)'}[field])
        if field in ['k','n']:
            assert np.all(p[field]>0);ax.set_yscale('log')
    for m in [1,2]:
        z=cond[cond.model==NAMES[m]];assert np.all(z.primary_shared_alpha>0)
        axs[1,0].plot(z.L_km,z.conditional_alpha/z.primary_shared_alpha,color=COLORS[m],ls=STYLES[m],lw=1.2)
    axs[1,0].axhline(1,color='#777777',ls=':',lw=.7)
    axs[1,0].set(xlabel='Reach length (km)',ylabel='Conditional α / shared α')
    order=['event_half_even','event_half_odd','length_10_60km','length_70_130km','length_140_200km']
    for m in [1,2]:
        zz=gr[gr.model==NAMES[m]].set_index('group').loc[order]
        shared=float(p[p.model==NAMES[m]].alpha.iloc[0]);assert shared>0
        axs[1,1].plot(range(5),zz.alpha/shared,color=COLORS[m],ls=STYLES[m],marker='o',ms=3,lw=1)
    axs[1,1].axhline(1,color='#777777',ls=':',lw=.7)
    axs[1,1].set(xticks=range(5),xticklabels=['Even','Odd','Short','Mid','Long'],ylabel='Conditional α / shared α')
    zz=gr[gr.model=='RWEC-2'].set_index('group').loc[order]
    axs[1,2].plot(range(5),zz.theta,color=COLORS[2],marker='o',ms=3,lw=1)
    theta=float(p[p.model=='RWEC-2'].theta.iloc[0]);axs[1,2].axhline(theta,color='#777777',ls=':',lw=.7)
    axs[1,2].set(xticks=range(5),xticklabels=['Even','Odd','Short','Mid','Long'],ylabel='200 km / characteristic length')
    for i,ax in enumerate(axs.flat):label(ax,chr(97+i),['Original parameter k','Complement parameter n','Lag parameter',
        'Conditional length sharing diagnostic','Conditional event and length groups','Conditional saturation scale'][i])
    legend(fig);export(fig,'03_parameters_and_sharing')

    z=load_np(OUT/'predictions/timestep_example.npz');t=z['time_h']
    j=next(j for j,r in enumerate(meta['development']) if r['event_id']=='E007');N=int(d['ns'][j])
    fig,axs=plt.subplots(1,3,figsize=(W,88/25.4),sharey=True)
    fig.subplots_adjust(left=.09,right=.985,bottom=.20,top=.75,wspace=.24)
    for m,ax in enumerate(axs):
        ax.plot(t,d['obs'][j,:N],color='#111111',lw=1,ls=':',label='HEC-RAS')
        for fac,col,sty in zip([1,2,4],['#777777','#4477AA','#CC6677'],['--','-','-.']):
            ax.plot(t,z[f'model{m}_factor{fac}'],color=col,ls=sty,lw=1.1,label=f'{1/fac:g} h')
        ax.set(xlabel='Time since flood onset (h)',title=NAMES[m]);label(ax,chr(97+m),'Fixed-parameter timestep sensitivity with hourly aggregation')
    axs[0].set_ylabel('Outlet flow (m³/s)')
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,ncol=4,loc='upper center',bbox_to_anchor=(.55,.98),frameon=False)
    export(fig,'04_timestep_sensitivity')

    profile=pd.read_csv(OUT/'tables/saturation_product_profile.csv').sort_values('ell_m')
    summary=read(OUT/'summary.json');a1=summary['models'][1]['alpha'];a2=summary['models'][2]['alpha'];ell=summary['models'][2]['ell_m']
    fig,axs=plt.subplots(1,2,figsize=(W,88/25.4));fig.subplots_adjust(left=.11,right=.98,bottom=.20,top=.85,wspace=.32)
    xx=np.arange(10,201,10)*1000
    axs[0].plot(xx/1000,a1*xx/1000,color=COLORS[1],lw=1.3,label='RWEC-1')
    axs[0].plot(xx/1000,a2*ell*(-np.expm1(-xx/ell))/1000,color=COLORS[2],ls='-.',lw=1.3,label='RWEC-2')
    axs[0].set(xlabel='Reach length (km)',ylabel='Effective coefficient length (km)');axs[0].legend(frameon=False)
    axs[1].plot(profile.ell_m/1000,profile.delta_NSE_vs_primary*1e6,color=COLORS[2],marker='o',ms=3,lw=1.1)
    axs[1].axhline(0,color='#777777',lw=.6)
    axs[1].set(xlabel='Characteristic length (km)',ylabel='NSE difference (10⁻⁶)')
    label(axs[0],'a','Actual fitted length dependence of potential correction')
    label(axs[1],'b','Fixed-product conditional sensitivity; not a confidence interval')
    export(fig,'05_saturation_identifiability')
    save(OUT/'provenance/figure_qa.json',dict(font_files=fonts,backend='python',pdf_audits=AUDITS,panels=PANELS,
        source_data='Complete70development events;32controls; no test targets.',
        samples_omitted=False,uncertainty='10–90% event distributions in metric panels; phase RMSE median and25–75% event distribution; not confidence intervals.',
        width_mm=180,editable_text=True,TIFF_dpi=600,PNG_preview_dpi=300,
        user_font_priority='Times New Roman and SimSun override generic skill font defaults',visual_review='pending'))
    print('Rendered',len(AUDITS),'figures; all PDF glyph audits passed.',flush=True)

if __name__=='__main__':main()
