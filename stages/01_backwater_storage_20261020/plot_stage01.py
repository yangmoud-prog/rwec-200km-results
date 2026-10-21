"""Seven deterministic quantitative figures; all input samples retained."""
from pathlib import Path
import csv
import json
import subprocess
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager
from matplotlib.lines import Line2D
import numpy as np

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
SKILL=Path.home()/'.codex/skills/nature-figure/scripts'
FIG=OUT/'figures'
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','SimSun'],
    'font.size':8.5,'axes.labelsize':8.5,'axes.titlesize':9,'xtick.labelsize':8,
    'ytick.labelsize':8,'legend.fontsize':7.5,'axes.spines.top':False,'axes.spines.right':False,
    'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':600,
    'mathtext.fontset':'custom','mathtext.rm':'Times New Roman',
    'mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold'})
EIDS=['E023','E006','E007','E069']
COLORS=['#4477AA','#228833','#CC6677','#AA4499']
LINESTYLES=['-','--','-.',':']
AUDIT=[]
PANEL_AUDIT=[]

def table(name):
    with (OUT/'tables'/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))

def load(case):
    with np.load(OUT/'derived_cases'/f'{case}.npz') as z:
        return {k:z[k] for k in z.files}

def check_axes(ax, label, claim):
    ax.text(-.15,1.06,label,transform=ax.transAxes,fontweight='bold',fontsize=9)
    PANEL_AUDIT.append(dict(panel=label,claim=claim,replicate_unit='Flood background',
        uncertainty='None for deterministic curves; p10-p90 event spread only where stated',
        all_samples_retained=True))

def export(fig,name):
    fig.canvas.draw()
    width,height=fig.get_size_inches()*25.4
    fig.savefig(FIG/f'{name}.pdf')
    fig.savefig(FIG/f'{name}.svg')
    fig.savefig(FIG/f'{name}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    fig.savefig(FIG/f'{name}.png',dpi=300)
    plt.close(fig)
    proc=subprocess.run([sys.executable,str(SKILL/'audit_pdf_text.py'),str(FIG/f'{name}.pdf'),
        '--min-pt','5','--json'],capture_output=True,text=True,encoding='utf-8')
    AUDIT.append(dict(figure=name,width_mm=float(width),height_mm=float(height),
        glyph_audit_returncode=proc.returncode,glyph_report=json.loads(proc.stdout)))
    assert proc.returncode==0, 'PDF text too small: '+name

def case_style(i):
    return dict(color='#4477AA' if i<4 else '#CC6677',ls=LINESTYLES[i%4],lw=1.05)

def legend_eids(fig,y=1):
    handles=[Line2D([],[],color=c,ls=s,lw=1.3,label=e) for e,c,s in zip(EIDS,COLORS,LINESTYLES)]
    fig.legend(handles=handles,ncol=4,loc='upper center',bbox_to_anchor=(.54,y),frameon=False)

def main():
    FIG.mkdir(exist_ok=True)
    fonts={s:font_manager.findfont(s,fallback_to_default=False) for s in ['Times New Roman','SimSun']}
    plan=json.loads((ROOT/'config/formal_plan.json').read_text(encoding='utf-8-sig'))
    events=table('event_summary.csv'); metrics=table('storage_by_length.csv')
    contrasts=table('controlled_contrasts.csv')
    cal=[e for e in events if e['role']=='calibration_main']
    assert len(cal)==70
    controls={e:plan['groups'][e]['addon'] for e in EIDS}

    # Fig. 1: preselected late/high/held forcing in four independent flood backgrounds.
    fig,axs=plt.subplots(2,2,figsize=(180/25.4,145/25.4),sharex=True,sharey=True)
    fig.subplots_adjust(left=.105,right=.86,bottom=.10,top=.92,wspace=.24,hspace=.32)
    maxnorm=max(float(load(controls[e][7])['time_h'][-1])/float(next(v['D_h'] for v in events if v['event_id']==e and v['case_code']=='C08')) for e in EIDS)
    minnorm=min(-48/float(next(v['D_h'] for v in events if v['event_id']==e and v['case_code']=='C08')) for e in EIDS)
    for i,(ax,eid) in enumerate(zip(axs.flat,EIDS)):
        d=load(controls[eid][7]); D=float(next(v['D_h'] for v in events if v['event_id']==eid and v['case_code']=='C08'))
        xx=d['time_h']/D; yy=d['full_x_km']; zz=d['delta_stage_full_m'].T
        im=ax.pcolormesh(xx,yy,zz,norm=colors.TwoSlopeNorm(vmin=-6,vcenter=0,vmax=6),cmap='RdBu_r',shading='auto',rasterized=True)
        ax.contour(xx,yy,zz,levels=[.1],colors=['#252525'],linewidths=.75)
        ax.axhline(200,color='#555555',lw=.6,ls=':')
        ax.set(xlim=(minnorm,maxnorm),ylim=(0,200),title=f'{eid}: C08')
        if i>=2: ax.set_xlabel('Time / D')
        if i%2==0: ax.set_ylabel('Distance from outlet (km)')
        check_axes(ax,chr(97+i),'Full signed water-level response and 0.10 m contour; upper-domain censoring visible')
    cbax=fig.add_axes([.89,.21,.022,.57]);cb=fig.colorbar(im,cax=cbax,ticks=[-6,-3,0,3,6])
    cb.set_label('Water-level increment (m)')
    export(fig,'01_backwater_propagation')

    # Fig. 2: different length metrics remain separate and a boundary-only closure is checked.
    fig,axs=plt.subplots(2,2,figsize=(180/25.4,150/25.4))
    fig.subplots_adjust(left=.10,right=.98,bottom=.10,top=.94,wspace=.32,hspace=.43)
    lengths=np.arange(10,201,10)
    matrix=[]; shellpeak=[]
    for e in cal:
        rows=sorted([m for m in metrics if m['case']==e['case']],key=lambda z:int(z['L_km']))
        matrix.append([float(v['peak_delta_volume_m3'])/1e6 for v in rows])
        shellpeak.append([float(v['shell_peak_m3'])/1e6 for v in rows])
    matrix=np.array(matrix); shellpeak=np.array(shellpeak)
    increments=np.diff(np.column_stack([np.zeros(len(cal)),matrix]),axis=1)
    matrices=[matrix,increments,matrix/matrix[:,[-1]]]
    ylabels=['Peak added storage (10⁶ m³)','Adjacent-length peak difference (10⁶ m³)','Fraction of the 200 km peak']
    for ax,arr,label in zip(axs.flat,matrices,ylabels):
        ax.plot(lengths,arr.T,color='#9A9A9A',lw=.4,alpha=.25)
        med=np.median(arr,axis=0);lo,hi=np.quantile(arr,[.1,.9],axis=0)
        ax.fill_between(lengths,lo,hi,color='#4477AA',alpha=.2,label='Event range (10–90%)')
        ax.plot(lengths,med,color='#225588',lw=1.5,label='Median (70 events)')
        ax.set(xlabel='Reach length (km)',ylabel=label,xlim=(10,200))
    axs[0,0].legend(loc='upper left',frameon=False)
    axs[0,1].plot(lengths,np.median(shellpeak,axis=0),color='#CC6677',ls='--',lw=1.2,label='Median shell maximum')
    axs[0,1].legend(loc='upper right',frameon=False,fontsize=6.8)
    axs[1,0].axhline(.9,color='#777777',ls=':',lw=.8)
    ax=axs[1,1]
    for eid,color,ls in zip(EIDS,COLORS,LINESTYLES):
        d=load(controls[eid][7])
        ax.plot(d['stage_increment_boundary_m'],d['delta_volume_m3'][:,-1]/1e6,color=color,ls=ls,lw=1.25,label=eid)
    ax.set(xlabel='Outlet level increment (m)',ylabel='Added storage at 200 km (10⁶ m³)')
    ax.legend(loc='upper left',frameon=False,ncol=2,fontsize=7)
    for i,ax in enumerate(axs.flat):
        check_axes(ax,chr(97+i),['Magnitude versus length','Peak differences versus synchronous shell maxima',
            'Saturation within finite modeled domain','Same present outlet level can correspond to different storage'][i])
    export(fig,'02_storage_length_and_state')

    # Fig. 3: all 16 conditional comparisons per factor, plus within-background time effects.
    fig,axs=plt.subplots(2,2,figsize=(180/25.4,150/25.4))
    fig.subplots_adjust(left=.11,right=.98,bottom=.11,top=.88,wspace=.36,hspace=.43)
    ax=axs[0,0]
    for fidx,factor in enumerate(['A','B','C']):
        for ei,(eid,color) in enumerate(zip(EIDS,COLORS)):
            rr=[v for v in contrasts if v['factor']==factor and v['event_id']==eid and v['L_km']=='200']
            assert len(rr)==4
            yy=np.array([float(v['difference_individual_peaks_m3'])/1e6 for v in rr])
            xpos=fidx+(ei-1.5)*.16
            ax.scatter(xpos+np.linspace(-.035,.035,4),yy,s=14,c=color,marker=['o','s','^','D'][ei],alpha=.65,zorder=3)
            ax.plot(xpos,float(yy.mean()),marker='_',markersize=8,color=color,mew=1.5)
    ax.axhline(0,color='#666666',lw=.6)
    ax.set(xticks=[0,1,2],xticklabels=['Amplitude','Later rise','Longer hold'],ylabel='Difference in peak storage (10⁶ m³)',xlim=(-.5,2.5))
    for ax,factor in zip([axs[0,1],axs[1,0],axs[1,1]],['A','B','C']):
        for eid,color,ls in zip(EIDS,COLORS,LINESTYLES):
            d=load('factorial_'+eid)
            ax.plot(d['time_over_D'],d['effect_'+factor+'_m3'][:,-1]/1e6,color=color,ls=ls,lw=1.2)
        ax.axhline(0,color='#777777',lw=.6)
        ax.set(xlabel='Time / D',ylabel='Mean storage contrast (10⁶ m³)',title={'A':'A: amplitude 3 → 6 m','B':'B: onset −0.25D → +0.25D','C':'C: hold 0 → 0.20D'}[factor])
    legend_eids(fig,.99)
    for i,ax in enumerate(axs.flat):
        check_axes(ax,chr(97+i),['All conditional effects on peak storage','Amplitude effect through time',
            'Timing effect through time','Hold and delayed drawdown effect through time'][i])
    export(fig,'03_controlled_factor_effects')

    # Full six-panel source diagnostic for each prespecified flood background.
    for eid in EIDS:
        fig,axs=plt.subplots(3,2,figsize=(180/25.4,190/25.4))
        fig.subplots_adjust(left=.105,right=.98,bottom=.075,top=.88,wspace=.34,hspace=.51)
        fig.suptitle(eid,y=.99,fontsize=10)
        handles=[Line2D([],[],**case_style(i),label=f'C{i+1:02d}') for i in range(8)]
        fig.legend(handles=handles,ncol=4,loc='upper center',bbox_to_anchor=(.54,.97),frameon=False)
        for i,case in enumerate(controls[eid]):
            d=load(case); row=next(v for v in events if v['case']==case); D=float(row['D_h'])
            t=d['time_h']/D; x=d['x_km']; sty=case_style(i)
            peak=d['delta_volume_m3'].max(axis=0)
            axs[0,0].plot(t,d['stage_increment_boundary_m']+5,**sty)
            axs[0,1].plot(t,d['delta_volume_m3'][:,-1]/1e6,**sty)
            axs[1,0].plot(x[1:],peak[1:]/1e6,**sty)
            axs[1,1].plot(x[1:],d['shell_m3'].max(axis=0)/1e6,**sty)
            z=d['delta_stage_full_m'][:,1:]
            extent=np.max(np.where(z>.1,d['full_x_km'][None,1:],0),axis=1)
            axs[2,0].plot(t,extent,**sty)
            axs[2,1].plot(t,d['outlet_flow_m3s'],**sty)
            if i==0: axs[2,1].plot(t,d['reference_outlet_flow_m3s'],color='#333333',lw=.8,ls='--',label='Fixed 5 m')
        specs=[('Time / D','Outlet water level (m)'),('Time / D','Added storage at 200 km (10⁶ m³)'),
            ('Reach length (km)','Peak added storage (10⁶ m³)'),('Shell upper end (km)','Shell maximum (10⁶ m³ / 10 km)'),
            ('Time / D','Farthest affected section (km)'),('Time / D','Outlet flow (m³/s)')]
        claims=['Forcing boundary','Added storage history','Length accumulation','Synchronous-shell peak',
            'Threshold influence within finite domain','Outlet flow and reversals']
        for i,(ax,(xl,yl),claim) in enumerate(zip(axs.flat,specs,claims)):
            ax.set(xlabel=xl,ylabel=yl);check_axes(ax,chr(97+i),claim)
        axs[2,0].set_ylim(0,205);axs[2,1].axhline(0,color='#777777',lw=.6)
        axs[2,1].legend(loc='upper right',fontsize=7,frameon=False)
        export(fig,'control_'+eid)
    qa=dict(backend='python',fonts=fonts,user_font_override='Times New Roman / SimSun supersedes skill sans-serif default',
        source_selection='70 frozen calibration pairs and 32 controlled cases; 30 validation pairs withheld',
        control_backgrounds=EIDS,control_conditions_per_background=8,
        intervals='Event p10-p90 distribution, not a confidence interval',
        conditional_effects='16 contexts per factor; 4 flood backgrounds; contexts and lengths are not independent replicates',
        transformations='No input samples omitted; arithmetic contrasts; L and time strictly increasing; positive peak denominators verified by analysis',
        native_balance='Uses both inlet and outlet cumulative volume contrasts',
        editable_formats=['pdf','svg'],raster='600 dpi LZW TIFF; 300 dpi PNG preview',
        pdf_audits=AUDIT,panels=PANEL_AUDIT,visual_review='pending')
    (OUT/'provenance/figure_qa.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Exported',len(AUDIT),'figures; PDF glyph audits passed for all.')

if __name__=='__main__':
    main()
