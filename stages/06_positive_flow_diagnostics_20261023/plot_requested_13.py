"""Display all thirteen user-requested positive-flow flood hydrographs at20 km."""
from pathlib import Path
import json,hashlib,re,sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager

OUT=Path(__file__).resolve().parent
A=OUT.parents[1]
FIG=OUT/'figures/requested13_x20_20261024';FIG.mkdir(exist_ok=True)
fontManager.addfont(r'C:\Windows\Fonts\times.ttf');fontManager.addfont(r'C:\Windows\Fonts\simsun.ttc')
mpl.rcParams.update({'font.family':['Times New Roman','SimSun'],'font.size':8,
 'axes.labelsize':8,'axes.titlesize':8.5,'xtick.labelsize':7,'ytick.labelsize':7,
 'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.65,'pdf.fonttype':42,'svg.fonttype':'none',
 'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic'})
EXPECTED={'E006':957,'E010':676,'E014':948,'E026':730,'E032':892,'E037':49,
 'E042':661,'E051':90,'E066':1000,'E067':934,'E068':482,'E081':501,'E094':453}
EVENTS=list(EXPECTED)
COLORS={'truth':'#222222','RWEC0':'#686080','linear_joint':'#367E83'}
LABELS={'truth':'HEC-RAS参考','RWEC0':'RWEC-0','linear_joint':'RWEC-1（已有联合版）'}
STYLES={'truth':'-','RWEC0':'--','linear_joint':'-.'}
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def load(p):
 with np.load(p,allow_pickle=False) as z:return {k:z[k] for k in z.files}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
plan=read(A/'config/formal_plan.json');cases={c['event_id']:c for c in plan['formal_cases'] if c['kind']=='var' and c['role']=='main' and c['event_id'] in EVENTS}
data={};rows=[];series=[];manifest=[]
for eid in EVENTS:
 cache=OUT/'data'/f'{eid}.npz';d=load(cache)
 rawfile=A/'results'/cases[eid]['case']/'reach_timeseries.npz';raw=load(rawfile)
 assert np.allclose(d['time'],raw['time_since_flood_h'],atol=1e-6)
 assert np.array_equal(d['truth'][:,0],raw['flow_m3s'][:,2])
 q=d['truth'][:,0];i=int(q.argmin());minimum=float(q[i])
 assert minimum>0 and round(minimum)==EXPECTED[eid]
 rows.append(dict(event=eid,x_from_outlet_km=20,routing_length_km=180,n_hours=len(q),
   time_start_h=float(d['time'][0]),time_end_h=float(d['time'][-1]),minimum_HEC_flow_m3s=minimum,
   minimum_time_h=float(d['time'][i]),rounded_minimum_m3s=round(minimum),screenshot_minimum_m3s=EXPECTED[eid],
   screenshot_matches=True,negative_samples=int((q<0).sum()),HEC_peak_m3s=float(q.max())))
 for k,t in enumerate(d['time']):series.append(dict(event=eid,time_h=float(t),HEC_m3s=float(q[k]),RWEC0_m3s=float(d['RWEC0'][k,0]),RWEC1_joint_m3s=float(d['linear_joint'][k,0])))
 manifest.append(dict(event=eid,raw_sha256=sha(rawfile),cached_prediction_sha256=sha(cache),raw_time_series_matches_cache=True))
 data[eid]=d
pd.DataFrame(rows).to_csv(OUT/'tables/requested_13_x20_minima.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(series).to_csv(OUT/'tables/requested_13_x20_series.csv',index=False,encoding='utf-8-sig')

def draw(ax,eid,small=False):
 d=data[eid];days=d['time']/24
 for name in COLORS:ax.plot(days,d[name][:,0],color=COLORS[name],ls=STYLES[name],lw=1.25 if name=='truth' else 1.05,label=LABELS[name])
 ax.axhline(1000,color='#C8C8C8',lw=.5,ls=':',zorder=0)
 ax.axvline(0,color='#C8C8C8',lw=.5,ls=':',zorder=0)
 ax.set_ylim(bottom=0);ax.set_xlim(float(days[0]),float(days[-1]))
 ax.set_xlabel('距上游洪水起点时间 (天)');ax.set_ylabel('流量 (m³/s)')
 v=next(r for r in rows if r['event']==eid)
 ax.set_title(f'{eid}  最低流量 {v["minimum_HEC_flow_m3s"]:.1f} m³/s',loc='left',pad=7,fontweight='bold')
 if not small:ax.legend(frameon=False,ncol=3,loc='upper right',fontsize=7.5)
 return ax
def export(fig,name,highres=True):
 fig.savefig(FIG/f'{name}.png',dpi=300)
 fig.savefig(FIG/f'{name}.pdf')
 fig.savefig(FIG/f'{name}.svg')
 if highres:fig.savefig(FIG/f'{name}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
 plt.close(fig)

if '--overviews-only' not in sys.argv:
 for eid in EVENTS:
  fig,ax=plt.subplots(figsize=(180/25.4,78/25.4));fig.subplots_adjust(left=.1,right=.98,bottom=.2,top=.82)
  draw(ax,eid);export(fig,f'{eid}_x20_hydrograph')
  print('Exported',eid,flush=True)
for part,group in enumerate([EVENTS[:7],EVENTS[7:]],1):
 nrows=int(np.ceil(len(group)/2));height_mm=nrows*45+18;fig,axs=plt.subplots(nrows,2,figsize=(180/25.4,height_mm/25.4),squeeze=False)
 fig.subplots_adjust(left=.09,right=.98,bottom=.065,top=1-18/height_mm,wspace=.25,hspace=.64)
 for ax,eid in zip(axs.flat,group):draw(ax,eid,small=True)
 for ax in list(axs.flat)[len(group):]:ax.set_visible(False)
 handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,ncol=3,loc='upper center',bbox_to_anchor=(.52,.995),frameon=False)
 export(fig,f'overview_{part:02d}')

text=['# 距下游出口20 km：用户指定13场洪水过程线','',
 '上游输入位于200 km，演算长度180 km。13场次与所给清单逐一对应；逐小时HEC最低流量取整后全部一致，所有保存小时的参考流量均为正。',
 '', '黑色实线为HEC-RAS参考，紫色虚线为RWEC-0，青色点划线为已有RWEC-1联合版。后者使用真实未来下游水位，且沿用此前有限优化预算的参数；本次只提取并绘制已有结果，没有重新率定。',
 '', '横轴0为上游洪水起点，负时间保留洪水前的已保存过程；没有删点、裁剪流量或重启模型状态。各图纵轴按本场流量范围显示，峰值比较应读取刻度。',
 '', '![前7场总览](figures/requested13_x20_20261024/overview_01.png)','',
 '![后6场总览](figures/requested13_x20_20261024/overview_02.png)','',
 '| 场次：点击打开单场图 | 最低流量/(m³/s) | 最低值时刻/h | HEC洪峰/(m³/s) |','|---|---:|---:|---:|']
for r in rows:
 eid=r['event'];text.append(f'| [{eid}](figures/requested13_x20_20261024/{eid}_x20_hydrograph.png) | {r["minimum_HEC_flow_m3s"]:.2f} | {r["minimum_time_h"]:.0f} | {r["HEC_peak_m3s"]:.2f} |')
text+=['','最低值来自1 h保存点，不表示未保存的60 s求解时刻的绝对极值。',
 '', '[最低流量核对表](tables/requested_13_x20_minima.csv) · [逐小时流量数据](tables/requested_13_x20_series.csv)',
 '', '所有单场与总览均提供PNG、PDF、SVG及600 dpi TIFF，位于同名图件目录。原始计算和30场封存测试均未改变。']
(OUT/'REQUESTED_13_HYDROGRAPHS.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
qa=dict(client_date='2026-10-24',requested_events=EVENTS,input_event_count=13,displayed_event_count=13,omitted_events=0,
 scenario='old channel main dynamic-stage cases',x_from_outlet_km=20,routing_length_km=180,
 all_screenshot_minima_match=True,all_saved_HEC_flows_positive=True,samples_before=sum(r['n_hours'] for r in rows),
 samples_after=sum(r['n_hours'] for r in rows),sample_exclusions=0,raw_truth_clipped=False,new_parameter_fits=0,new_HEC_runs=0,sealed_test_results_opened=0,
 figure_contract={'claim':'Show all requested positive hydrographs and their model differences, without selection by accuracy.',
 'archetype':'quantitative grid','backend':'Python','width_mm':180,'fonts':['Times New Roman','SimSun'],
 'uncertainty':'Exact deterministic saved time series; no aggregate confidence intervals are appropriate.'},sources=manifest)
(OUT/'provenance/requested13_x20_qa.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
print('Completed all13 events; minima checked; fifteen figure sets exported',flush=True)
