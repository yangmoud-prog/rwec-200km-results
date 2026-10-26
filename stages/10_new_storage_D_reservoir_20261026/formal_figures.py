from figures10 import *
DOMAINS={'DUF':'均匀河槽＋固定水位','DUR':'均匀河槽＋动态库区','DCF':'复杂河槽＋固定水位','DCR':'复杂河槽＋动态库区'}
def main():
 df=pd.read_csv(P/'tables/D_event_metrics.csv');v=df.query("subset=='validation'");means=v.groupby(['domain','model','x_from_outlet_km'])[['NSE','RMSE','KGE']].mean()
 fig,axs=plt.subplots(4,3,figsize=(11.5,10.8),layout='constrained')
 for i,domain in enumerate(DOMAINS):
  for j,metric in enumerate(['NSE','RMSE','KGE']):
   ax=axs[i,j]
   for name in NAMES:
    d=means.loc[(domain,name)];ax.plot(200-d.index,d[metric],label=NAMES[name],color=COLORS[name],ls=STYLES[name],marker='o',markersize=3)
   ax.set_xticks([80,100,120,140,160,180,200],['120/80','100/100','80/120','60/140','40/160','20/180','0/200'],rotation=30);ax.set_xlabel('距出口 / 距入口（km）');ax.set_ylabel(metric+('（m³/s）' if metric=='RMSE' else ''));ax.set_title(f'{chr(97+i*3+j)}  {DOMAINS[domain]}');ax.grid(alpha=.15)
 axs[0,0].legend(loc='best',fontsize=7);export(fig,'08_D_four_conditions_mean_skill')
 cfg=read(D/'config/D_frozen.json')
 for domain,geom in [('DUR','uniform'),('DCR','complex')]:
  fig,axs=plt.subplots(2,3,figsize=(11.5,6.1),layout='constrained')
  for i,e in enumerate(['D005','D008']):
   d=load(D/'models'/f"{e}_{geom}_reservoir_{cfg['version']}"/'extracted.npz');mask=d['time']>=0
   for j,x in enumerate([120,20,0]):
    ax=axs[i,j];idx=int(np.argmin(abs(d['x']-x)));ax.plot(d['time'][mask]/24,d['q'][mask,idx],color='#111111',lw=1.4,label='HEC-RAS')
    for name in NAMES:
     pp=load(P/'predictions'/f'{domain}_{e}_{name}.npz');ax.plot(pp['time']/24,pp['q'][int((200-x)/20)],color=COLORS[name],ls=STYLES[name],lw=1.0,label=NAMES[name])
    ax.set_title(f'{e}：距出口 {x} km / 距入口 {200-x} km');ax.set_xlabel('预报起点后时间（天）');ax.set_ylabel('流量（m³/s）');ax.grid(alpha=.15)
  axs[0,0].legend(fontsize=7);export(fig,f'09_{domain}_full_hydrographs')
 paired=pd.read_csv(P/'tables/D_paired_improvements.csv').query("subset=='validation'");fig,axs=plt.subplots(2,2,figsize=(9.5,6.4),layout='constrained')
 for ax,domain in zip(axs.flat,DOMAINS):
  d=paired.query('domain==@domain');events=list(d.event.unique());yy=np.arange(len(events))
  for j,name in enumerate(['new','general','NLM']):
   vals=d.groupby('event')[name+'_minus_original'].mean().reindex(events);ax.plot(vals,yy+(j-1)*.17,'o',color=COLORS[name],label=NAMES[name],markersize=4)
  ax.axvline(0,color='.5',lw=.7);ax.set_yticks(yy,events);ax.set_title(DOMAINS[domain]);ax.set_xlabel('每场七断面平均 NSE 改善');ax.grid(axis='x',alpha=.15)
 axs[0,0].legend(fontsize=7);export(fig,'10_D_event_paired_gains')
 factorial=pd.read_csv(P/'tables/D_factorial_optimized_error.csv').query("subset=='validation' and x_from_outlet_km==20");physical=pd.read_csv(P/'tables/D_factorial_physical_discharge.csv').query("subset=='validation' and x_from_outlet_km==20");components=['geometry_at_fixed','backwater_uniform','backwater_complex','interaction'];labels=['固定边界下\n几何变化','均匀河槽\n动态回水','复杂河槽\n动态回水','几何与回水\n交互项'];fig,axs=plt.subplots(1,2,figsize=(10.5,3.6),layout='constrained')
 vals=physical.groupby('component').RMS_m3s.mean().reindex(components);axs[0].bar(np.arange(4),vals,color=['#8193A7','#94B8C5','#217F9A','#BC5A45']);axs[0].set_xticks(np.arange(4),labels);axs[0].set_ylabel('流量差的均方根（m³/s）');axs[0].set_title('水动力效应：同一固定断面的配对差')
 for j,name in enumerate(NAMES):
  vals=factorial.query('model==@name')[components].mean();axs[1].plot(np.arange(4),vals,color=COLORS[name],marker='o',ls=STYLES[name],label=NAMES[name])
 axs[1].axhline(0,color='.5',lw=.7);axs[1].set_xticks(np.arange(4),labels);axs[1].set_ylabel('分别率定后的 RMSE 差（m³/s）');axs[1].set_title('模型误差：距出口 20 km / 距入口 180 km');axs[1].legend(fontsize=7);export(fig,'11_D_geometry_backwater_attribution')
 phase=pd.read_csv(P/'tables/D_storage_release_phase_errors.csv').query("subset=='validation' and domain=='DCR'");fig,axs=plt.subplots(1,2,figsize=(9.4,3.6),layout='constrained')
 for ax,x in zip(axs,[20,0]):
  for name in NAMES:
   vals=phase.query('x_from_outlet_km==@x and model==@name').groupby('phase').RMSE_m3s.mean().reindex(['storage','release','weak_or_recovered']);ax.plot(range(3),vals,c=COLORS[name],ls=STYLES[name],marker='o',label=NAMES[name])
  ax.set_xticks(range(3),['蓄入阶段','释放阶段','弱变化／恢复']);ax.set_title(f'距出口 {x} km / 距入口 {200-x} km');ax.set_ylabel('平均 RMSE（m³/s）');ax.grid(alpha=.15)
 axs[0].legend(fontsize=7);export(fig,'12_D_storage_release_phase_errors')
 sens=pd.read_csv(P/'tables/D_fixed_parameter_dt_sensitivity.csv');fig,axs=plt.subplots(2,2,figsize=(9.5,6.0),layout='constrained')
 for ax,domain in zip(axs.flat,DOMAINS):
  for name in NAMES:
   d=sens.query('domain==@domain and model==@name').groupby('dt_factor').NSE.mean();base=v.query('domain==@domain and model==@name').NSE.mean();ax.plot([.5,1,2],[d.loc[.5],base,d.loc[2]],color=COLORS[name],ls=STYLES[name],marker='o',label=NAMES[name])
  ax.set_xticks([.5,1,2]);ax.set_xlabel('各子段计算时段的统一倍数');ax.set_ylabel('验证洪水平均 NSE');ax.set_title(DOMAINS[domain]);ax.grid(alpha=.15)
 axs[0,0].legend(fontsize=7);export(fig,'13_D_fixed_parameter_step_sensitivity')
if __name__=='__main__':main()
