from figures10 import *
def main():
 fig,axs=plt.subplots(2,3,figsize=(11.5,6.0),layout='constrained');event='D008';domain='DCR'
 for i,target in enumerate([20,0]):
  jj=int(np.flatnonzero(TARGETS==target)[0]);reach=int((200-target)/20)-1
  for name in NAMES:
   d=load(P/'predictions'/f'storage_{domain}_{event}_{name}.npz');t=d['time']/24;vol=np.cumsum(d['reach_volume_bookkeeping'],axis=0)[reach];anom=vol-vol[0];hec=d['HEC_geometric_storage'][jj];hec=hec-hec[0];ht=d['HEC_time']/24
   if name=='original':
    axs[i,0].plot(ht,hec/1e6,c='#111111',lw=1.5,label='HEC 几何体积');axs[i,1].plot(ht[1:],np.diff(hec)/3600,c='#111111',lw=1.5,label='HEC 几何体积')
   axs[i,0].plot(t,anom/1e6,c=COLORS[name],ls=STYLES[name],label=NAMES[name]);axs[i,1].plot(ht[1:],np.diff(anom[::4])/3600,c=COLORS[name],ls=STYLES[name],label=NAMES[name]);proj=np.cumsum(d['inter_reach_resampling_volume'],axis=0)[reach];axs[i,2].plot(t,proj/1000,c=COLORS[name],ls=STYLES[name],label=NAMES[name])
  for j,label in enumerate(['有效储量异常与几何储量对照','完整过程储量变化率','时段衔接的数值体积偏差']):
   axs[i,j].set_title(f'{label}\n距出口 {target} km / 距入口 {200-target} km',fontsize=9);axs[i,j].set_xlabel('预报起点后时间（天）');axs[i,j].grid(alpha=.15)
  axs[i,0].set_ylabel('储量异常（百万 m³）');axs[i,1].set_ylabel('储量变化率（m³/s）');axs[i,2].set_ylabel('累计数值体积（千 m³）')
 axs[0,0].legend(fontsize=7);export(fig,'15_D_effective_storage_and_numerical_interfaces')
if __name__=='__main__':main()
