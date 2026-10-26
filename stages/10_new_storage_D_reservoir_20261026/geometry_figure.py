from figures10 import *
import sys
sys.path.insert(0,str(D/'scripts'))
from build_D import section
def main():
 fig,axs=plt.subplots(2,2,figsize=(10.6,6.1),layout='constrained');x=np.linspace(0,200000,401)
 axs[0,0].plot(x/1000,.000075*x,color='#8193A7',label='均匀河槽')
 axs[0,0].plot(x/1000,.000075*x+(.000025*200000/(2*np.pi))*np.sin(2*np.pi*x/200000),color='#217F9A',label='复杂河槽')
 axs[0,0].set_ylabel('河床高程（m）');axs[0,0].set_title('a  坡度连续变化：0.05—0.10‰');axs[0,0].legend()
 axs[0,1].plot(x/1000,500+150*np.sin(2*np.pi*x/200000),color='#217F9A');axs[0,1].axhline(500,color='#8193A7',ls='--');axs[0,1].set_ylabel('主槽底宽（m）');axs[0,1].set_title('b  底宽连续变化：350—650 m')
 for ax in axs[0]:ax.set_xticks([0,50,100,150,200],['0/200','50/150','100/100','150/50','200/0']);ax.set_xlabel('距出口 / 距入口（km）')
 for pos,style in [(50000,'-'),(150000,'--')]:
  s=section(pos,True);pts=np.array(s['points']);centre=.5*(pts[0,0]+pts[-1,0]);axs[1,0].plot(pts[:,0]-centre,pts[:,1]-s['z'],ls=style,label=f'距出口 {pos/1000:g} km / 距入口 {200-pos/1000:g} km')
 axs[1,0].set_ylim(0,8);axs[1,0].set_xlabel('距断面中心（m）');axs[1,0].set_ylabel('距主槽底高程（m）');axs[1,0].set_title('c  主槽 n=0.025；两侧滩地 n=0.050');axs[1,0].legend(fontsize=7)
 hh=np.linspace(3,9,100);bw=read(D/'models/D001_complex_reservoir_v3/D_design.json')['weir_length_m'];axs[1,1].plot(hh,1.7*bw*(hh-2)**1.5,color='#BC5A45',label='已知自由溢流规则');axs[1,1].scatter([5],[1000],s=18,c='#111111');axs[1,1].set_xlabel('库水位（m）');axs[1,1].set_ylabel('预定泄量（m³/s）');axs[1,1].set_title('d  V=6×10⁷H；堰顶 2 m，堰宽 %.3f m'%bw);axs[1,1].legend(fontsize=7)
 for ax in axs.flat:ax.grid(alpha=.15)
 export(fig,'14_D_geometry_and_reservoir')
if __name__=='__main__':main()
