"""Native HEC level-pool reservoir plus a free overflow connection; all inputs causal."""
from hec200 import *

def section(x,complex_channel):
 if not complex_channel:
  z=.000075*x;w=500.;bank=12.;flood=0.;nc=.025;nf=.025
  pts=[(0,z+24),(36,z+12),(72,z),(572,z),(608,z+12),(644,z+24)];banks=(36,608)
 else:
  # Smooth variation has the same mean slope and end elevations as the uniform control.
  z=.000075*x+(.000025*200000/(2*np.pi))*np.sin(2*np.pi*x/200000)
  w=500+150*np.sin(2*np.pi*x/200000);bank=3.5;flood=600.;nc=.025;nf=.05
  left=72+flood+3*bank;right=left+w
  pts=[(0,z+27.5),(72,z+bank),(72+flood,z+bank),(left,z),(right,z),(right+3*bank,z+bank),(right+3*bank+flood,z+bank),(right+3*bank+flood+72,z+27.5)]
  banks=(72+flood,right+3*bank)
 pts=[(round(a,2),round(b,4)) for a,b in pts];banks=tuple(round(v,2) for v in banks)
 return dict(x=x,z=round(z,4),width=w,bank=bank,n_channel=nc,n_floodplain=nf,points=pts,banks=banks)

def storage_area(name,area_m2,minimum=0.):
 # In a metric HEC .g file elevation-volume tabulation uses 1000 m^3.
 hs=np.arange(minimum,21,1.);pairs=[(h,area_m2*(h-minimum)/1000) for h in hs]
 g=[f'Storage Area={name:16s},-1000,0','Storage Area Surface Line= 5 ']
 for x,y in [(-1000,-1000),(-2000,-1000),(-2000,1000),(-1000,1000),(-1000,-1000)]:g.append(f'{x:16g}{y:16g}                ')
 g+=['Storage Area Type= 1 ','Storage Area Area=','Storage Area Min Elev=',f'Storage Area Vol Elev= {len(hs)} ']
 for i in range(0,len(pairs),5):g.append(''.join(f'{h:8g}{v:8g}' for h,v in pairs[i:i+5]))
 g+=['Storage Area Is2D=0','Storage Area Point Generation Data=,,,','Storage Area 2D Points= 0 ','Storage Area Mannings=0.06','']
 return g

def build_D(case,e,geom='complex',boundary_kind='reservoir',dx=2000,dt=60,end_h=552,area=60000000.,crest=3.5):
 folder=ROOT/'models'/case
 if (folder/'D_design.json').exists():return folder
 folder=build(case,e,'ref',dx=dx,dt=dt,end_h=end_h)
 assert not (folder/'run.json').exists()
 g=[f'Geom Title={case}','Program Version=6.31','Viewing Rectangle=-5000,200000,5000,-5000','','River Reach=Ideal200        ,Main            ','Reach XY= 2 ',f'{0:16g}{0:16g}{200000:16g}{0:16g}','Reverse River Text=0']
 if boundary_kind=='reservoir':g+=['Reach Downstream Storage Area=Reservoir       ']
 g+=[''];sections=[]
 for x in np.arange(200000,-1,-dx):
  s=section(float(x),geom=='complex');sections.append(s);length=dx if x>0 else 0;pts=s['points']
  g+=[f'Type RM Length L Ch R = 1 ,{x},{length},{length},{length}',f'#Sta/Elev= {len(pts)} ']
  for i in range(0,len(pts),5):g.append(''.join(f'{a:8.2f}{b:8.4f}' for a,b in pts[i:i+5]))
  g+=['#Mann= 3 ,0,0',''.join(f'{a:8g}{n:8g}{0:8g}' for a,n in [(0,s['n_floodplain']),(s['banks'][0],s['n_channel']),(s['banks'][1],s['n_floodplain'])]),f"Bank Sta={s['banks'][0]},{s['banks'][1]}",'XS Rating Curve=0,0',f"XS HTab Starting El and Incr={s['z']+.01:.4f},.025,1001",'XS HTab Horizontal Distribution=5,5,5','Exp/Cntr=0.3,0.1','']
 length=1000/(1.7*(5-crest)**1.5)
 if boundary_kind=='reservoir':
  g+=storage_area('Reservoir',area)+storage_area('TailReceiver',1e10,-1.)
  g+=['Connection=ReleaseWeir     ,-1500,0','Connection Desc=Fixed passive broad crested free overflow','Connection Line=2',f'{-1500:16g}{0:16g}{-2500:16g}{0:16g}','Connection Up SA=Reservoir       ','Connection Dn SA=TailReceiver    ','Conn Routing Type= 3 ','Conn Use RC Family=False','Conn OverFlow Method 2D=False','Conn Weir WD=10','Conn Weir Coef=1.7','Conn Weir Is Ogee= 0 ','Conn Simple Spill Pos Coef=0.05','Conn Simple Spill Neg Coef=0.05','Conn Weir SE= 2 ',f'{0:8g}{crest:8g}{length:8.3f}{crest:8g}','Conn HTab HWMax=20','Conn HTab TWMax=20','']
 g+=['Use User Specified Reach Order=0','GIS Units=METERS',''];text(folder/f'{NAME}.g01','\n'.join(g))
 if boundary_kind=='reservoir':
  u=(folder/f'{NAME}.u01').read_text();start=u.index('Boundary Location=Ideal200        ,Main            ,0,');end=u.index('Met Point Raster Parameters=',start)
  u=u[:start]+u[end:];u=u.replace('Use Restart=0','Use Restart=0\nInitial Storage Elev=Reservoir       ,5\nInitial Storage Elev=TailReceiver    ,0')
  text(folder/f'{NAME}.u01',u)
 cfg=json.loads((folder/'settings.json').read_text(encoding='utf-8'));cfg.update(design_id='compound_reservoir_D_v1',kind=boundary_kind,geometry_kind=geom,htab_dh_m=.025,htab_rows=1001)
 cfg['input_hashes']={f'{NAME}{s}':sha(folder/f'{NAME}{s}') for s in ['.g01','.p01','.prj','.u01']};dump(folder/'settings.json',cfg)
 dump(folder/'D_design.json',dict(geometry=geom,boundary=boundary_kind,sections=sections,reservoir_area_m2=area,storage_curve='V_R=area*H_R; SI volume table scaled by 1000',release_C=1.7,weir_length_m=round(length,3),weir_crest_m=crest,release_rule='Qrelease=1.7*L*max(H_R-crest,0)^1.5; free overflow',tail_receiver_area_m2=1e10,tail_initial_H=0,tail_must_remain_below_crest=True,dx_m=dx,dt_s=dt,end_h=end_h))
 return folder

if __name__=='__main__':
 e=event(qp=5000,duration=132,event_id='DP01')
 f=build_D('DP01_complex_reservoir_v2',e,end_h=552);print(f,flush=True);print(run(f),flush=True)
