from common import *
import h5py
d=load(OUT/'data/A_E006.npz')
for hour in [-48,0,24,72,120]:
 i=np.flatnonzero(GRID==hour)[0];print(hour,'zvarfirstlast',d['z'][i,[0,-1]],'zref',d['zref'][i,[0,-1]],'dz',d['delta_stage'][i,[0,-1]])
plan=read(A/'config/formal_plan.json')
for c in plan['formal_cases']:
 if c['event_id']=='E006' and c['role']=='main':
  print(c['kind'],c['case']);raw=load(A/'results'/c['case']/'reach_timeseries.npz');print('raw',raw['x_from_outlet_km'][[0,-1]],'H0',raw['stage_m'][[0,48,72,120],0]);print('rawt',raw['time_since_flood_h'][[0,48,72,120]])
  with h5py.File(A/'models'/c['case']/'Backwater200.p01.hdf') as f:
   a=f['Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Cross Sections/Water Surface'][:];print('hdfshape',a.shape,'rawlast',a[[0,48,72,120],-1]);print('hdf time',f['Results/Unsteady/Output/Output Blocks/Base Output/Unsteady Time Series/Time'][:][[0,48,72,120]])
