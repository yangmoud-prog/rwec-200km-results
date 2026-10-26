from common10 import *
import h5py,re
def main():
 cfg=read(D/'config/D_frozen.json');rows=[]
 for e in cfg['events']:
  for geom in cfg['geometries']:
   for bc in cfg['boundaries']:
    folder=D/'models'/f"{e['event_id']}_{geom}_{bc}_{cfg['version']}"
    with h5py.File(folder/'Backwater200.p01.hdf') as f:
     base=f['Results/Unsteady/Output/Output Blocks/Base Output'];summary=base['Summary Output/Cross Sections'];iterations=summary['Maximum Iteration'][:];z=base['Unsteady Time Series/Cross Sections/Water Surface'][:];q=base['Unsteady Time Series/Cross Sections/Flow'][:];info=f['Geometry/Cross Sections/Property Tables/XSEC Info'][:];vals=f['Geometry/Cross Sections/Property Tables/XSEC Value'][:];maxfr=0.;miny=1e6
     for j,(st,n,_) in enumerate(info):
      table=vals[st:st+n];A=np.interp(z[:,j],table[:,0],table[:,1:4].sum(axis=1));T=np.interp(z[:,j],table[:,0],table[:,16]);maxfr=max(maxfr,float(np.sqrt((q[:,j]**2*T/(9.81*A**3))).max()));miny=min(miny,float(z[:,j].min()-table[0,0]))
     saiter=int(base['Summary Output/Storage Areas/Maximum Iteration'][:].max()) if bc=='reservoir' else 0
    logfile=(folder/'Backwater200.p01.computeMsgs.txt').read_text(errors='replace');warnings=[line.strip() for line in logfile.splitlines() if re.search(r'error|warning|unstable|failed',line,re.I) and 'Volume Accounting Error' not in line]
    row=dict(case=folder.name,event=e['event_id'],geometry=geom,boundary=bc,maximum_XS_iteration=int(iterations.max()),maximum_SA_iteration=saiter,maximum_Froude_output=maxfr,minimum_depth_above_HTab_base_m=miny,log_warning_lines=len(warnings),warning_excerpt=' | '.join(warnings[:4]),gate=bool(iterations.max()<30 and saiter<30 and miny>0 and maxfr<1))
    rows.append(row)
 frame('D_independent_numerical_audit.csv',rows);save(P/'D_NUMERICAL_AUDIT.json',dict(cases=len(rows),all_passed=all(r['gate'] for r in rows),maximum_iterations=max(r['maximum_XS_iteration'] for r in rows),maximum_Froude=max(r['maximum_Froude_output'] for r in rows),minimum_depth=min(r['minimum_depth_above_HTab_base_m'] for r in rows)))
 print(read(P/'D_NUMERICAL_AUDIT.json'))
if __name__=='__main__':main()
