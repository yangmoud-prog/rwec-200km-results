from common import *
from models import rwec,muskingum
from moving_balance import audit
sys.path.insert(0,str(P.parent/'07_dynamic_endpoint_research_20261024'))
import dynamic_core as old
dt=.5;t=np.arange(0,240+dt,dt);e=hec.event();u=hec.boundary(e,t,'ref')[0];b=60*(1-np.cos(2*np.pi*t/t[-1]));y0=np.full(10,4.);v0=np.full(10,20000*(500*4+3*16))
new=rwec(u,b,dt,5.,1.5,1.,0.,0.,.000075,y0,v0);reference=old.route(u,np.zeros(len(t)),dt,5.,1.5,1.,0.)
assert new[5];qerror=float(np.max(abs(new[0].T-reference[0])));assert qerror<1.e-8
nested=rwec(u,b,dt,5.,1.5,1.,0.,12.,.000075,y0,v0);assert np.max(abs(nested[0]-new[0]))==0
steady=rwec(np.full(len(t),1000.),np.zeros(len(t)),dt,5.,1.5,1.,0.,0.,.000075,y0,v0);assert np.max(abs(steady[0]-1000))<1.e-8
m=muskingum(u,dt,2.,.05,1.);lm=old.muskingum(u,dt,2.,.05);assert m[2];merror=float(np.max(abs(m[0].T-lm)));assert merror<1.e-7
mb=audit(new[0],new[4],b,dt);assert mb['max_balance_m3']<1.e-4
result=dict(original_RWEC_equivalence_max_m3s=qerror,linear_Muskingum_equivalence_max_m3s=merror,
 constant_flow_max_error_m3s=float(np.max(abs(steady[0]-1000))),alpha_zero_nested_exact=True,
 fixed_volume_balance_max_m3=new[7],moving_balance_max_m3=mb['max_balance_m3'],omitting_sweep_RMS_m3s=mb['rms_no_sweep_m3s'],
 forecast_API_inputs=['upstream_flow','predicted_endpoint','fixed_geometry','known_initial_state','fixed_parameters'],future_HEC_H_or_storage_argument=False,
 no_clipping_negative_output=True,invalid_candidate_returns_failure=True)
save(P/'provenance/model_contract_checks.json',result);print(json.dumps(result,indent=2))
