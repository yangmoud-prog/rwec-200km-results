from common10 import *
def main():
 r=pd.read_csv(P/'tables/D_pilot_numerical_quality.csv');p=pd.read_csv(P/'tables/D_pilot_pairs.csv');rows=[]
 for version in ['v2','v3']:
  rr=r[(r.version==version)&(r.boundary=='reservoir')];pp=p[(p.version==version)&(p.geometry=='complex')];fixed=r[r.boundary=='fixed']
  eligible=len(rr)==8 and len(pp)==4 and bool(rr.numerical_gate.all() and fixed.numerical_gate.all()) and pp.endpoint_range_km.max()>20 and pp.max_suppression_20_80_m3s.max()>100 and pp.max_recovery_20_80_m3s.max()>100
  score=float(np.mean(np.minimum(pp.max_suppression_20_80_m3s,pp.max_recovery_20_80_m3s)))
  rows.append(dict(version=version,eligible=bool(eligible),physical_strength=score))
 selected=max([q for q in rows if q['eligible']],key=lambda q:q['physical_strength'])['version'];assert selected==read(D/'config/D_frozen.json')['version']
 save(P/'D_SELECTION_REPLAY.json',dict(replayed_from_corrected_extraction=rows,selected=selected,original_frozen_record=read(D/'config/physical_selection.json'),selected_design_unchanged=True,geometry_or_boundary_changed=False,NSE_used=False,interpretation='The original v2 eligibility flag was false while v3 was selected. Complete re-extraction makes v2 eligible too; the predeclared physical-strength criterion still selects v3. Both records are retained.'))
 print(read(P/'D_SELECTION_REPLAY.json'))
if __name__=='__main__':main()
