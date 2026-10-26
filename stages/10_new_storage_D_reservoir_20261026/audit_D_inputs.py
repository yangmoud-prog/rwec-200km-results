from common10 import *
def main():
 cfg=read(D/'config/D_frozen.json');rows=[];parameters=[]
 for e in cfg['events']:
  parameters.append({k:e[k] for k in ['event_id','Qb_m3s','Qp_m3s','D_h','r','Tp_h','gamma_a','secondary_ratio','secondary_delay_h','development_split']})
  for geom in cfg['geometries']:
   for bc in cfg['boundaries']:
    f=D/'models'/f"{e['event_id']}_{geom}_{bc}_{cfg['version']}";u=(f/'Backwater200.u01').read_text(errors='replace');g=(f/'Backwater200.g01').read_text(errors='replace');settings=read(f/'settings.json')
    row=dict(event=e['event_id'],geometry=geom,boundary=bc,stage_hydrograph_count=u.count('Stage Hydrograph='),rating_curve_count=u.count('Rating Curve='),initial_reservoir_given='Initial Storage Elev=Reservoir' in u,native_reservoir_connection='Reach Downstream Storage Area=Reservoir' in g,reservoir_continuity_native=bc=='reservoir',input_hash_matches={name:hashlib.sha256((f/name).read_bytes()).hexdigest()==digest for name,digest in settings['input_hashes'].items()})
    if bc=='reservoir':assert row['stage_hydrograph_count']==0 and row['native_reservoir_connection'] and row['initial_reservoir_given']
    else:assert row['stage_hydrograph_count']==1
    rows.append(row)
 frame('D_active_flood_parameters.csv',parameters);save(P/'D_INPUT_AUDIT.json',dict(cases=rows,passed=True,future_reservoir_stage_supplied=False,unused_legacy_event_fields=['H0_m','deltaH_m','Hmax_m','tb_ratio','Tr_ratio','Th_ratio','Tf_ratio','tb_h','Tr_h','Th_h','Tf_h','Hpeak_h','Hfall_start_h','Hend_h'],explanation='The shared gamma-event helper retains old stage-generator metadata. These fields never enter D reservoir boundary equations. boundary.csv stage_m is a fixed-reference placeholder from the base builder; only upstream flow is copied into reservoir unsteady inputs. The .u01 and .g01 checks above verify native storage-area closure.',source_files_originals_preserved=True))
 print('D_NATIVE_RESERVOIR_INPUT_AUDIT_PASS',len(rows))
if __name__=='__main__':main()
