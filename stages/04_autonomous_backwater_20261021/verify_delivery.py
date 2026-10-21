from common import *
import subprocess,sys,re,datetime

def run():
    v=read(OUT/'provenance/verification.json');manifest=read(B/'config/pilot_input_manifest.json')
    audits=[]
    for case in manifest['cases']:
        folder=B/'models'/case['case'];rewritten=[]
        for fn,digest in case['input_hashes'].items():
            actual=sha(folder/fn)
            if actual!=digest:
                assert fn=='Backwater200.prj'
                text=(folder/fn).read_text();required={'Current Plan':'p01','Geom File':'g01','Unsteady File':'u01','Plan File':'p01'}
                for k,val in required.items():assert f'{k}={val}' in text
                assert 'SI Units' in text
                rewritten.append(dict(file=fn,pre_compute_sha256=digest,post_HEC_save_sha256=actual,
                    interpretation='HEC Project_Save added axes/DSS/GIS display metadata; hydraulic geometry, plan and boundary files unchanged'))
        quality=read(B/'results'/case['case']/'quality.json')
        assert quality['basic_quality_pass'] and quality['tail_recovered'] and quality['initialization_pass']
        audits.append(dict(case=case['case'],hydraulic_g_p_u_inputs_unchanged=True,HEC_project_metadata_changes=rewritten,
            native_balance_max_fraction=float(max(quality['native_max_balance_residual_fraction_storage_by_length'][1:])),all_pilot_basic_gates_passed=True))
    v['B_input_audits']=audits;v['B_HEC_time_space_sensitivity_complete']=False
    save(OUT/'provenance/verification.json',v)
    # Basic formula-source checks catch escaped control characters and missing delimiters.
    for name in ['EQUATIONS.md','REPORT.md']:
        text=(OUT/name).read_text(encoding='utf-8')
        assert not any(ord(c)<32 and c not in '\n\r\t' for c in text)
        assert '\t' not in text,'Unexpected tab in LaTeX source'
        assert text.count('$$')%2==0
    # Use the project-wide font override; retain source-validator warning evidence.
    skill=Path(r'C:\Users\mryan\.codex\skills\nature-figure')
    cmd=[sys.executable,str(skill/'scripts/validate_figure.py'),str(OUT/'figures.py'),'--json']
    result=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8')
    source=__import__('json').loads(result.stdout)
    assert [f['check_id'] for f in source['findings'] if f['level']=='FAIL']==['FONT-FAMILY']
    qa=dict(source_validator=source,font_exception='User explicitly requires Times New Roman and SimSun instead of skill sans-serif default',
        width_warning='Static parser reads 180/25.4 as180 inches; actual PDF MediaBox independently checked180mm',
        Chinese_capitalization_warning='English capitalization is inapplicable to Chinese label',
        final_size_visual_inspection='Every static panel and complete figure checked; deterministic hydrographs have no random confidence intervals',
        backend='Python only',source_data='all70 A development pairs plus32 controls; all4 B pilot backgrounds; no event exclusions',
        animation_requested=False,delivery_format='static figures, tables and editable LaTeX')
    glyph=[]
    for p in sorted((OUT/'figures').glob('*.pdf')):
        result=subprocess.run([sys.executable,str(skill/'scripts/audit_pdf_text.py'),str(p),'--min-pt','5','--json'],capture_output=True,text=True,encoding='utf-8')
        assert result.returncode==0,result.stdout
        glyph.append(__import__('json').loads(result.stdout))
    qa['PDF_glyph_audits']=glyph
    save(OUT/'provenance/figure_qa.json',qa)
    print('Scientific, source, boundary-freeze, and PDF glyph checks complete',flush=True)
if __name__=='__main__':run()
