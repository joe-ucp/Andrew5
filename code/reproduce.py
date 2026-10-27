"""Public reproduction entry point; scientific inputs in, external outputs out.

check: verify public scientific files without data or writes.
demo: run a small synthetic diagnostic, not research data.
influence: reproduce four retained-table arithmetic summaries.
stages: list the ordered scientific computation and prerequisites.
stage NAME/all: execute the same code tree against explicit external inputs.
"""
from pathlib import Path
import argparse
import datetime
import importlib
import importlib.metadata
import json
import os
import platform
import sys
import time

sys.dont_write_bytecode=True
os.environ['PYTHONDONTWRITEBYTECODE']='1'
from data_paths import resolve_data_root
from public_runtime import ROOT, Execution, registry, scientific_file, sha, parameters_sha256

# Each stage has declared input groups, compact prerequisites, owned outputs and
# dependencies. A single-stage run may reuse registered retained prerequisites;
# all executes the chain in order and retains newly generated prerequisites.
STAGES={
 'reference':('reference_analysis','main',[],['reference','primary'],[],
              'calibration_results calibration_gate', ['original_results/reference_feature_selection.tsv','original_results/reference_cell_scores.tsv','original_results/calibration_simulations.tsv.gz'],[]),
 'calibration':('calibration_diagnostics','main',['reference'],['reference'],[],
                'calibration_operating_points reference_feature_capture reference_capture_summary',[],[]),
 'primary':('primary_analysis','main',['reference'],['primary'],['evidence/calibration_status.json'],
            'primary_results countermodel_results sample_score_summary sign_flip_reference',['original_results/primary_cell_scores.tsv.gz'],[]),
 'challenge':('biological_challenge','main',['reference'],['challenge'],[],
              'biological_challenge_results biological_challenge_strata',['original_results/biological_challenge_cell_scores.tsv.gz'],['data_primary/GSE217837']),
 'identifiability':('verify_identifiability','construct',['primary'],[],[],
                    'identifiability_witness.csv',[],[]),
 'depth':('depth_stability','main',['primary'],['primary'],[],
          'depth_stability_summary',['original_results/depth_stability_replicates.tsv','original_results/depth_stability_sample_means.tsv'],[]),
 'countermodels':('countermodel_stability','main',['primary'],['primary'],['evidence/calibration_status.json','evidence/primary_results.tsv','evidence/countermodel_results.tsv'],
                  'countermodel_stability_summary countermodel_stability_donor_checks',['original_results/countermodel_stability_pools.tsv.gz','original_results/countermodel_stability_group_draws.tsv','original_results/countermodel_stability_batches.tsv'],[]),
 'operating-points':('revision_operating_points','main',['reference','primary'],[],[],
                     'operating_points roc_frontiers observed_spearman observed_depth_summary',[],['revision_01/operating_points']),
 'mechanism':('revision_mechanism','main',['primary'],['primary'],[],
              'population_audit observed_components observed_component_contrasts observed_depth_relationships observed_depth_bins thinning_sample_summaries thinning_donor_contrasts thinning_stability_summary',[],[]),
 'benchmark':('revision_benchmark','run_primary',['reference','primary','mechanism'],['reference','primary'],[],
              'benchmark_feature_stability benchmark_gene_membership benchmark_reference benchmark_primary_samples benchmark_primary_contrasts',[],['revision_01/benchmark']),
 'independent':('revision_independent','main',['benchmark'],['baron'],['evidence/revision_01/benchmark_programs.json','evidence/revision_01/baron_input_manifest.json'],
                'independent_native_components independent_state_benchmark independent_feature_coverage independent_missing_universe independent_thinning independent_component_thinning',[],['revision_01/independent']),
 'detection':('revision_detection_expectation','main',['mechanism','independent'],['primary','baron'],['evidence/revision_01/thinning_sample_summaries.tsv','evidence/revision_01/independent_component_thinning.tsv','evidence/revision_01/baron_input_manifest.json'],
              'detection_expectation_samples',[],['revision_01/detection_expectation']),
 'nuisance':('revision_nuisance','main',['primary','countermodels'],['primary'],[],
             'nuisance_precision_batches nuisance_precision_summary nuisance_pool_drift nuisance_coupling_differences nuisance_population_audit',[],['revision_01/nuisance']),
 'transport':('revision02_decisions','transport',['independent'],['reference','baron'],['evidence/revision_01/baron_input_manifest.json','evidence/revision_01/independent_thinning.tsv'],
              'source_operating_points target_native_thresholds threshold_summary',[],['revision_02/threshold_transport']),
 'primary-diagnostics':('revision02_decisions','primary',['mechanism','benchmark'],['primary'],['evidence/revision_01/benchmark_programs.json'],
                        'zero_baseline_samples zero_baseline_contrasts ucell_floor_samples ucell_floor_contrasts fixed_cohort_summary',[],['revision_02/primary_diagnostics']),
 'genes':('revision02_gene_specificity','main',['mechanism'],['primary'],[],
          'gene_sample_summaries gene_donor_contrasts gene_contribution_summary gene_concentration matched_candidate_support matched_baseline_quality matched_program_comparison',[],['revision_02/gene_specificity']),
 'influence':('revision03_influence','main',['genes','primary-diagnostics','mechanism'],[],[],
              'influence_single influence_topk influence_selected influence_summary',[],[]),
 'figures':('revision_figures','main',['calibration','countermodels','operating-points','mechanism','benchmark','independent','nuisance','transport','primary-diagnostics','genes'],[],[], '',[],[]),
}


def evidence_folder(stage):
    if stage in list(STAGES)[:7]:return 'evidence'
    if stage in {'transport','primary-diagnostics','genes'}:return 'evidence/revision_02'
    if stage=='influence':return 'evidence/revision_03'
    return 'evidence/revision_01'


def table_path(stage,name):return evidence_folder(stage)+'/'+(name if name.endswith('.csv') else name+'.tsv')


def source_group(relative):
    if relative.startswith('data_reference/'):return 'reference'
    if relative.startswith('revision_01/sources/'):return 'baron'
    if relative.startswith('data_primary/GSE217837/'):return 'challenge'
    return 'primary'


def environment():
    packages={}
    for line in (ROOT/'code/requirements.txt').read_text().splitlines():
        if '==' not in line:continue
        name,expected=line.strip().split('==')
        try:actual=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:actual=None
        packages[name]={'expected':expected,'installed':actual,'matches':actual==expected}
    return {'python':platform.python_version(),'platform':platform.platform(),'packages':packages}


def check():
    records=registry()
    for relative in records['compact']:scientific_file(relative)
    result={'status':'PASS','scope':'Registered compact scientific inputs only; no data loaded and no files written.',
            'compact_inputs_verified':len(records['compact']),'parameters_sha256':parameters_sha256(),'environment':environment()}
    if (ROOT/'PUBLIC_FILES.json').is_file():
        from public_inventory import verify_manifest
        result['public_files']=verify_manifest()
        result['scope']='Published public-file inventory and registered compact scientific inputs; no data loaded and no files written.'
    if any(v['installed'] is None for v in result['environment']['packages'].values()):
        raise RuntimeError('Install dependencies: python -m pip install -r code/requirements.txt')
    return result


def fetch_sources(data_root,groups):
    from fetch_sources import download
    from public_runtime import _verify
    source=resolve_data_root(data_root,for_write=True)
    for relative,record in registry()['data'].items():
        if record['kind']!='source' or source_group(relative) not in groups:continue
        path=source/relative
        if not path.exists():download({**record,'local_path':relative},path,source)
        _verify(path,record,relative)


def seed_inputs(context,stage,strict_chain):
    requested=['evidence/frozen_score.json']+list(STAGES[stage][4])
    if stage=='figures':requested += list(context.records['compact'])
    # Stage-local EV reads refer to the external tree. A stand-alone stage may
    # reuse verified earlier summaries; a complete chain must use its own work.
    if strict_chain:
        for dependency in STAGES[stage][2]:
            receipt=context.output/'receipts'/(dependency+'.json')
            if not receipt.is_file() or json.loads(receipt.read_text()).get('status')!='PASS':
                raise FileNotFoundError('Complete chain lacks completed prerequisite stage: '+dependency)
        for previous in list(STAGES)[:list(STAGES).index(stage)]:
            receipt=context.output/'receipts'/(previous+'.json')
            if not receipt.is_file():continue
            context.required_generated.update(table_path(previous,name) for name in STAGES[previous][5].split())
            context.required_generated.update('data/'+name for name in STAGES[previous][6])
            context.required_data_prefixes+=tuple(STAGES[previous][7])
            if previous=='benchmark':context.required_generated.add('evidence/revision_01/benchmark_programs.json')
            if previous=='reference':context.required_generated.add('evidence/calibration_status.json')
            if previous=='mechanism':
                for rel in context.records['data']:
                    if rel.startswith('revision_01/GSM'):context.required_generated.add('data/'+rel)
    for relative in dict.fromkeys(requested):
        if relative.startswith('evidence/'):
            context.seed(relative)


def compare_native_log_rank_metrics(retained, generated, score_loader):
    """Check only the native log scorer's rank metrics against certified zero ties.

    Every other value retains the ordinary strict comparison. No score or
    published table is rounded, replaced or refitted by this comparison.
    """
    import numpy as np
    import pandas as pd
    from rank_metric_comparison import compare_rank_metric, log_score_error_radii
    pd.testing.assert_index_equal(retained.index,generated.index)
    pd.testing.assert_index_equal(retained.columns,generated.columns)
    adjusted=generated.copy()
    audits=[]
    columns=['beta_vs_alpha_AUROC','alpha_score_depth_r','beta_score_depth_r']
    for index,row in generated.iterrows():
        old=retained.loc[index]
        if old['method']!='signed_log1p_CPT10K' or row['method']!=old['method']:continue
        for column in columns:
            expected,actual=float(old[column]),float(row[column])
            if np.isclose(expected,actual,rtol=1e-12,atol=2e-14,equal_nan=True):continue
            bundle=score_loader(int(row['donor']),str(row['mode']))
            programs=list(bundle['programs']);methods=list(bundle['methods'])
            pi=programs.index(row['program']);mi=methods.index(row['method'])
            scores=np.asarray(bundle['scores'][pi,mi],dtype=float)
            radii=np.asarray(bundle['log_roundoff_radius'][pi],dtype=float)
            positive=np.asarray(bundle['log_state_is_beta'])
            depths=np.asarray(bundle['total_counts'],dtype=float)
            if positive.dtype!=np.dtype(bool) or positive.shape!=scores.shape or depths.shape!=scores.shape:
                raise ValueError('Native score/state/depth metadata do not align')
            if int(positive.sum())!=int(row['beta_n']) or int((~positive).sum())!=int(row['alpha_n']):
                raise ValueError('Native score counts disagree with the benchmark table')
            # The unchanged programs have at most 100 signed terms and L1
            # magnitude <=2*log1p(10000). This ceiling rejects arbitrary radii.
            ceiling=float(log_score_error_radii(2*np.log1p(10000),100,True))
            if radii.shape!=scores.shape or not np.isfinite(radii).all() or np.any(radii<0) or np.any(radii>ceiling):
                raise ValueError('Invalid exact-cancellation roundoff budget')
            if column=='beta_vs_alpha_AUROC':
                check=compare_rank_metric(expected,actual,scores,radii,metric='auc',positive=positive)
            else:
                selected=positive if column=='beta_score_depth_r' else ~positive
                check=compare_rank_metric(expected,actual,scores[selected],radii[selected],
                                          metric='spearman',covariate=depths[selected])
            if not check['accepted']:
                raise AssertionError('Native log rank metric exceeds its exact-cancellation envelope: '+json.dumps(
                    {'donor':int(row['donor']),'mode':str(row['mode']),'program':str(row['program']),
                     'column':column,'comparison':check}))
            audits.append({'donor':int(row['donor']),'mode':str(row['mode']),
                           'program':str(row['program']),'method':str(row['method']),
                           'column':column,'retained':expected,'generated':actual,**check})
            adjusted.loc[index,column]=expected
    pd.testing.assert_frame_equal(retained,adjusted,check_dtype=False,check_exact=False,rtol=1e-12,atol=2e-14)
    return audits


def compare_outputs(context,stage):
    import pandas as pd
    comparisons=[]
    paths=[context.output/table_path(stage,name) for name in STAGES[stage][5].split()]
    for path in paths:
        if not path.is_file():raise FileNotFoundError('Scientific stage did not produce declared table: '+str(path))
        relative=path.relative_to(context.output).as_posix()
        if relative not in context.records['compact']:raise ValueError('No registered retained comparison: '+relative)
        retained=ROOT/relative
        if sha(retained)!=context.records['compact'][relative]['sha256']:
            raise ValueError('Retained scientific comparison changed: '+relative)
        separator=',' if path.suffix=='.csv' else '\t'
        a=pd.read_csv(retained,sep=separator,float_precision='round_trip')
        b=pd.read_csv(path,sep=separator,float_precision='round_trip')
        ignored=['elapsed_seconds'] if path.name=='population_audit.tsv' else []
        for column in ignored:
            a=a.drop(columns=column,errors='ignore');b=b.drop(columns=column,errors='ignore')
        audits=[]
        if stage=='independent' and path.name=='independent_state_benchmark.tsv':
            import numpy as np
            bundles={}
            def load_scores(donor,mode):
                key=(donor,mode)
                if key not in bundles:
                    if mode not in {'original_universe','shared_features'}:raise ValueError('Unknown native universe')
                    filename=f'data/revision_01/independent/donor{donor}_{mode}_scores.npz'
                    with np.load(context.generated_file(filename),allow_pickle=False) as source:
                        bundles[key]={name:source[name] for name in ['scores','programs','methods',
                            'total_counts','log_roundoff_radius','log_state_is_beta']}
                return bundles[key]
            audits=compare_native_log_rank_metrics(a,b,load_scores)
        else:
            pd.testing.assert_frame_equal(a,b,check_dtype=False,check_exact=False,rtol=1e-12,atol=2e-14)
        comparisons.append({'path':relative,'rows':len(a),'byte_exact':sha(path)==sha(retained),
                            'numeric_rtol':1e-12,'numeric_atol':2e-14,'excluded_columns':ignored,
                            'exact_zero_rank_metric_comparisons':audits})
    return comparisons


def run_stage(stage,data_root,output_dir,fetch=False,strict_chain=False):
    if stage not in STAGES:raise ValueError('Unknown stage: '+stage)
    module_name,entry,deps,groups,prereqs,table_names,names,prefixes=STAGES[stage]
    if fetch:fetch_sources(data_root,groups)
    if stage=='mechanism':
        import csv
        with (ROOT/'primary_sample_manifest.tsv').open(newline='',encoding='utf-8') as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                if row['treatment_code'] in ['D','H_GLP']:
                    names=names+[f"revision_01/{row['accession']}_cell_components.tsv.gz"]+[f"revision_01/{row['accession']}_depth{d}_thinning.npz" for d in [1000,2000,5000,10000]]
    context=Execution(data_root,output_dir,stage,prefixes,names)
    receipt_path=context.output/'receipts'/(stage+'.json')
    if receipt_path.exists():raise FileExistsError('Stage receipt exists; choose a new output directory or a not-yet-run stage: '+stage)
    for name in table_names.split():
        if (context.output/table_path(stage,name)).exists():raise FileExistsError('Stage outputs already exist: '+name)
    if groups and context.source is None:raise ValueError('This stage requires --data-root')
    os.environ['MPLBACKEND']='Agg';os.environ['MPLCONFIGDIR']=str(context.output_path('.cache/matplotlib'))
    if context.source:os.environ['STRIKE007_CACHE']=str(context.source)
    for folder in ['evidence','evidence/revision_01','evidence/revision_02','evidence/revision_03','manuscript/figures','receipts']:
        context.output_path(folder).mkdir(parents=True,exist_ok=True)
    # Import numerical modules before binding path aliases, so a later module
    # cannot inherit another stage's external EV/OUT globals through imports.
    module=importlib.import_module(module_name)
    if stage=='figures':
        for name in ['make_figures','revision02_gene_figure','revision02_figures']:importlib.import_module(name)
    if stage=='challenge':
        sys.path.insert(0,str(ROOT/'scripts_primary'));importlib.import_module('audit_challenge')
    started=time.monotonic()
    with context:
        seed_inputs(context,stage,strict_chain)
        for relative,record in context.records['data'].items():
            if record['kind']=='source' and source_group(relative) in groups:context.verified_data(relative)
        context.bind_modules()
        if stage=='challenge':
            required=['submitted_expression.npz','features.tsv','barcodes.tsv']
            if not all((context.source/'data_primary/GSE217837'/name).is_file() for name in required):
                context.output_path('data/data_primary/GSE217837').mkdir(parents=True,exist_ok=True)
                sys.path.insert(0,str(ROOT/'scripts_primary'))
                prepare=importlib.import_module('audit_challenge');context.bind_modules();prepare.main()
        if stage=='influence':result=module.main(context.output/evidence_folder(stage))
        elif stage=='identifiability':
            args=argparse.Namespace(input=context.evidence_path('primary_cell_scores.tsv.gz'),input_scope='selected_scored_cells',output_dir=context.output/'evidence',delimiter='\t',donor_column=None,treatment_column=None,annotation_column=None,score_column=None,alpha_labels=['alpha','alpha_related','alpha-related','alpha_like','alpha-like'])
            result=module.construct(args)
        elif stage=='figures':
            result=module.main()
            for name in ['make_figures','revision02_gene_figure','revision02_figures']:
                drawing=importlib.import_module(name);context.bind_modules()
                if hasattr(drawing,'main'):drawing.main()
                else:drawing.mechanism();drawing.thresholds();drawing.supplement()
        else:result=getattr(module,entry)()
        comparisons=compare_outputs(context,stage)
        report={'status':'PASS','stage':stage,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'seconds':time.monotonic()-started,'scope':'Same code tree; generated outputs external; retained scientific results not overwritten.',
                'complete_chain_mode':strict_chain,'input_origins':context.reads,'compared_tables':comparisons,
                'parameters_sha256':parameters_sha256(),'code_sha256':sha(ROOT/'code'/f'{module_name}.py')}
        receipt_path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        context.inventory()
    return report


def demo(output_dir):
    import signature_diagnostic as diagnostic
    output=resolve_data_root(output_dir,for_write=True)
    old=sys.argv
    sys.argv=['signature_diagnostic','--example','--data-root',str(output),'--output-prefix','example']
    try:diagnostic.main()
    finally:sys.argv=old
    return {'status':'PASS','scope':'96-cell synthetic example only; no research data or biological validation.','output_dir':str(output)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    sub.add_parser('check');sub.add_parser('stages')
    for command in ['demo','influence']:
        p=sub.add_parser(command);p.add_argument('--output-dir',type=Path,required=True)
    p=sub.add_parser('stage');p.add_argument('name',choices=['all',*STAGES]);p.add_argument('--data-root',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--fetch',action='store_true',help='Download missing registered source files only; immutable checksum mismatches stop')
    args=parser.parse_args()
    if args.command=='check':report=check()
    elif args.command=='stages':report={name:{'depends_on':spec[2],'source_groups':spec[3],'retained_tables_compared':spec[5].split()} for name,spec in STAGES.items()}
    elif args.command=='demo':report=demo(args.output_dir)
    elif args.command=='influence':report=run_stage('influence',None,args.output_dir)
    elif args.name=='all':report=[run_stage(name,args.data_root,args.output_dir,args.fetch,True) for name in STAGES]
    else:report=run_stage(args.name,args.data_root,args.output_dir,args.fetch)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
