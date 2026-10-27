"""Describe a supplied signed signature; no corrected score or acceptance decision."""
from pathlib import Path
import argparse
import hashlib
import html
import json
import sys

import numpy as np
import pandas as pd
from scipy import sparse

from data_paths import resolve_data_root, checked_data_destination, data_path
from primary_analysis import map_universe
from reference_analysis import auc
from revision_score import components_sparse
from revision02_math import METHODS, score_bundle, zero_centered_components


def identifiers(values, label):
    values=list(values)
    if not values or any(not isinstance(v,str) or not v.strip() for v in values):
        raise ValueError(label+' must contain nonempty string identifiers')
    if len(set(values))!=len(values):raise ValueError(label+' contains duplicate identifiers')
    return values


def clean_json(value):
    if isinstance(value,dict):return {str(k):clean_json(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [clean_json(v) for v in value]
    if isinstance(value,np.generic):return clean_json(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def groups(frame, columns):
    if not columns:return [({},frame.index.to_numpy())]
    answer=[]
    for key,part in frame.groupby(columns,dropna=False,sort=True):
        if not isinstance(key,tuple):key=(key,)
        answer.append((dict(zip(columns,key)),part.index.to_numpy()))
    return answer


def threshold_counts(values, threshold, minimum):
    values=np.asarray(values)
    finite=np.isfinite(values);n=int(finite.sum());v=values[finite]
    return {'evaluable_cells':n,'non_evaluable_cells':int((~finite).sum()),
            'threshold':threshold,'rule':'score > supplied threshold',
            'positive_n':int((v>threshold).sum()),'not_positive_n':int((v<=threshold).sum()),
            'tie_n':int((v==threshold).sum()),'rate_eligible':n>=minimum,
            'positive_fraction':float((v>threshold).mean()) if n>=minimum else None}


def diagnose(counts, genes, cells, signature, metadata=None, *, thresholds=None,
             baseline_treatment=None, comparison_treatment=None,
             state_positive=None, state_negative=None, minimum_group_size=20):
    """Pure diagnostic on aligned identifiers; input matrix may contain extra genes."""
    if not isinstance(signature,dict):raise ValueError('Signature must be a JSON object')
    genes=identifiers(genes,'genes');cells=identifiers(cells,'cells')
    universe=identifiers(signature.get('universe',[]),'signature universe')
    up=identifiers(signature.get('up',[]),'up block');down=identifiers(signature.get('down',[]),'down block')
    if set(up)&set(down):raise ValueError('Up and down blocks must be disjoint')
    if not (set(up)|set(down))<=set(universe):raise ValueError('Every program gene must belong to the explicit universe')
    if isinstance(minimum_group_size,bool) or not isinstance(minimum_group_size,(int,np.integer)) or minimum_group_size<1:
        raise ValueError('minimum_group_size must be a positive integer')
    x=sparse.csc_matrix(counts,copy=True)
    if x.shape!=(len(genes),len(cells)):raise ValueError('Matrix dimensions disagree with ordered gene/cell identifiers')
    if x.data.dtype.kind not in 'iuf' or not np.isfinite(x.data).all() or (x.data<0).any() or (x.data!=np.floor(x.data)).any():
        raise ValueError('Input must contain finite nonnegative integer counts')
    x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
    matrix_totals=np.asarray(x.sum(axis=0),dtype=float).ravel()
    meta=pd.DataFrame({'cell_id':cells}) if metadata is None else metadata.copy()
    if 'cell_id' not in meta:raise ValueError('Metadata requires cell_id')
    identifiers(meta.cell_id.tolist(),'metadata cell_id')
    if set(meta.cell_id)!=set(cells):raise ValueError('Metadata cell IDs must match matrix cell IDs exactly')
    reordered=meta.cell_id.tolist()!=cells
    meta=meta.set_index('cell_id').loc[cells].reset_index()
    if 'total_counts' in meta and pd.api.types.is_bool_dtype(meta.total_counts):
        raise ValueError('Whole-library totals must be integers, not booleans')
    totals=matrix_totals if 'total_counts' not in meta else pd.to_numeric(meta.total_counts,errors='raise').to_numpy(dtype=float)
    if not np.isfinite(totals).all() or (totals<0).any() or (totals!=np.floor(totals)).any() or (totals<matrix_totals).any():
        raise ValueError('Whole-library totals must be finite integers at least as large as matrix totals')
    eligible=totals>0
    if not eligible.any():raise ValueError('No positive-total cells are eligible')
    mapped=map_universe(x[:,eligible],pd.DataFrame({'symbol':genes}),universe)
    lookup={g:i for i,g in enumerate(universe)}
    comp=zero_centered_components(components_sparse(mapped,[lookup[g] for g in up],[lookup[g] for g in down]),len(universe))
    scores=score_bundle(mapped,totals[eligible],universe,{'supplied':{'up':up,'down':down}})[0]
    np.testing.assert_allclose(scores[0],comp['score'],atol=2e-15,rtol=0)
    table=meta.loc[eligible].reset_index(drop=True)
    table['total_counts']=totals[eligible]
    # Reserved output names never overwrite user metadata silently.
    reserved=set(METHODS)|set(comp)|{'up_detection_fraction','down_detection_fraction','ucell_zero_floor_fraction'}
    if (set(meta.columns)-{'cell_id','total_counts'})&reserved:
        raise ValueError('Metadata columns collide with diagnostic output names')
    for key,value in comp.items():table[key]=value
    table['up_detection_fraction']=comp['k_up']/len(up)
    table['down_detection_fraction']=comp['k_down']/len(down)
    for method,value in zip(METHODS,scores):table[method]=value
    table['ucell_zero_floor_fraction']=np.where(np.isfinite(scores[2]),(scores[2]==0).astype(float),np.nan)
    thresholds={} if thresholds is None else thresholds
    if not isinstance(thresholds,dict):raise ValueError('Thresholds must be a method-to-number object')
    if set(thresholds)-set(METHODS):raise ValueError('Threshold keys must name declared scoring methods')
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not np.isfinite(v) for v in thresholds.values()):
        raise ValueError('Supplied thresholds must be finite numbers')
    columns=[c for c in ['group','donor','treatment','state'] if c in table]
    for column in columns:table[column]=table[column].fillna('<missing>').astype(str)
    metrics=['total_counts','k_positive','k_up','k_down','up_detection_fraction','down_detection_fraction',
             'zero_rank_baseline','up_rank_contribution','down_rank_contribution','centered_up','centered_down',
             'occupancy_component','occupancy_up','occupancy_down','positive_rank_residual','residual_bound','ucell_zero_floor_fraction',*METHODS]
    summaries=[]
    for labels,index in groups(table,columns):
        part=table.loc[index]
        summaries.append({'group':labels,'n':len(part),'minimum_group_size':minimum_group_size,
                          'rate_eligible':len(part)>=minimum_group_size,
                          'means':{k:float(part[k].mean()) for k in metrics},
                          'method_evaluable_cells':{m:int(np.isfinite(part[m]).sum()) for m in METHODS},
                          'supplied_thresholds':{m:threshold_counts(part[m],t,minimum_group_size) for m,t in thresholds.items()}})
    contrasts=[]
    if (baseline_treatment is None)!=(comparison_treatment is None):raise ValueError('Supply both treatment labels')
    if baseline_treatment is not None:
        if baseline_treatment==comparison_treatment:raise ValueError('Treatment labels must differ')
        if 'treatment' not in table:raise ValueError('Treatment comparison requires treatment metadata')
        for labels,index in groups(table,[c for c in columns if c!='treatment']):
            part=table.loc[index];a=part[part.treatment==baseline_treatment];b=part[part.treatment==comparison_treatment]
            enough=len(a)>=minimum_group_size and len(b)>=minimum_group_size
            delta={m:float(b[m].mean()-a[m].mean()) for m in metrics} if len(a) and len(b) else {}
            paired=None
            if 'pair_id' in table and len(a) and len(b):
                identifiers(a.pair_id.tolist(),'baseline pair_id');identifiers(b.pair_id.tolist(),'comparison pair_id')
                complete=set(a.pair_id)==set(b.pair_id)
                paired={'complete_matching':complete,'baseline_unpaired_n':len(set(a.pair_id)-set(b.pair_id)),
                        'comparison_unpaired_n':len(set(b.pair_id)-set(a.pair_id)),'transitions':{}}
                if complete:
                    b_aligned=b.set_index('pair_id').loc[a.pair_id]
                    for method,threshold in thresholds.items():
                        prior=a[method].to_numpy();after=b_aligned[method].to_numpy();valid=np.isfinite(prior)&np.isfinite(after)
                        n=int(valid.sum());lost=int(((prior>threshold)&(after<=threshold)&valid).sum())
                        gained=int(((prior<=threshold)&(after>threshold)&valid).sum())
                        paired['transitions'][method]={'paired_evaluable_n':n,'negative_to_positive_n':gained,'positive_to_negative_n':lost,
                            'negative_to_positive_fraction':gained/n if n>=minimum_group_size else None,
                            'positive_to_negative_fraction':lost/n if n>=minimum_group_size else None}
            contrasts.append({'stratum':labels,'baseline':baseline_treatment,'comparison':comparison_treatment,
                              'baseline_n':len(a),'comparison_n':len(b),'rate_and_flag_eligible':enough,
                              'mean_change_comparison_minus_baseline':delta,
                              'paired_observations':paired,
                              'directional_flags':{
                                  'score_rises_while_down_detection_falls':bool(delta['signed_rank']>0 and delta['down_detection_fraction']<0) if enough else None,
                                  'score_rises_while_up_detection_falls':bool(delta['signed_rank']>0 and delta['up_detection_fraction']<0) if enough else None},
                              'scope':'Descriptive arithmetic directions, not tests of mechanism or treatment efficacy.'})
    states=[]
    if (state_positive is None)!=(state_negative is None):raise ValueError('Supply both state labels')
    if state_positive is not None:
        if state_positive==state_negative:raise ValueError('State labels must differ')
        if 'state' not in table:raise ValueError('State discrimination requires state metadata')
        for labels,index in groups(table,[c for c in columns if c!='state']):
            part=table.loc[index];pos=part[part.state==state_positive];neg=part[part.state==state_negative]
            for method in METHODS:
                p=pos[method].to_numpy();n=neg[method].to_numpy();p=p[np.isfinite(p)];n=n[np.isfinite(n)]
                enough=len(p)>=minimum_group_size and len(n)>=minimum_group_size
                states.append({'stratum':labels,'method':method,'positive_state':state_positive,'negative_state':state_negative,
                               'positive_evaluable_n':len(p),'negative_evaluable_n':len(n),'rate_eligible':enough,
                               'AUROC':auc(p,n) if enough else None,
                               'positive_threshold_counts':threshold_counts(p,thresholds[method],minimum_group_size) if method in thresholds else None,
                               'negative_threshold_counts':threshold_counts(n,thresholds[method],minimum_group_size) if method in thresholds else None})
    result={'schema_version':1,'purpose':'Detection and score diagnostics for one explicitly supplied signature.',
            'scope':'No corrected score, pass/fail conclusion, conversion probability, ancestry assignment or functional claim is produced.',
            'eligibility':{'input_cells':len(cells),'positive_total_cells':int(eligible.sum()),'excluded_zero_total_cells':int((~eligible).sum()),
                           'metadata_reordered_by_cell_id':reordered,'minimum_group_size':minimum_group_size,
                           'counts_retained_below_minimum':True,'rates_and_directional_flags_suppressed_below_minimum':True},
            'signature':{'universe_n':len(universe),'up_n':len(up),'down_n':len(down),
                         'universe_genes_absent_from_input':sorted(set(universe)-set(genes)),
                         'up_genes_absent_from_input':sorted(set(up)-set(genes)),
                         'down_genes_absent_from_input':sorted(set(down)-set(genes)),
                         'missing_gene_policy':'Zero-filled inside the explicit universe; original block denominators retained.'},
            'definitions':{'score':'occupancy_component + positive_rank_residual',
                           'occupancy_component':'0.5*(detected_up/up_n - detected_down/down_n)',
                           'zero_rank_baseline':'(universe_n - detected_universe_genes)/(2*universe_n)',
                           'centered_blocks':'Signed block contributions after removing their cancelling zero-rank baselines.',
                           'ucell_form':'Fixed cap1500 and frozen universe; signed block difference clipped at0. This is not a package-equivalence claim.',
                           'unclipped_ucell_difference':'Explanatory pre-clipping difference, not a proposed replacement scorer.'},
            'groups':summaries,'treatment_contrasts':contrasts,'state_discrimination':states}
    return clean_json(result),table


def report_html(result):
    def safe(value):return html.escape(str(value),quote=True)
    def formatted(value):
        if value is None:return 'not available'
        if isinstance(value,float):return f'{value:.6g}'
        if isinstance(value,dict):return '; '.join(str(k)+'='+str(v) for k,v in value.items())
        return str(value)
    def table(rows):
        if not rows:return '<p>No comparison was requested or available.</p>'
        keys=list(rows[0]);head=''.join('<th>'+safe(k.replace('_',' '))+'</th>' for k in keys)
        body=''.join('<tr>'+''.join('<td>'+safe(formatted(row.get(k)))+'</td>' for k in keys)+'</tr>' for row in rows)
        return '<div class="scroll"><table><thead><tr>'+head+'</tr></thead><tbody>'+body+'</tbody></table></div>'
    def group_table(metrics):
        return table([{'group':r['group'],'n':r['n'],**{k:r['means'][k] for k in metrics}} for r in result['groups']])
    thresholds=[{'group':r['group'],'method':m,**v} for r in result['groups'] for m,v in r['supplied_thresholds'].items()]
    contrast_rows=[];paired_rows=[]
    for r in result['treatment_contrasts']:
        delta=r['mean_change_comparison_minus_baseline']
        contrast_rows.append({'stratum':r['stratum'],'baseline':r['baseline'],'comparison':r['comparison'],'baseline_n':r['baseline_n'],
                              'comparison_n':r['comparison_n'],**{'change_'+k:delta.get(k) for k in ['signed_rank','up_detection_fraction','down_detection_fraction','occupancy_component','positive_rank_residual']},
                              **r['directional_flags']})
        p=r.get('paired_observations')
        if p:
            for method,values in p['transitions'].items():paired_rows.append({'stratum':r['stratum'],'method':method,**values})
    state_rows=[{k:v for k,v in r.items() if not k.endswith('_threshold_counts')} for r in result['state_discrimination']]
    sections=[('Eligibility', '<pre>'+safe(json.dumps(result['eligibility'],indent=2))+'</pre>'),
              ('Signature and definitions','<pre>'+safe(json.dumps({'signature':result['signature'],'definitions':result['definitions']},indent=2))+'</pre>'),
              ('Group detection means',group_table(['total_counts','k_positive','k_up','k_down','up_detection_fraction','down_detection_fraction'])),
              ('Zero baselines and centered blocks',group_table(['zero_rank_baseline','up_rank_contribution','down_rank_contribution','centered_up','centered_down'])),
              ('Exact score decomposition',group_table(['occupancy_component','positive_rank_residual','signed_rank','residual_bound'])),
              ('Declared scoring methods',group_table([*METHODS,'ucell_zero_floor_fraction'])),('Supplied strict thresholds',table(thresholds)),
              ('Treatment or labelled comparison directions',table(contrast_rows)),('Explicitly paired threshold transitions',table(paired_rows)),('State discrimination',table(state_rows)),
              ('Provenance','<pre>'+safe(json.dumps(result.get('provenance',{}),indent=2))+'</pre>')]
    return '<!doctype html><html lang="en"><meta charset="utf-8"><title>Signature diagnostic</title><style>body{font:16px system-ui,sans-serif;max-width:1250px;margin:32px auto;padding:0 20px;color:#171717}h1,h2{line-height:1.2}.scope{padding:14px;background:#eef2f5}table{border-collapse:collapse;font-size:13px}th,td{border:1px solid #ccc;padding:7px;text-align:left;vertical-align:top}th{background:#f3f3f3}.scroll{overflow-x:auto}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:13px}</style><h1>Signature diagnostic</h1><p class="scope">'+safe(result['scope'])+'</p>'+''.join('<h2>'+safe(title)+'</h2>'+body for title,body in sections)+'</html>'


def example_inputs(directory):
    """Small entirely synthetic example: no biological data or interpretation."""
    if directory.exists():raise FileExistsError('Example input directory already exists: '+str(directory))
    directory.mkdir(parents=True)
    genes=[f'g{i}' for i in range(1603)];cells=[f'synthetic_{i:03}' for i in range(96)]
    rng=np.random.default_rng(7007);x=np.zeros((len(genes),len(cells)),dtype=np.int64)
    state=np.repeat(['state_low','state_high'],48);treatment=np.tile(np.repeat(['vehicle','exposure'],24),2)
    for j in range(len(cells)):
        x[:3,j]=rng.integers(1,7,3) if state[j]=='state_high' else rng.integers(0,2,3)
        x[3:6,j]=rng.integers(1,7,3) if state[j]=='state_low' else rng.integers(0,2,3)
        if treatment[j]=='exposure':x[4:6,j]=0
        x[10,j]=10
    paths={k:directory/name for k,name in {'matrix':'counts.npz','genes':'genes.txt','cells':'cells.txt','metadata':'metadata.tsv','signature':'signature.json'}.items()}
    sparse.save_npz(paths['matrix'],sparse.csc_matrix(x))
    paths['genes'].write_text('\n'.join(genes)+'\n',encoding='utf-8');paths['cells'].write_text('\n'.join(cells)+'\n',encoding='utf-8')
    pd.DataFrame({'cell_id':cells,'donor':'synthetic_donor','state':state,'treatment':treatment}).to_csv(paths['metadata'],sep='\t',index=False)
    paths['signature'].write_text(json.dumps({'up':genes[:3],'down':genes[3:6],'universe':genes},indent=2)+'\n',encoding='utf-8')
    return paths


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path)
    parser.add_argument('--output-prefix',default='diagnostics/signature_report',help='Relative path beneath external data root')
    parser.add_argument('--example',action='store_true',help='Generate/run tiny synthetic inputs outside the repository')
    for flag in ['matrix','genes','cells','metadata','signature','thresholds','provenance']:parser.add_argument('--'+flag,type=Path)
    parser.add_argument('--threshold',type=float,help='Supplied strict signed-rank threshold; never estimated here')
    parser.add_argument('--minimum-group-size',type=int,default=20)
    for flag in ['baseline-treatment','comparison-treatment','state-positive','state-negative']:parser.add_argument('--'+flag)
    args=parser.parse_args()
    try:
        root=resolve_data_root(args.data_root,for_write=True)
        prefix=data_path(args.output_prefix,root,for_write=True)
        outputs={kind:checked_data_destination(Path(str(prefix)+suffix),root) for kind,suffix in [('json','.json'),('html','.html'),('cells','.cells.tsv')]}
        if any(p.exists() for p in outputs.values()):raise FileExistsError('Output already exists; choose a new --output-prefix')
        if args.example:
            if any(getattr(args,k) for k in ['matrix','genes','cells','metadata','signature']):raise ValueError('--example cannot be combined with input paths')
            inputs=example_inputs(checked_data_destination(prefix.parent/(prefix.name+'_inputs'),root))
            args.baseline_treatment=args.baseline_treatment or 'vehicle';args.comparison_treatment=args.comparison_treatment or 'exposure'
            args.state_positive=args.state_positive or 'state_high';args.state_negative=args.state_negative or 'state_low'
        else:
            if not all(getattr(args,k) for k in ['matrix','genes','cells','signature']):raise ValueError('Supply --matrix, --genes, --cells and --signature, or --example')
            inputs={k:getattr(args,k).resolve() for k in ['matrix','genes','cells','signature','metadata'] if getattr(args,k)}
            for key in ['matrix','genes','cells','metadata']:
                if key in inputs:checked_data_destination(inputs[key],root)
        if inputs['matrix'].suffix!='.npz':raise ValueError('--matrix must be a scipy sparse .npz file')
        signature=json.loads(inputs['signature'].read_text(encoding='utf-8-sig'))
        thresholds=json.loads(args.thresholds.read_text(encoding='utf-8-sig')) if args.thresholds else {}
        if args.threshold is not None:
            if 'signed_rank' in thresholds:raise ValueError('Specify signed-rank threshold only once')
            thresholds['signed_rank']=args.threshold
        result,cells=diagnose(sparse.load_npz(inputs['matrix']),inputs['genes'].read_text(encoding='utf-8-sig').splitlines(),
                              inputs['cells'].read_text(encoding='utf-8-sig').splitlines(),signature,
                              pd.read_csv(inputs['metadata'],sep='\t',dtype={'cell_id':str}) if 'metadata' in inputs else None,
                              thresholds=thresholds,baseline_treatment=args.baseline_treatment,comparison_treatment=args.comparison_treatment,
                              state_positive=args.state_positive,state_negative=args.state_negative,minimum_group_size=args.minimum_group_size)
        if args.thresholds:inputs['thresholds']=args.thresholds.resolve()
        if args.provenance:inputs['provenance']=args.provenance.resolve()
        result['provenance']={'synthetic_example':args.example,
                              'input_files':{k:{'path':str(p),'sha256':sha(p)} for k,p in inputs.items()},
                              'supplied':json.loads(args.provenance.read_text(encoding='utf-8-sig')) if args.provenance else None,
                              'threshold_origin':'Supplied by caller; not fitted or validated by this diagnostic.',
                              'code_sha256':{p.name:sha(p) for p in [Path(__file__),Path(__file__).with_name('revision02_math.py'),Path(__file__).with_name('revision_score.py')]}}
        result=clean_json(result)
        prefix.parent.mkdir(parents=True,exist_ok=True)
        outputs['json'].write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        outputs['html'].write_text(report_html(result),encoding='utf-8')
        cells.to_csv(outputs['cells'],sep='\t',index=False)
        print(json.dumps({'outputs':{k:str(p) for k,p in outputs.items()},'eligibility':result['eligibility']}))
    except (ValueError,OSError,KeyError,TypeError) as error:parser.error(str(error))

if __name__=='__main__':main()
