"""Exploratory arithmetic influence of retained centered gene contributions.

Subtracting a contribution does not remove a matrix row, rerank a cell, refit a
signature, or predict an intervention. Universe and 50/50 denominators are fixed.
The chosen allocation is occupancy plus centered positive-rank remainder; raw
baseline-containing individual ranks would not define the same subtraction.
"""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import math
import numpy as np
import pandas as pd
from data_paths import resolve_data_root
from public_runtime import scientific_file, parameters_sha256

ROOT=Path(__file__).resolve().parents[1]
EV=ROOT/'evidence/revision_03'
ATOL=2e-15
NEAR_ZERO=1e-12
OUTPUTS=['influence_single.tsv','influence_topk.tsv','influence_selected.tsv',
         'influence_summary.tsv','influence_status.json']


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sign(value):return 'positive' if value>0 else ('negative' if value<0 else 'zero')


def calculate_influence(frame,native=None,max_k=5):
    """Pure calculation; exact parsed-float ordering, alphabetical gene ties."""
    if not isinstance(max_k,int) or isinstance(max_k,bool) or max_k<1:raise ValueError('Positive integer max_k required')
    required=['donor','gene','block','delta_signed_occupancy','delta_signed_residual',
              'delta_signed_zero_baseline_rank']
    if not set(required)<=set(frame):raise ValueError('Missing retained gene contribution columns')
    data=frame[required].copy()
    if data.empty or data[['donor','gene','block']].isna().any().any():raise ValueError('Nonempty labeled contributions required')
    if data.duplicated(['donor','gene']).any():raise ValueError('Duplicate donor/gene contribution')
    if not data.block.isin(['up','down']).all():raise ValueError('Unknown score block')
    numeric=data[required[3:]].to_numpy(dtype=float)
    if not np.isfinite(numeric).all():raise ValueError('Finite retained contributions required')
    data['contribution']=[math.fsum([o,r]) for o,r in zip(data.delta_signed_occupancy,data.delta_signed_residual)]
    errors=np.abs(data.contribution-data.delta_signed_zero_baseline_rank)
    if (errors>ATOL).any():raise AssertionError('Canonical occupancy plus remainder differs from retained centered term')
    single=[];topk=[];selected=[];summary=[]
    for donor,part in data.groupby('donor',sort=True):
        part=part.sort_values('gene',kind='stable')
        total=math.fsum(part.contribution);occupancy=math.fsum(part.delta_signed_occupancy)
        remainder=math.fsum(part.delta_signed_residual)
        native_delta=float(native[donor]) if native is not None else total
        if abs(total-native_delta)>ATOL:raise AssertionError('Contribution sum differs from retained original donor delta')
        donor_single=[]
        for row in part.itertuples(index=False):
            residual=math.fsum([total,-row.contribution])
            after_occupancy=math.fsum([occupancy,-row.delta_signed_occupancy])
            after_remainder=math.fsum([remainder,-row.delta_signed_residual])
            if abs(residual-math.fsum([after_occupancy,after_remainder]))>ATOL:raise AssertionError('Single subtraction component identity failed')
            record={'donor':donor,'gene':row.gene,'block':row.block,'native_contrast':native_delta,
                    'component_sum':total,'occupancy_contribution':row.delta_signed_occupancy,
                    'rank_remainder_contribution':row.delta_signed_residual,'total_contribution':row.contribution,
                    'positive_term':row.contribution>0,'remaining_occupancy':after_occupancy,
                    'remaining_rank_remainder':after_remainder,'single_residual':residual,
                    'remaining_sign':sign(residual),'near_zero_1e_12':abs(residual)<NEAR_ZERO,
                    'strict_sign_reversal':native_delta*residual<0}
            donor_single.append(record);single.append(record)
        ranked=part[part.contribution>0].sort_values(['contribution','gene'],ascending=[False,True],kind='stable')
        positive_n=len(ranked);prefix=[];donor_top=[]
        for rank,row in enumerate(ranked.head(max_k).itertuples(index=False),1):
            prefix.append(row)
            removed=math.fsum(r.contribution for r in prefix)
            removed_o=math.fsum(r.delta_signed_occupancy for r in prefix)
            removed_r=math.fsum(r.delta_signed_residual for r in prefix)
            remaining=math.fsum([total,-removed])
            after_o=math.fsum([occupancy,-removed_o]);after_r=math.fsum([remainder,-removed_r])
            if abs(remaining-math.fsum([after_o,after_r]))>ATOL:raise AssertionError('Cumulative subtraction component identity failed')
            selected.append({'donor':donor,'rank':rank,'gene':row.gene,'block':row.block,
                             'occupancy_contribution':row.delta_signed_occupancy,
                             'rank_remainder_contribution':row.delta_signed_residual,
                             'total_contribution':row.contribution,'available_positive_terms':positive_n})
            record={'donor':donor,'k':rank,'status':'evaluated','available_positive_terms':positive_n,
                    'selected_genes':';'.join(r.gene for r in prefix),'selected_blocks':';'.join(r.block for r in prefix),
                    'native_contrast':native_delta,'component_sum':total,
                    'removed_occupancy':removed_o,'removed_rank_remainder':removed_r,'removed_total':removed,
                    'remaining_occupancy':after_o,'remaining_rank_remainder':after_r,
                    'topk_residual':remaining,'remaining_sign':sign(remaining),
                    'near_zero_1e_12':abs(remaining)<NEAR_ZERO,'strict_sign_reversal':native_delta*remaining<0}
            donor_top.append(record);topk.append(record)
        for requested_k in range(min(positive_n,max_k)+1,max_k+1):
            record={'donor':donor,'k':requested_k,'status':'insufficient_positive_terms',
                    'available_positive_terms':positive_n,
                    'selected_genes':';'.join(r.gene for r in prefix),'selected_blocks':';'.join(r.block for r in prefix),
                    'native_contrast':native_delta,'component_sum':total,
                    'removed_occupancy':np.nan,'removed_rank_remainder':np.nan,'removed_total':np.nan,
                    'remaining_occupancy':np.nan,'remaining_rank_remainder':np.nan,
                    'topk_residual':np.nan,'remaining_sign':'unavailable',
                    'near_zero_1e_12':None,'strict_sign_reversal':None}
            donor_top.append(record);topk.append(record)
        smallest=min(donor_single,key=lambda r:(r['single_residual'],r['gene']))
        five=next((r for r in donor_top if r['k']==5 and r['status']=='evaluated'),None)
        summary.append({'donor':donor,'native_contrast':native_delta,'component_sum':total,
                        'native_sum_absolute_error':abs(native_delta-total),'original_occupancy':occupancy,
                        'original_rank_remainder':remainder,'gene_terms':len(part),'available_positive_terms':positive_n,
                        'min_single_residual':smallest['single_residual'],'min_single_removed_gene':smallest['gene'],
                        'min_single_removed_block':smallest['block'],'min_single_sign':smallest['remaining_sign'],
                        'single_strict_sign_reversals':sum(r['strict_sign_reversal'] for r in donor_single),
                        'single_zero_residuals':sum(r['remaining_sign']=='zero' for r in donor_single),
                        'single_near_zero_residuals':sum(r['near_zero_1e_12'] for r in donor_single),
                        'top5_available':five is not None,'top5_genes':five['selected_genes'] if five else '',
                        'top5_blocks':five['selected_blocks'] if five else '',
                        'top5_removed_total':five['removed_total'] if five else np.nan,
                        'top5_residual':five['topk_residual'] if five else np.nan,
                        'top5_sign':five['remaining_sign'] if five else 'unavailable',
                        'top5_near_zero_1e_12':five['near_zero_1e_12'] if five else None,
                        'top5_remaining_occupancy':five['remaining_occupancy'] if five else np.nan,
                        'top5_remaining_rank_remainder':five['remaining_rank_remainder'] if five else np.nan})
    return pd.DataFrame(single),pd.DataFrame(topk),pd.DataFrame(selected),pd.DataFrame(summary)


def calculate_from_retained_inputs():
    """Return deterministic tables after validating membership and old totals.

    Registered scientific input hashes are verified without private execution records.
    Uses the active generated workspace when configured, otherwise retained tables.
    """
    genes=pd.read_csv(scientific_file('evidence/revision_02/gene_donor_contrasts.tsv'),sep='\t',float_precision='round_trip')
    conf=json.loads((scientific_file('evidence/frozen_score.json')).read_text())
    if len(conf['universe'])!=12051 or len(conf['up'])!=50 or len(conf['down'])!=50:raise ValueError('Original universe and50/50 denominators required')
    if sorted(genes.donor.unique())!=['D1','D2','D3','D4'] or len(genes)!=400:raise ValueError('All100 genes for each of four donors required')
    for donor,part in genes.groupby('donor'):
        for block in ['up','down']:
            if set(part.loc[part.block==block,'gene'])!=set(conf[block]) or int((part.block==block).sum())!=50:
                raise ValueError('Original gene membership changed: '+donor+' '+block)
    old=pd.read_csv(scientific_file('evidence/revision_01/observed_component_contrasts.tsv'),sep='\t',float_precision='round_trip')
    old=old[old.population=='alpha_broad'].pivot(index='donor',columns='component',values='delta_mean')
    baseline=pd.read_csv(scientific_file('evidence/revision_02/zero_baseline_contrasts.tsv'),sep='\t',float_precision='round_trip').pivot(index='donor',columns='component',values='delta_mean')
    for donor,part in genes.groupby('donor'):
        for component,column in [('occupancy_component','delta_signed_occupancy'),('positive_rank_residual','delta_signed_residual')]:
            if abs(math.fsum(part[column])-old.loc[donor,component])>ATOL:raise AssertionError('Original component total does not reconcile')
        if abs(old.loc[donor,'score']-baseline.loc[donor,'score'])>ATOL:raise AssertionError('Retained original score checks disagree')
        if abs(math.fsum(math.fsum([o,r]) for o,r in zip(part.delta_signed_occupancy,part.delta_signed_residual))-baseline.loc[donor,'score'])>ATOL:
            raise AssertionError('Canonical total differs from retained Revision02 score')
    return calculate_influence(genes,old.score.to_dict())


def main(output_dir):
    """Regenerate compact influence tables externally; never modify published results."""
    destination=resolve_data_root(output_dir,for_write=True)
    for name in OUTPUTS:
        if (destination/name).exists():raise FileExistsError('Output already exists: '+str(destination/name))
    tables=calculate_from_retained_inputs()
    destination.mkdir(parents=True,exist_ok=True)
    for frame,name in zip(tables,OUTPUTS[:4]):frame.to_csv(destination/name,sep='\t',index=False)
    inputs=['evidence/frozen_score.json','evidence/revision_02/gene_donor_contrasts.tsv',
            'evidence/revision_01/observed_component_contrasts.tsv','evidence/revision_02/zero_baseline_contrasts.tsv']
    status={'status':'COMPLETE_EXPLORATORY_ARITHMETIC_INFLUENCE',
            'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'parameters_sha256':parameters_sha256(),
            'input_sha256':{name:sha(scientific_file(name)) for name in inputs},
            'code_sha256':sha(Path(__file__)),
            'rows':dict(zip(['single','topk','selected','summary'],map(len,tables))),
            'maximum_native_sum_error':float(tables[3].native_sum_absolute_error.max()),
            'reconciliation_atol':ATOL,'reconciliation_rtol':0,
            'contribution':'delta_signed_occupancy + delta_signed_residual; zero-baseline-centered allocation',
            'remaining_delta':'math.fsum(all100 canonical terms) minus selected canonical contributions; no rerank or renormalization',
            'rank_tie_rule':'Exact parsed-float equality, then alphabetical symbol',
            'near_zero_flag':'abs(remaining) < 1e-12; strict signs are reported separately',
            'outcome_informed':True,
            'scope':'Descriptive outcome-selected arithmetic influence; not gene intervention, refitted signature, causal necessity or statistical significance.',
            'output_sha256':{name:sha(destination/name) for name in OUTPUTS[:4]}}
    (destination/'influence_status.json').write_text(json.dumps(status,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return status


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True,help='Absolute external output folder; existing outputs are never replaced')
    args=parser.parse_args()
    print(json.dumps(main(args.output_dir),indent=2))
