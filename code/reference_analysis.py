"""Frozen external program and explicitly scoped preparation-level calibration."""
from data_paths import evidence_path, prepare_original_results
from pathlib import Path
import gzip,tarfile,io,json,hashlib,datetime,os
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from data_paths import data_path,configure_data_root_cli
ROOT=Path(__file__).resolve().parents[1]
EV=ROOT/'evidence'
from public_runtime import frozen_score, scientific_file
SEED=7007
CELL_CYCLE='MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA CENPU PIMREG JPT1 CCNA2 CCNB1'.split()
EXCLUDE=set(CELL_CYCLE+['INS','GCG','PDX1','MAFA','IAPP','SST','PPY','GHRL','GFP','EGFP','EYFP','eGFP','ARX','IRX2','MAFB','TTR'])
PAIRS=[('IIDP1697','6_S6','12_S12'),('IIDP1738','7_S7','13_S13'),('IIDP1770','9_S9','15_S15'),('IIDP1776','10_S10','16_S16')]

def references():
    bulk=pd.read_csv(data_path('data_reference/GSE117454_RNA_seq_matrix.txt.gz',must_exist=True),sep='\t',index_col=0).groupby(level=0,sort=True).sum()
    sc={}
    with tarfile.open(data_path('data_reference/GSE123844_RAW.tar',must_exist=True)) as t:
        for m in t:
            key={'GSM3513969':'alpha','GSM3513970':'beta','GSM3513971':'PM','GSM3513972':'graft'}[m.name[:10]]
            sc[key]=pd.read_csv(io.BytesIO(gzip.decompress(t.extractfile(m).read())),sep='\t',index_col=0).groupby(level=0,sort=True).sum()
    return bulk,sc

def scores_dense(x,up,down):
    # x: genes by cells, universe in frozen order; average ranks for tied zero counts.
    out=[]
    for start in range(0,x.shape[1],128):
        h=np.rint(2*rankdata(x[:,start:start+128],axis=0,method='average')).astype(np.int64)
        numerator=h[up].sum(axis=0)*len(down)-h[down].sum(axis=0)*len(up)
        out.extend((numerator/(2*x.shape[0]*len(up)*len(down))).tolist())
    return np.asarray(out)

def auc(pos,neg):
    both=np.r_[pos,neg];r=rankdata(both)
    return float((r[:len(pos)].sum()-len(pos)*(len(pos)+1)/2)/(len(pos)*len(neg)))

def bootstrap_ci(pos,neg,rng,n=500):
    vals=[auc(rng.choice(pos,len(pos)),rng.choice(neg,len(neg))) for _ in range(n)]
    return np.quantile(vals,[.025,.975]).tolist()

def main():
    prepare_original_results()
    bulk,sc=references()
    EV.mkdir(exist_ok=True)
    sm=pd.read_csv(ROOT/'primary_sample_manifest.tsv',sep='\t').iloc[0]
    pf=pd.read_csv(data_path(sm.features_path,must_exist=True),sep='\t',header=None)[1].astype(str)
    cols=[c for _,a,b in PAIRS for c in (a,b)]
    cpm=bulk[cols].div(bulk[cols].sum(),axis=1)*1e6
    universe=sorted(g for g in bulk.index.intersection(sc['alpha'].index).intersection(pf) if g not in EXCLUDE and not g.startswith(('MT-','RPL','RPS')) and not any(v in g.upper() for v in ['EGFP','EYFP','VECTOR']) and (cpm.loc[g]>=1).sum()>=4)
    log=np.log2(cpm.loc[universe]+1)
    diff=np.stack([(log[b]-log[a]).values for _,a,b in PAIRS],axis=1)
    med=np.median(diff,axis=1)
    up_idx=sorted(np.flatnonzero((diff>0).all(axis=1)),key=lambda i:(-med[i],universe[i]))[:50]
    down_idx=sorted(np.flatnonzero((diff<0).all(axis=1)),key=lambda i:(med[i],universe[i]))[:50]
    config={'seed':SEED,'universe':universe,'up':[universe[i] for i in up_idx],'down':[universe[i] for i in down_idx],'pairs':PAIRS,'exclusions':sorted(EXCLUDE),'score':'mean_fractional_rank_up-minus-mean_fractional_rank_down','plan_sha256':frozen_score()['plan_sha256']}
    # Preserve historical score identity; program construction precedes calibration in this execution.
    cpath=EV/'frozen_score.json'; payload=json.dumps(config,sort_keys=True,indent=2)
    retained=scientific_file('evidence/frozen_score.json')
    if retained.read_text()!=payload or (cpath.exists() and cpath.read_text()!=payload):
        raise RuntimeError('Frozen score mismatch; explicit amendment required')
    # Preserve the registered bytes across platforms; serialization is not a refit.
    if not cpath.exists():cpath.write_bytes(retained.read_bytes())
    table=pd.DataFrame({'gene':universe,'median_paired_log2_CPM_change':med,'selected':['up' if i in up_idx else 'down' if i in down_idx else '' for i in range(len(universe))]})
    for i,(d,a,b) in enumerate(PAIRS):table[d+'_log2_CPM_change']=diff[:,i]
    table.to_csv(evidence_path('reference_feature_selection.tsv'),sep='\t',index=False)
    rng=np.random.default_rng(SEED);cell_rows=[];native={}
    for key,d in sc.items():
        arr=d.reindex(universe,fill_value=0).to_numpy()
        v=scores_dense(arr,up_idx,down_idx);native[key]=v
        for cell,score,depth in zip(d.columns,v,d.sum()):cell_rows.append({'preparation':key,'cell':cell,'depth':int(depth),'score':score})
    pd.DataFrame(cell_rows).to_csv(evidence_path('reference_cell_scores.tsv'),sep='\t',index=False)
    summary=[{'comparison':'native_PM_vs_alpha','depth':'native','auc':auc(native['PM'],native['alpha']),'positive_n':len(native['PM']),'negative_n':len(native['alpha'])}]
    sim_rows=[];gate_rows=[]
    # Full count vectors preserve hormones and total counts in multinomial simulation;
    # excluded genes consume depth but do not contribute to the score.
    genes=sc['alpha'].index; uix=genes.get_indexer(universe)
    probs={k:v.reindex(genes,fill_value=0).to_numpy(dtype=float) for k,v in sc.items()}
    for k,v in probs.items():probs[k]=v/np.maximum(v.sum(axis=0),1)
    for depth in [2000,5000,10000]:
        sims={}
        for key in ['alpha','PM']:
            picks=rng.integers(0,probs[key].shape[1],400)
            x=np.stack([rng.multinomial(depth,probs[key][:,i]) for i in picks],axis=1)
            sims[key]=scores_dense(x[uix],up_idx,down_idx)
        for frac in [0,.01,.05,.1,.2,.5]:
            ai=rng.integers(0,probs['alpha'].shape[1],400)
            bi=rng.integers(0,probs['beta'].shape[1],400)
            x=np.stack([rng.multinomial(depth,(1-frac)*probs['alpha'][:,a]+frac*probs['beta'][:,b]) for a,b in zip(ai,bi)],axis=1)
            sims[f'mixture_{frac:g}']=scores_dense(x[uix],up_idx,down_idx)
        null=np.concatenate([v for k,v in sims.items() if k!='PM'])
        threshold=float(np.quantile(null,.95));sensitivity=float(np.mean(sims['PM']>threshold))
        aucs=[];fprs=[]
        for key,v in sims.items():
            for j,s in enumerate(v):sim_rows.append({'depth':depth,'model':key,'draw':j,'score':s})
            if key!='PM':
                a=auc(sims['PM'],v);aucs.append(a);fpr=float(np.mean(v>threshold));fprs.append(fpr)
                summary.append({'comparison':'PM_vs_'+key,'depth':depth,'auc':a,'positive_n':400,'negative_n':400,'null_false_positive_rate':fpr,'threshold':threshold,'positive_sensitivity':sensitivity})
        passed=min(aucs)>=.8 and sensitivity>=.8
        gate_rows.append({'depth':depth,'minimum_auc':min(aucs),'threshold':threshold,'positive_sensitivity':sensitivity,'maximum_individual_null_FPR':max(fprs),'pass':passed})
    pd.DataFrame(summary).to_csv(EV/'calibration_results.tsv',sep='\t',index=False)
    pd.DataFrame(sim_rows).to_csv(evidence_path('calibration_simulations.tsv.gz'),sep='\t',index=False)
    pd.DataFrame(gate_rows).to_csv(EV/'calibration_gate.tsv',sep='\t',index=False)
    status={'technical_calibration_pass':bool(all(r['pass'] for r in gate_rows)),'biological_specificity':'UNRESOLVED: reference does not label adverse beta-origin lineage','calibration_unit':'preparation, not donor-independent cell replication','universe_n':len(universe),'up_n':len(up_idx),'down_n':len(down_idx),'native_auc':auc(native['PM'],native['alpha']),'gate':gate_rows}
    (EV/'calibration_status.json').write_text(json.dumps(status,indent=2));print(json.dumps(status,indent=2))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    from public_runtime import execute_stage
    execute_stage('reference', args.data_root, args.output_dir)
