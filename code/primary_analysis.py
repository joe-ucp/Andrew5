"""Unintegrated donor-wise diagnostic analysis under the frozen calibration gate."""
from data_paths import evidence_path, prepare_original_results
from pathlib import Path
import sys,json,itertools,hashlib
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata
from reference_analysis import CELL_CYCLE, ROOT, EV
sys.path.insert(0,str(ROOT/'scripts_primary'))
from audit_primary import load_sample

def score_sparse(x,weights):
    x=x.tocsc(); n=x.shape[0];out=np.empty(x.shape[1])
    nu=int((weights>0).sum());nd=int((weights<0).sum())
    iw=np.where(weights>0,nd,np.where(weights<0,-nu,0)).astype(np.int64)
    for j in range(x.shape[1]):
        a,b=x.indptr[j:j+2];idx=x.indices[a:b];vals=x.data[a:b]
        nz=vals>0;idx=idx[nz];vals=vals[nz];z=n-len(vals)
        half=np.rint(2*rankdata(vals,method='average')).astype(np.int64)+2*z
        numerator=np.dot(iw[idx],half-(z+1))
        out[j]=numerator/(2*n*nu*nd)
    return out

def map_universe(x,feat,universe):
    look={g:i for i,g in enumerate(universe)};rr=[];cc=[]
    for j,g in enumerate(feat.symbol):
        if g in look:rr.append(look[g]);cc.append(j)
    a=sparse.csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(universe),x.shape[0]))
    return (a@x).tocsc()

def full_gene(x,feat,g):return np.asarray(x[np.flatnonzero(feat.symbol.to_numpy()==g)].sum(axis=0)).ravel()

def simulate_scores(profiles,depths,weights,rng,n=400,beta=None,frac=0):
    """profiles universe+one outside category, each column sums to one."""
    cols=[];rows=[];vals=[]
    ai=rng.integers(0,profiles.shape[1],n)
    bi=rng.integers(0,beta.shape[1],n) if beta is not None else None
    for j,a in enumerate(ai):
        p=profiles[:,a].toarray().ravel()
        if beta is not None:p=(1-frac)*p+frac*beta[:,bi[j]].toarray().ravel()
        p=np.maximum(p,0);p/=p.sum()
        c=rng.multinomial(int(rng.choice(depths)),p)[:-1]
        ii=np.flatnonzero(c);rows.extend(ii);cols.extend([j]*len(ii));vals.extend(c[ii])
    x=sparse.csc_matrix((vals,(rows,cols)),shape=(len(weights),n))
    return score_sparse(x,weights)

def prob_profiles(x,totals):
    mat=x@sparse.diags(1/np.maximum(totals,1))
    other=np.maximum(0,1-np.asarray(mat.sum(axis=0)).ravel())
    return sparse.vstack([mat,sparse.csc_matrix(other.reshape(1,-1))],format='csc')

def group_margin(pool,base,n1,n0,rng,reps=200):
    draws=np.empty(reps)
    for j in range(reps):draws[j]=rng.choice(pool,n1).mean()-rng.choice(base,n0).mean()
    return float(np.quantile(draws,.95)),float(np.mean(draws)),float(np.std(draws,ddof=1))

def main():
    prepare_original_results()
    calibration=json.loads((EV/'calibration_status.json').read_text())
    assert calibration['technical_calibration_pass'] is False, 'This diagnostic branch requires the actual failed calibration receipt'
    cfg=json.loads((EV/'frozen_score.json').read_text());universe=cfg['universe'];weights=np.zeros(len(universe));idx={g:i for i,g in enumerate(universe)}
    weights[[idx[g] for g in cfg['up']]]=1/len(cfg['up']);weights[[idx[g] for g in cfg['down']]]=-1/len(cfg['down'])
    rows=pd.read_csv(ROOT/'primary_sample_manifest.tsv',sep='\t').to_dict('records')
    rng=np.random.default_rng(7008);cells=[];summaries=[];pools={};profiles={};all_scores={}
    for row in rows:
        x,feat,bc=load_sample(row);lib=np.asarray(x.sum(axis=0)).ravel();ng=np.asarray((x>0).sum(axis=0)).ravel();mt=np.asarray(x[feat.symbol.str.startswith('MT-').to_numpy()].sum(axis=0)).ravel()/np.maximum(lib,1)
        qc=(lib>=1000)&(ng>=400)&(mt<.2)
        markers={g:full_gene(x,feat,g) for g in ['INS','GCG','SST','PPY','GHRL','ARX','IRX2','SMOC1']}
        alpha=qc&(markers['GCG']>=3)
        beta=qc&(markers['INS']>=3)
        for g in ['INS','SST','PPY','GHRL']:alpha &=markers['GCG']>markers[g]
        for g in ['GCG','SST','PPY','GHRL']:beta &=markers['INS']>markers[g]
        chosen=alpha|beta;ix=np.flatnonzero(chosen)
        xu=map_universe(x[:,ix],feat,universe)
        scores=score_sparse(xu,weights)
        cyc=np.asarray(x[np.flatnonzero(feat.symbol.isin(CELL_CYCLE))].sum(axis=0)).ravel()/np.maximum(lib,1)
        counts=np.asarray((xu>0)[np.flatnonzero(weights)].sum(axis=0)).ravel()
        trim=lib<=np.quantile(lib[qc],.99)
        ai=alpha[ix];bi=beta[ix]
        arr=pd.DataFrame({'accession':row['accession'],'donor':row['donor'],'treatment':row['treatment_code'],'barcode':bc.barcode.iloc[ix].values,'annotation':np.where(ai,'alpha_related','beta_like'),'score':scores,'library_counts':lib[ix],'detected_genes':ng[ix],'score_genes_detected':counts,'cycle_fraction':cyc[ix],'GCG':markers['GCG'][ix],'INS':markers['INS'][ix],'SMOC1':markers['SMOC1'][ix],'strict_alpha_marker':(markers['ARX'][ix]+markers['IRX2'][ix])>0,'GCG_ge10':markers['GCG'][ix]>=10,'below_99pct_UMI':trim[ix]})
        # Fixed-depth observed-cell resampling, retaining entire library mass via outside category.
        pp=prob_profiles(xu,lib[ix]);resampled=[]
        for start in range(0,len(ix),400):
            sub=pp[:,start:start+400];subrows=[];subcols=[];subvals=[]
            for j in range(sub.shape[1]):
                p=sub[:,j].toarray().ravel();p/=p.sum();c=rng.multinomial(2000,p)[:-1];ii=np.flatnonzero(c)
                subrows.extend(ii);subcols.extend([j]*len(ii));subvals.extend(c[ii])
            xx=sparse.csc_matrix((subvals,(subrows,subcols)),shape=(len(universe),sub.shape[1]))
            resampled.extend(score_sparse(xx,weights))
        arr['score_depth2000']=resampled
        cells.append(arr)
        for label,mask in [('alpha_related',ai),('beta_like',bi)]:
            for sensitivity,extra in [('primary',np.ones(len(ix),bool)),('strict_marker',arr.strict_alpha_marker.to_numpy()),('GCG_ge10',arr.GCG_ge10.to_numpy()),('trim99',arr.below_99pct_UMI.to_numpy())]:
                if label=='beta_like' and sensitivity!='primary':continue
                use=mask&extra
                if use.sum():summaries.append({'accession':row['accession'],'donor':row['donor'],'treatment':row['treatment_code'],'annotation':label,'sensitivity':sensitivity,'n':int(use.sum()),'mean_score':float(scores[use].mean()),'median_score':float(np.median(scores[use])),'mean_depth2000':float(arr.score_depth2000.to_numpy()[use].mean()),'median_library_counts':float(np.median(lib[ix][use])),'mean_detected_score_genes':float(counts[use].mean())})
        key=(row['donor'],row['treatment_code']);all_scores[key]=arr
        if row['treatment_code']=='D':
            profiles[row['donor']]={'alpha':pp[:,ai],'beta':pp[:,bi],'depths':lib[ix][ai],'cycle':cyc[ix][ai],'scores':scores[ai]}
        print(row['accession'],row['donor'],row['treatment_code'],'alpha',int(ai.sum()),'beta',int(bi.sum()),flush=True)
        del x,xu,pp
    c=pd.concat(cells,ignore_index=True);c.to_csv(evidence_path('primary_cell_scores.tsv.gz'),sep='\t',index=False)
    s=pd.DataFrame(summaries);s.to_csv(EV/'sample_score_summary.tsv',sep='\t',index=False)
    result=[]
    for donor in sorted(s.donor.unique()):
        for label in ['alpha_related','beta_like']:
            for sensitivity in ['primary','strict_marker','GCG_ge10','trim99']:
                ss=s[(s.donor==donor)&(s.annotation==label)&(s.sensitivity==sensitivity)]
                if len(ss)==0:continue
                ctl=ss[ss.treatment=='D'].iloc[0]
                for _,tr in ss[ss.treatment!='D'].iterrows():
                    result.append({'donor':donor,'contrast':tr.treatment+'-D','annotation':label,'sensitivity':sensitivity,'n_control':int(ctl.n),'n_treated':int(tr.n),'mean_control':ctl.mean_score,'mean_treated':tr.mean_score,'delta':tr.mean_score-ctl.mean_score,'delta_depth2000':tr.mean_depth2000-ctl.mean_depth2000,'interpretation':'diagnostic_only_calibration_failed'})
    rr=pd.DataFrame(result)
    nullrows=[]
    for di,(donor,p) in enumerate(profiles.items()):
        local=np.random.default_rng(7100+di)
        base=simulate_scores(p['alpha'],p['depths'],weights,local)
        scenarios={}
        for frac in [0,.01,.05,.10,.20,.50]:scenarios['beta_RNA_'+str(frac)]=simulate_scores(p['alpha'],p['depths'],weights,local,beta=p['beta'],frac=frac)
        for frac in [.01,.05,.10]:
            mix=base.copy();hit=local.choice(len(mix),int(round(frac*len(mix))),replace=False);mix[hit]=local.choice(scenarios['beta_RNA_0.5'],len(hit));scenarios['doublet_fraction_'+str(frac)]=mix
        cycle=p['cycle'];q=np.quantile(cycle,[.25,.75]);low=cycle<=q[0];high=(cycle>=q[1])&(cycle>0)&~low
        cycle_valid=(cycle>0).sum()>=20 and low.sum()>=20 and high.sum()>=20
        if cycle_valid:
            lowpool=simulate_scores(p['alpha'][:,low],p['depths'],weights,local);highpool=simulate_scores(p['alpha'][:,high],p['depths'],weights,local)
            for frac in [0,.25,.5,.75,1]:scenarios['cycle_high_weight_'+str(frac)]=np.r_[local.choice(lowpool,int(400*(1-frac))),local.choice(highpool,int(400*frac))]
        for contrast in ['H_GLP-D','H-D','H_LY-D']:
            r=rr[(rr.donor==donor)&(rr.contrast==contrast)&(rr.annotation=='alpha_related')&(rr.sensitivity=='primary')].iloc[0]
            for model,pool in scenarios.items():
                eps,mean,sd=group_margin(pool,base,int(r.n_treated),int(r.n_control),local)
                nullrows.append({'donor':donor,'contrast':contrast,'model':model,'technical_margin_q95':eps,'simulated_mean_delta':mean,'simulated_SD_group_delta':sd,'observed_delta':r.delta,'exceeds_margin':bool(r.delta>eps),'actual_n_control':r.n_control,'actual_n_treated':r.n_treated,'synthetic_pool_n':400,'MC_group_replicates':200,'cycle_calibrated':bool(cycle_valid),'assumption_status':'declared_sensitivity_not_empirically_calibrated_contamination'})
        print(donor,'nulls complete','cycling',cycle_valid,flush=True)
    nr=pd.DataFrame(nullrows);nr.to_csv(EV/'countermodel_results.tsv',sep='\t',index=False)
    rr.to_csv(EV/'primary_results.tsv',sep='\t',index=False)
    ds=rr[(rr.contrast=='H_GLP-D')&(rr.annotation=='alpha_related')&(rr.sensitivity=='primary')].sort_values('donor').delta.to_numpy()
    signs=np.array(list(itertools.product([-1,1],repeat=4)));vals=(signs*ds).mean(axis=1);observed=ds.mean();tail=float(np.mean(vals>=observed-1e-14))
    pd.DataFrame({'signs':[' '.join(map(str,x)) for x in signs],'mean_signed_delta':vals,'at_least_observed':vals>=observed-1e-14}).to_csv(EV/'sign_flip_reference.tsv',sep='\t',index=False)
    (EV/'primary_status.json').write_text(json.dumps({'delta_by_donor':dict(zip(sorted(s.donor.unique()),ds.tolist())),'mean_delta':float(observed),'descriptive_one_sided_sign_flip_tail':tail,'randomization_test':False,'technical_calibration_pass':calibration['technical_calibration_pass'],'biological_specificity':'UNRESOLVED','joint_prediction':'UNRESOLVED','primary_interpretation':'prespecified diagnostic after failed external calibration'},indent=2))
    print((EV/'primary_status.json').read_text())

if __name__ == '__main__':
    import argparse
    from pathlib import Path
    from public_runtime import execute_stage
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    execute_stage('primary', args.data_root, args.output_dir)
