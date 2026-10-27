"""Pure gene summaries and vehicle-matched programs; no file access or fitting."""
import math
import numpy as np
from scipy import sparse
from scipy.stats import rankdata


def summarize_genes(counts, totals):
    """Frozen-universe genes by selected cells, with original whole-library totals.

    Positive-expression summaries condition on detection separately in each arm.
    centered_rank_sum is the exact integer sum of 2*r_positive-K-1 for each gene.
    """
    x=sparse.csc_matrix(counts,copy=True);x.sum_duplicates();x.eliminate_zeros();x.sort_indices()
    totals=np.asarray(totals)
    if x.shape[1]==0 or totals.shape!=(x.shape[1],):raise ValueError('Nonempty cells and aligned totals required')
    if not np.isfinite(x.data).all() or (x.data<0).any() or (x.data!=np.floor(x.data)).any():raise ValueError('Nonnegative integer counts required')
    if not np.isfinite(totals).all() or (totals<=0).any() or (totals!=np.floor(totals)).any():raise ValueError('Positive whole-library integer totals required')
    if (np.asarray(x.sum(axis=0)).ravel()>totals).any():raise ValueError('Universe counts exceed whole-library totals')
    g,n=x.shape;detected=np.diff(x.tocsr().indptr).astype(np.int64)
    raw_sum=np.asarray(x.sum(axis=1)).ravel().astype(np.int64)
    norm=x @ sparse.diags(10000/totals)
    norm_sum=np.asarray(norm.sum(axis=1)).ravel()
    norm.data=np.log1p(norm.data)
    log_sum=np.asarray(norm.sum(axis=1)).ravel()
    centered=np.zeros(g,dtype=np.int64)
    for j in range(n):
        a,b=x.indptr[j:j+2];values=x.data[a:b];indices=x.indices[a:b]
        twice=np.rint(2*rankdata(values,method='average')).astype(np.int64)
        centered[indices]+=twice-len(values)-1
    def conditional(values):
        return np.divide(values,detected,out=np.full(g,np.nan),where=detected>0)
    return {'n':np.int64(n),'G':np.int64(g),'detected_n':detected,'detection_fraction':detected/n,
            'mean_CP10K':norm_sum/n,'positive_mean_raw':conditional(raw_sum),
            'positive_mean_CP10K':conditional(norm_sum),'positive_mean_log1p_CP10K':conditional(log_sum),
            'centered_rank_sum':centered}


def program_numerators(summary,programs):
    """50-up then 50-down indices; mean score numerator / (100*G*n)."""
    programs=np.asarray(programs,dtype=np.int64)
    if programs.ndim!=2 or programs.shape[1]!=100:raise ValueError('Exactly 50 up and 50 down genes required')
    if (programs<0).any() or (programs>=int(summary['G'])).any():raise ValueError('Program index outside frozen universe')
    if any(len(np.unique(row))!=100 for row in programs):raise ValueError('Each program must have 100 disjoint genes')
    def contrast(values):return values[programs[:,:50]].sum(axis=1)-values[programs[:,50:]].sum(axis=1)
    detection=contrast(summary['detected_n']);residual=contrast(summary['centered_rank_sum'])
    return {'occupancy':detection,'residual':residual,'score':int(summary['G'])*detection+residual}


def donor_contrasts(control,treated,programs):
    """Integer numerators preserve all ties under common per-donor denominators."""
    if control['G']!=treated['G']:raise ValueError('Universe mismatch')
    a=program_numerators(control,programs);b=program_numerators(treated,programs)
    n0,n1,g=int(control['n']),int(treated['n']),int(control['G'])
    return {key:{'numerators':b[key]*n0-a[key]*n1,
                 'denominator':100*n0*n1*(1 if key=='occupancy' else g)} for key in a}


def equal_donor_mean(contrasts,component):
    denominators=[int(c[component]['denominator']) for c in contrasts]
    common=math.lcm(*denominators)
    values=[c[component]['numerators'] for c in contrasts]
    numerators=np.array([sum(int(values[d][j])*(common//denominators[d]) for d in range(len(values)))
                         for j in range(len(values[0]))],dtype=object)
    return {'numerators':numerators,'denominator':common*len(contrasts)}


def match_candidates(genes,original,p,L):
    """Fixed baseline calipers, distance and support; treatment information absent."""
    genes=np.asarray(genes);original=np.asarray(original,dtype=np.int64)
    p=np.asarray(p,float);L=np.asarray(L,float)
    if p.shape!=L.shape or p.shape!=(4,len(genes)) or len(original)!=100:raise ValueError('Four donors and 100 original slots required')
    if not np.isfinite(p).all() or not np.isfinite(L).all():raise ValueError('Finite matching covariates required')
    allowed=np.ones(len(genes),bool);allowed[original]=False
    pbar=p.mean(axis=0);Lbar=L.mean(axis=0);candidates=[];support=[]
    for slot,gene in enumerate(original):
        eligible=np.flatnonzero(allowed&(np.abs(pbar-pbar[gene])<=.025)&(np.abs(Lbar-Lbar[gene])<=.5))
        distance=(((p[:,eligible]-p[:,gene,None])/.05)**2+(L[:,eligible]-L[:,gene,None])**2).sum(axis=0)
        ordering=np.lexsort((genes[eligible],distance));chosen=eligible[ordering[:200]]
        candidates.append(chosen)
        support.append({'slot':slot,'gene':str(genes[gene]),'block':'up' if slot<50 else 'down',
                        'eligible_candidates':len(eligible),'retained_candidates':len(chosen),
                        'support_pass':len(chosen)>=10})
    return candidates,support


def draw_matched_programs(original,candidates,p,L,seed=8202,target=500,max_attempts=50000):
    """Baseline-only rejection sampling, with explicit failure and no relaxation."""
    original=np.asarray(original,dtype=np.int64);p=np.asarray(p);L=np.asarray(L)
    if len(original)!=100 or len(candidates)!=100:raise ValueError('Exactly 100 labeled slots required')
    if any(len(c)<10 for c in candidates):
        return np.empty((0,100),dtype=np.int64),[],{'status':'FAILED_CANDIDATE_SUPPORT','attempts':0,'accepted':0,'assignment_rejections':0,'quality_rejections':0}
    rng=np.random.default_rng(seed);accepted=[];quality=[];assignment_reject=quality_reject=0
    baseline_p=np.stack([p[:,original[:50]].mean(axis=1),p[:,original[50:]].mean(axis=1)],axis=1)
    baseline_L=np.stack([L[:,original[:50]].mean(axis=1),L[:,original[50:]].mean(axis=1)],axis=1)
    for attempt in range(1,max_attempts+1):
        selected=np.full(100,-1,dtype=np.int64);used=np.zeros(p.shape[1],bool)
        for slot in rng.permutation(100):
            candidate=np.asarray(candidates[slot]);available=candidate[~used[candidate]]
            if not len(available):break
            pick=int(rng.choice(available));selected[slot]=pick;used[pick]=True
        if (selected<0).any():assignment_reject+=1;continue
        candidate_p=np.stack([p[:,selected[:50]].mean(axis=1),p[:,selected[50:]].mean(axis=1)],axis=1)
        candidate_L=np.stack([L[:,selected[:50]].mean(axis=1),L[:,selected[50:]].mean(axis=1)],axis=1)
        delta_p=candidate_p-baseline_p;delta_L=candidate_L-baseline_L
        if (np.abs(delta_p)>.025).any() or (np.abs(delta_L)>.25).any():quality_reject+=1;continue
        accepted.append(selected)
        for donor in range(4):
            for block in range(2):
                quality.append({'program':len(accepted),'donor':'D'+str(donor+1),'block':['up','down'][block],
                                'original_detection':baseline_p[donor,block],'matched_detection':candidate_p[donor,block],
                                'delta_detection':delta_p[donor,block],'original_log_abundance':baseline_L[donor,block],
                                'matched_log_abundance':candidate_L[donor,block],'delta_log_abundance':delta_L[donor,block]})
        if len(accepted)==target:break
    programs=np.stack(accepted) if accepted else np.empty((0,100),dtype=np.int64)
    return programs,quality,{'status':'COMPLETE' if len(accepted)==target else 'FAILED_TARGET_NOT_REACHED',
                            'attempts':attempt,'accepted':len(accepted),'assignment_rejections':assignment_reject,
                            'quality_rejections':quality_reject,'seed':seed,'target':target,'max_attempts':max_attempts,
                            'unique_ordered_programs':len({tuple(x) for x in programs})}
