"""Deterministic scientific figures from full-precision retained results."""
from data_paths import evidence_path
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reference_analysis import ROOT,EV
OUT=ROOT/'manuscript/figures'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
C={'alpha':'#49667c','PM':'#b55a24','beta':'#668a50','graft':'#9a8c9a'}

def save(fig,name):
    fig.savefig(OUT/(name+'.png'),bbox_inches='tight',metadata={'Software':'STRIKE007 deterministic plotting'})
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight',metadata={'Creator':'UCP Technology LLC','CreationDate':None,'ModDate':None})
    plt.close(fig)

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    r=pd.read_csv(evidence_path('reference_cell_scores.tsv'),sep='\t');cp=pd.read_csv(EV/'calibration_operating_points.tsv',sep='\t');cr=pd.read_csv(EV/'calibration_results.tsv',sep='\t')
    fig=plt.figure(figsize=(7.3,6.4));grid=fig.add_gridspec(2,2);ax=[fig.add_subplot(grid[0,:]),fig.add_subplot(grid[1,0]),fig.add_subplot(grid[1,1])]
    for k in ['alpha','PM','beta']:
        x=np.sort(r[r.preparation==k].score);ax[0].step(x,np.arange(1,len(x)+1)/len(x),where='post',label=f'{k} (n={len(x)})',color=C[k])
    ax[0].set(xlabel='Frozen rank score',ylabel='Cumulative fraction of captured cells',title='A  Native preparations');ax[0].legend(frameon=False,fontsize=8)
    for name,lab,style in [('alpha_only','Alpha control','-'),('low_mixtures','Control and 0-10% beta RNA','--'),('stress_inclusive','Including 20-50% stress',':')]:
        sub=cp[cp.null_range==name];ax[1].plot(sub.depth,sub.PM_sensitivity,style,marker='o',label=lab)
    ax[1].axhline(.8,color='gray',ls='--',lw=1);ax[1].set(xlabel='Simulated total count depth',ylabel='PM sensitivity at null 95th percentile',ylim=(0,1),title='B  Operating point fails');ax[1].legend(frameon=False,fontsize=7)
    for model,lab in [('alpha','Alpha control'),('mixture_0.1','10% beta RNA'),('mixture_0.5','50% beta RNA stress')]:
        sub=cr[cr.comparison=='PM_vs_'+model];ax[2].plot(pd.to_numeric(sub.depth),sub.auc,marker='o',label=lab)
    ax[2].axhline(.8,color='gray',ls='--',lw=1);ax[2].set(xlabel='Simulated total count depth',ylabel='AUROC',ylim=(.5,1),title='C  Ranking and sensitivity differ');ax[2].legend(frameon=False,fontsize=7)
    fig.tight_layout();save(fig,'figure1_calibration')
    p=pd.read_csv(EV/'primary_results.tsv',sep='\t');s=p[(p.annotation=='alpha_related')&(p.contrast=='H_GLP-D')]
    fig,ax=plt.subplots(2,1,figsize=(7.2,6.3));x=np.arange(4)
    for j,(sens,lab) in enumerate([('primary','Main alpha definition'),('strict_marker','ARX or IRX2 detected'),('GCG_ge10','GCG >=10'),('trim99','Upper 1% counts removed')]):
        z=s[s.sensitivity==sens].sort_values('donor');ax[0].plot(x,z.delta,marker=['o','s','^','D'][j],label=lab)
    z=s[s.sensitivity=='primary'].sort_values('donor');depth=pd.read_csv(EV/'depth_stability_summary.tsv',sep='\t').sort_values('donor');ax[0].errorbar(x,depth.mean_delta,yerr=depth.SD_delta,marker='x',ls='--',color='black',label='2000 counts: 20 draws, mean +/- MC SD')
    ax[0].axhline(0,color='gray',lw=1);ax[0].set(xticks=x,xticklabels=z.donor,xlabel='Biological donor',ylabel='Mean score difference',title='A  H + GLP-1 minus vehicle');ax[0].legend(frameon=False,fontsize=7)
    n=pd.read_csv(EV/'countermodel_stability_summary.tsv',sep='\t');n=n[n.contrast=='H_GLP-D']
    for j,models in enumerate([['beta_RNA_0','beta_RNA_0.01','beta_RNA_0.05','beta_RNA_0.1'],['doublet_fraction_0.01','doublet_fraction_0.05','doublet_fraction_0.1'],[m for m in n.model.unique() if m.startswith('cycle_')]]):
        if not models:continue
        vals=[n[(n.donor==d)&n.model.isin(models)].pooled_MC_q95.max() for d in z.donor]
        ax[1].plot(x,vals,marker='o',ls='--',label=['0-10% beta RNA MC envelope','1-10% doublet MC envelope','Cycling composition MC envelope'][j])
    ax[1].plot(x,z.delta,color='black',marker='D',lw=2,label='Observed main difference');ax[1].axhline(0,color='gray',lw=1)
    ax[1].set(xticks=x,xticklabels=z.donor,xlabel='Biological donor',ylabel='Score difference / simulated q95',title='B  Assumed countermodels with pool variability');ax[1].legend(frameon=False,fontsize=8)
    fig.tight_layout();save(fig,'figure2_donors')
    b=pd.read_csv(evidence_path('biological_challenge_cell_scores.tsv.gz'),sep='\t')
    fig,ax=plt.subplots(2,1,figsize=(7.2,5.8))
    for grp,lab,color in [('SMOC1_high','SMOC1 positive (219)','#b55a24'),('SMOC1_low','SMOC1 zero (2683)','#49667c')]:
        a=np.sort(b[b.SMOC1_beta_group==grp].conversion_frozen_universe)
        ax[0].step(a,np.arange(1,len(a)+1)/len(a),where='post',label=lab,color=color)
    ax[0].set(xlabel='Frozen score with declared missing-feature handling',ylabel='Cumulative fraction',title='A  Observational beta-like states');ax[0].legend(frameon=False,fontsize=8)
    cap=pd.read_csv(EV/'reference_capture_summary.tsv',sep='\t');ax[1].bar(cap.preparation,cap.median_detected_score_genes,color=[C[k] for k in cap.preparation]);ax[1].set(ylabel='Median detected score genes\n(of 100)',title='B  Reference feature capture');ax[1].set_ylim(0,20)
    fig.tight_layout();save(fig,'figure3_biological_and_capture')
    print('Wrote three figures as PNG and PDF.')

if __name__=='__main__':
    import argparse
    from public_runtime import execute_stage
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();execute_stage('figures',args.data_root,args.output_dir)
