"""Revision 02 presentation of retained decisions, components and fixed cohorts."""
from pathlib import Path
from data_paths import logical_evidence_key
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from revision02_common import ROOT, EV, OLD, sha

OUT=ROOT/'manuscript/figures'
DONORS=['D1','D2','D3','D4']
COLORS=['#355C7D','#C96F32','#35826A','#8D4E8A']
records={}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
                     'axes.labelsize':8,'legend.fontsize':7,'xtick.labelsize':7,'ytick.labelsize':7})


def read(name,old=False):return pd.read_csv((OLD if old else EV)/(name+'.tsv'),sep='\t')


def panel(ax,label,title):
    ax.set_title(label+'  '+title,loc='left',fontweight='bold',pad=10)
    ax.grid(axis='y',color='#E6E9EC',lw=.5,zorder=0)


def save(fig,name,sources):
    fig.savefig(OUT/(name+'.png'),dpi=220,bbox_inches='tight',facecolor='white')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight',facecolor='white',metadata={'CreationDate':None,'ModDate':None})
    records[name]={'sources':{logical_evidence_key(p):sha(p) for p in sources},
                   'png_sha256':sha(OUT/(name+'.png')),'pdf_sha256':sha(OUT/(name+'.pdf'))}
    plt.close(fig)


def mechanism():
    native=read('zero_baseline_contrasts').pivot(index='donor',columns='component',values='delta_mean').loc[DONORS]
    thin=read('thinning_stability_summary',True);thin=thin[thin.population=='alpha_broad']
    fig,axes=plt.subplots(3,2,figsize=(6.65,7.7));a,b=axes[0];x=np.arange(4)
    for offset,key,label,color in [(-.24,'occupancy_up','Up detection / 100',COLORS[2]),
                                  (0,'occupancy_down','Minus down detection / 100',COLORS[3]),
                                  (.24,'positive_rank_residual','Rank remainder','#AA9438')]:
        a.bar(x+offset,native[key],width=.22,color=color,label=label,zorder=3)
    a.scatter(x,native.score,color='#172D3D',marker='D',s=22,label='Total score',zorder=4)
    a.set(xticks=x,xticklabels=DONORS,ylim=(-.01,.043),ylabel='Treated minus vehicle')
    a.legend(frameon=False,loc='upper left',fontsize=6.8);a.axhline(0,color='gray',lw=.7)
    panel(a,'A','Detection imbalance is primary')
    for offset,key,label,color in [(-.18,'centered_up','Up, zero baseline removed',COLORS[2]),
                                  (.18,'centered_down','Down, zero baseline removed',COLORS[3])]:
        b.bar(x+offset,native[key],width=.34,color=color,label=label,zorder=3)
    b.scatter(x,native.score,color='#172D3D',marker='D',s=22,zorder=4)
    b.set(xticks=x,xticklabels=DONORS,ylim=(-.01,.043),ylabel='Treated minus vehicle')
    b.legend(frameon=False,fontsize=6.8,loc='upper left');b.axhline(0,color='gray',lw=.7)
    panel(b,'B','Separate zero-centered blocks')
    for i,(donor,ax) in enumerate(zip(DONORS,axes[1:].ravel())):
        sub=thin[thin.donor==donor].sort_values('depth')
        ax.plot(sub.depth,sub.original_same_eligible_delta,'s--',color='#63727E',ms=3.5,label='Native, same eligible cells')
        ax.errorbar(sub.depth,sub.mean_delta,yerr=sub.SD_delta,fmt='o-',color=COLORS[0],ms=3.5,capsize=2,label='Thinned mean ± technical SD')
        ax.axhline(0,color='gray',lw=.7)
        labels=[f'{int(r.depth/1000)}k\n{int(r.control_eligible_n)}:{int(r.treated_eligible_n)}' for r in sub.itertuples()]
        ax.set(xscale='log',xticks=sub.depth,xticklabels=labels,ylim=(-.01,.032),xlabel='Counts; vehicle:treated eligible n')
        ax.minorticks_off()
        if i==0:ax.legend(frameon=False,fontsize=6.4);ax.set_ylabel('Score contrast')
        panel(ax,chr(ord('C')+i),donor+': same-eligible comparison')
    fig.subplots_adjust(left=.12,right=.98,bottom=.08,top=.95,hspace=.64,wspace=.38)
    save(fig,'rev2_mechanism',[EV/'zero_baseline_contrasts.tsv',OLD/'thinning_stability_summary.tsv'])


def thresholds():
    table=read('threshold_summary');rank=table[(table.method=='signed_rank')&(table.cohort.isin(['all_native','depth_eligible']))]
    fig,axes=plt.subplots(3,2,figsize=(6.65,7.7))
    panels=[('source_fixed','alpha','A','Source threshold 0: alpha FPR',(0,3)),
            ('source_fixed','beta','B','Source threshold 0: beta detection',(50,101)),
            ('target_native_fixed','alpha','C','Native target q95 held fixed: alpha FPR',(0,65)),
            ('target_depth_recalibrated','alpha','D','Depth-recalibrated alpha FPR',(0,7))]
    for ax,(rule,state,letter,title,ylim) in zip(axes[:2].ravel(),panels):
        for donor,color in zip(range(1,5),COLORS):
            sub=rank[(rank.donor==donor)&(rank.rule==rule)&(rank.cell_type==state)&(rank.depth>0)].sort_values('depth')
            ax.errorbar(sub.depth,100*sub.positive_fraction,yerr=100*sub.positive_fraction_technical_SD,
                        color=color,fmt='o-',ms=3,capsize=2,label='human'+str(donor))
        ax.set(xscale='log',xticks=[1000,2000,5000,10000],xticklabels=['1k','2k','5k','10k'],ylim=ylim,
               ylabel='Positive fraction (%)',xlabel='Thinned counts')
        ax.minorticks_off();panel(ax,letter,title)
        if letter=='A':ax.legend(frameon=False,ncol=2,fontsize=6.5,loc='upper right')
        if letter=='D':ax.axhline(5,color='gray',ls=':',lw=.7)
    e,f=axes[2]
    for donor,color in zip(range(1,5),COLORS):
        sub=rank[(rank.donor==donor)&(rank.rule=='source_fixed')&(rank.cell_type=='alpha')&(rank.depth>0)].sort_values('depth')
        e.errorbar(sub.depth,sub.AUROC,yerr=sub.AUROC_technical_SD,fmt='o-',ms=3,color=color,capsize=2)
    e.set(xscale='log',xticks=[1000,2000,5000,10000],xticklabels=['1k','2k','5k','10k'],ylim=(.96,1.004),ylabel='Beta-versus-alpha AUROC',xlabel='Thinned counts');e.minorticks_off()
    panel(e,'E','Ranking remains strong')
    paired=rank[(rank.rule=='source_fixed')&(rank.cell_type=='beta')&(rank.depth==2000)].sort_values('donor')
    for i,row in enumerate(paired.itertuples()):
        f.plot([i-.16,i+.16],[100*row.native_positive_fraction,100*row.positive_fraction],color=COLORS[i],lw=1.5)
        f.scatter(i-.16,100*row.native_positive_fraction,marker='s',facecolors='none',edgecolors=COLORS[i],s=28,label='Native, same cells' if i==0 else None)
        f.errorbar(i+.16,100*row.positive_fraction,yerr=100*row.positive_fraction_technical_SD,color=COLORS[i],fmt='o',ms=4,capsize=2,label='2k thinned' if i==0 else None)
    f.set(xticks=range(4),xticklabels=['human1','human2','human3','human4'],ylim=(68,100),ylabel='Beta detection (%)')
    f.legend(frameon=False,fontsize=6.5,loc='upper right');panel(f,'F','Source rule: paired detection loss')
    fig.subplots_adjust(left=.12,right=.98,bottom=.07,top=.95,hspace=.68,wspace=.40)
    save(fig,'rev2_threshold_transport',[EV/'threshold_summary.tsv'])


def supplement():
    native=read('zero_baseline_contrasts').pivot(index='donor',columns='component',values='delta_mean').loc[DONORS]
    fixed=read('fixed_cohort_summary');floor=read('ucell_floor_contrasts');floor=floor[floor.program=='original'].set_index('donor').loc[DONORS]
    fig,axes=plt.subplots(2,2,figsize=(6.65,6.3));a,b,c,d=axes.ravel();x=np.arange(4)
    for offset,key,label,color in [(-.24,'up_rank_contribution','Raw up rank',COLORS[2]),
                                  (0,'zero_rank_baseline','Changing zero baseline','#AA9438'),
                                  (.24,'centered_up','Up after baseline removal',COLORS[0])]:
        a.bar(x+offset,native[key],width=.22,label=label,color=color,zorder=3)
    a.set(xticks=x,xticklabels=DONORS,ylim=(-.01,.065),ylabel='Treated minus vehicle');a.axhline(0,color='gray',lw=.7)
    a.legend(frameon=False,fontsize=6.8);panel(a,'A','Raw up rank includes moving baseline')
    for donor,color in zip(DONORS,COLORS):
        sub=fixed[fixed.donor==donor].sort_values('depth')
        b.errorbar(sub.depth,sub.delta_thinned_score,yerr=sub.delta_thinned_score_technical_SD,fmt='o-',ms=3,color=color,capsize=2,label=donor)
        b.axhline(sub.delta_native_score.iloc[0],color=color,ls=':',lw=.8)
    b.set(xscale='log',xticks=[1000,2000,5000,10000],xticklabels=['1k','2k','5k','10k'],ylabel='Score contrast',xlabel='Thinned counts',ylim=(-.012,.032));b.minorticks_off();b.axhline(0,color='gray',lw=.7)
    b.legend(frameon=False,ncol=4,loc='upper left',fontsize=6.5);panel(b,'B','Fixed cohort: native N ≥ 10k')
    c.bar(x-.18,100*floor.control_floor_fraction,width=.34,color=COLORS[0],label='Vehicle',zorder=3)
    c.bar(x+.18,100*floor.treated_floor_fraction,width=.34,color=COLORS[1],label='Treated',zorder=3)
    c.set(xticks=x,xticklabels=DONORS,ylim=(97,100.25),ylabel='Cells at clipped floor (%)');c.legend(frameon=False,fontsize=6.8);panel(c,'C','Original UCell-form clipping floor')
    for i,donor in enumerate(DONORS):
        d.scatter(i-.13,floor.loc[donor,'delta_clipped_mean'],color=COLORS[i],marker='s',s=30,label='Clipped' if i==0 else None)
        d.scatter(i+.13,floor.loc[donor,'delta_unclipped_mean'],facecolors='none',edgecolors=COLORS[i],marker='o',s=30,label='Unclipped diagnostic' if i==0 else None)
    d.axhline(0,color='gray',lw=.7);d.set(xticks=x,xticklabels=DONORS,ylim=(-.009,.022),ylabel='Treated minus vehicle');d.legend(frameon=False,fontsize=6.8);panel(d,'D','Clipping changes mean contrasts')
    fig.subplots_adjust(left=.12,right=.98,bottom=.10,top=.93,hspace=.58,wspace=.42)
    save(fig,'rev2_supplement_baseline_floor',[EV/(n+'.tsv') for n in ['zero_baseline_contrasts','fixed_cohort_summary','ucell_floor_contrasts']])


def main():
    mechanism();thresholds();supplement()
    (EV/'figure_build.json').write_text(json.dumps({'code_sha256':sha(Path(__file__)),'figures':records},indent=2)+'\n')


if __name__=='__main__':
    import argparse
    from public_runtime import execute_stage
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();execute_stage('figures',args.data_root,args.output_dir)
