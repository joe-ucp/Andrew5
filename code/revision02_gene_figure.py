"""Two paired heatmaps of retained gene results; no fitting or new inference."""
from pathlib import Path
from data_paths import logical_evidence_key
import argparse
import datetime
import hashlib
import json
import os
import textwrap
import numpy as np
import pandas as pd
from data_paths import data_path,resolve_data_root,ENV_NAME

ROOT=Path(__file__).resolve().parents[1]
EV=ROOT/'evidence/revision_02'
FIG=ROOT/'manuscript/figures'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    os.environ['MPLCONFIGDIR']=str(data_path('revision_02/gene_figure_runtime',for_write=True))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'ps.fonttype':42})
    summary_path=EV/'gene_contribution_summary.tsv';contrast_path=EV/'gene_donor_contrasts.tsv'
    summary=pd.read_csv(summary_path,sep='\t');contrasts=pd.read_csv(contrast_path,sep='\t')
    rows=[]
    for block in ['up','down']:
        chosen=summary[(summary.block==block)&summary.display_top10].sort_values('within_block_display_rank')
        if len(chosen)!=10 or chosen.within_block_display_rank.tolist()!=list(range(1,11)):
            raise ValueError('Exactly the retained first10 ranks of each block are required')
        rows.extend(chosen.gene.tolist())
    donor_order=['D1','D2','D3','D4']
    occupancy=contrasts.pivot(index='gene',columns='donor',values='delta_signed_occupancy').loc[rows,donor_order].to_numpy()*1000
    positive=contrasts.pivot(index='gene',columns='donor',values='delta_positive_mean_log1p_CP10K').loc[rows,donor_order].to_numpy()
    fig,axes=plt.subplots(1,2,figsize=(6.65,7.6),gridspec_kw={'left':.18,'right':.97,'bottom':.18,'top':.89,'wspace':.30})
    titles=['A  Signed occupancy contribution\nScore units × 1,000','B  Conditional-positive expression\nΔ mean log1p(CP10K)']
    for k,(ax,matrix,title) in enumerate(zip(axes,[occupancy,positive],titles)):
        limit=float(np.nanmax(np.abs(matrix))) if np.isfinite(matrix).any() else 1.
        if not limit:limit=1.
        norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit)
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#ECECEC')
        image=ax.imshow(matrix,aspect='auto',cmap=cmap,norm=norm,interpolation='nearest')
        ax.set_xticks(range(4),donor_order,fontsize=8.5)
        ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=4)
        ax.set_yticks(range(20),rows if k==0 else ['']*20,fontsize=8.3)
        ax.tick_params(axis='y',length=0,pad=5)
        ax.set_title(title,loc='left',fontweight='bold',fontsize=9,pad=24)
        ax.set_xticks(np.arange(-.5,4,1),minor=True);ax.set_yticks(np.arange(-.5,20,1),minor=True)
        ax.grid(which='minor',color='white',linewidth=.7);ax.tick_params(which='minor',length=0)
        ax.axhline(9.5,color='#263744',linewidth=1.6)
        for spine in ax.spines.values():spine.set_visible(False)
        for i in range(20):
            for j in range(4):
                value=matrix[i,j]
                label='NA' if not np.isfinite(value) else f'{value:+.2f}'
                color='white' if np.isfinite(value) and abs(value)>.61*limit else '#18232D'
                ax.text(j,i,label,ha='center',va='center',fontsize=7.4,color=color)
        position=ax.get_position();color_ax=fig.add_axes([position.x0,.113,position.width,.018])
        colorbar=fig.colorbar(image,cax=color_ax,orientation='horizontal')
        colorbar.ax.tick_params(labelsize=7,length=2,pad=2)
        colorbar.set_label('Signed score contribution' if k==0 else 'Conditional-positive expression change',fontsize=7.2,labelpad=2)
    axes[0].text(-.43,.755,'UP BLOCK',transform=axes[0].transAxes,rotation=90,
                 ha='center',va='center',fontsize=8,fontweight='bold',color='#35826A')
    axes[0].text(-.43,.245,'DOWN BLOCK',transform=axes[0].transAxes,rotation=90,
                 ha='center',va='center',fontsize=8,fontweight='bold',color='#8D4E8A')
    note=('H_GLP minus vehicle. The same top 10 genes per block are shown for every donor, ranked by equal-donor '
          'mean absolute occupancy contribution. Positive-expression means condition separately in each arm. '
          'All 100 genes are retained in the source tables; no pathway or lineage inference is implied.')
    fig.text(.02,.019,textwrap.fill(note,123),ha='left',va='bottom',fontsize=6.8,color='#435460')
    FIG.mkdir(parents=True,exist_ok=True)
    png=FIG/'rev2_gene_drivers.png';pdf=FIG/'rev2_gene_drivers.pdf'
    fig.savefig(png,dpi=240,facecolor='white',metadata={'Software':'STRIKE007 revision02 gene figure'})
    fig.savefig(pdf,facecolor='white',metadata={'Creator':'STRIKE007 revision02 gene figure','CreationDate':None,'ModDate':None})
    plt.close(fig)
    status={'status':'RENDERED_VISUAL_REVIEW_PENDING','updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'figure_inches':[6.65,7.6],'plotted_genes':rows,'donors':donor_order,
            'sources':{logical_evidence_key(p):sha(p) for p in [summary_path,contrast_path]},
            'code_sha256':sha(Path(__file__)),'png_sha256':sha(png),'pdf_sha256':sha(pdf),
            'panel_A':'signed occupancy donor contrast multiplied by1000; equal50 blocks coefficient±1/100',
            'panel_B':'conditional-positive mean log1p(CP10K) donor contrast for the same fixed gene rows',
            'visual_review':'PENDING'}
    (EV/'gene_figure_status.json').write_text(json.dumps(status,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(status,indent=2))


if __name__=='__main__':
    import argparse
    from public_runtime import execute_stage
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();execute_stage('figures',args.data_root,args.output_dir)
