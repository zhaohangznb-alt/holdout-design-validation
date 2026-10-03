"""Scientific result plots; no fitted models or labels changed."""
from pathlib import Path
import argparse,csv,json,sys,os
HERE=Path(__file__).resolve().parent
if (HERE/'plot_dependencies').exists():sys.path.insert(0,str(HERE/'plot_dependencies'))
os.environ['MPLCONFIGDIR']=str(HERE/'plot_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
MODELS=['rf','ridge','tanimoto1nn'];NAMES=['Random forest','Ridge (alpha = 1)','Tanimoto 1-NN']
TARGETS={'kinase':['alk','braf','egfr','erbb2','kdr','met'],'nonkinase':['ache','dpp4','hdac1','pde4d','ptgs2']}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--analysis',type=Path,default=HERE/'analysis_v1');ap.add_argument('--output',type=Path,default=HERE/'figures');args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.analysis/'target_route_contrasts.csv').open(encoding='utf8',newline='') as f:rows=list(csv.DictReader(f))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':11,'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white'})
    for scenario,label in [('revised','Broad'),('single_protein','Operational')]:
        fig,axs=plt.subplots(2,3,figsize=(10.8,6.3),gridspec_kw={'height_ratios':[6,5]})
        for j,model in enumerate(MODELS):
            subset=[r for r in rows if r['view']=='main' and r['scenario']==scenario and r['model']==model]
            allvals=[float(r[key]) for r in subset for key in ['temporal_exact_minus_source_exact','source_scaffold_minus_source_exact']]
            left=min(allvals+[0]);right=max(allvals+[0]);pad=.10*(right-left+.1)
            for i,panel in enumerate(['kinase','nonkinase']):
                ax=axs[i,j];ts=TARGETS[panel];y=np.arange(len(ts));data={r['target']:r for r in subset if r['panel']==panel}
                a=[float(data[t]['temporal_exact_minus_source_exact']) for t in ts];b=[float(data[t]['source_scaffold_minus_source_exact']) for t in ts]
                ax.axvline(0,color='#555555',lw=.8,ls='--',zorder=0)
                for yy,x1,x2 in zip(y,a,b):ax.plot([x1,x2],[yy,yy],color='#c3c9ce',lw=1,zorder=1)
                ax.scatter(a,y,s=35,color='#2166ac',marker='o',label='Year exact minus source exact',zorder=3)
                ax.scatter(b,y,s=31,color='#d9822b',marker='s',label='Source scaffold minus source exact',zorder=2)
                ax.set_yticks(y,[t.upper() for t in ts]);ax.invert_yaxis();ax.set_xlim(left-pad,right+pad)
                ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0,pad=5);ax.spines['bottom'].set_color('#aaaaaa')
                ax.grid(axis='x',color='#e8e8e8',lw=.6);ax.set_axisbelow(True)
                if i==0:ax.set_title(NAMES[j],loc='left',pad=12,fontweight='bold')
                if j==0:ax.set_ylabel('Six kinase targets' if i==0 else 'Five non-kinase targets',labelpad=15)
                if i==1:ax.set_xlabel('Paired change in skill')
        handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.53,.01),fontsize=9)
        fig.suptitle(label+' cohorts: per-target changes under fixed model specifications',x=.06,y=.995,ha='left',fontsize=13,fontweight='bold')
        fig.subplots_adjust(left=.11,right=.985,top=.89,bottom=.14,wspace=.29,hspace=.25)
        for ext in ['png','pdf','svg']:fig.savefig(args.output/('paired_changes_'+scenario+'.'+ext),dpi=200)
        plt.close(fig)
    (args.output/'FIGURE_NOTES.md').write_text('Each point is a within-target route contrast. RF uses the median of three seed-level skills; deterministic models occur once. Targets have equal weight in reported summaries. Axis limits differ between model columns and are shared between target panels within a column. These are descriptive route contrasts; no route confidence intervals are plotted. The six original kinase targets and five independent non-kinase target identities are shown separately.\n',encoding='utf8')
    (args.output/'plot_environment.json').write_text(json.dumps(dict(python=sys.version,matplotlib=matplotlib.__version__,training_or_scoring_executed=False),indent=2),encoding='utf8')
if __name__=='__main__':main()
