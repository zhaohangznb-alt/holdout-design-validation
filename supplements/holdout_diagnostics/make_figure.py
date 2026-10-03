from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();R=a.output.resolve();(R/'02_editable_source/figures').mkdir(parents=True,exist_ok=True)
T='alk braf egfr erbb2 kdr met ache dpp4 ptgs2 hdac1 pde4d'.split()
c=pd.read_csv(R/'analysis/context_diagnostics.csv',float_precision='round_trip');p=pd.read_csv(R/'analysis/target_route_diagnostics.csv',float_precision='round_trip')
p=p[p.route.eq('temporal_exact')&p.reference_route.eq('source_exact')].set_index('target').loc[T]
a=c[c.route.eq('source_exact')].set_index('target').loc[T]
b=c[c.route.eq('temporal_exact')].set_index('target').loc[T]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False,'axes.linewidth':.6,'xtick.major.size':3,'ytick.major.size':0})
fig,axes=plt.subplots(1,2,figsize=(7.05,4.7),sharey=True,gridspec_kw={'width_ratios':[1,1.22]})
y=np.array([0,1,2,3,4,5,6.6,7.6,8.6,9.6,10.6]);blue='#28749A';orange='#B75A27';ink='#292F33'
for ax in axes:
    ax.set_ylim(11.25,-.65);ax.set_yticks(y,labels=[t.upper() for t in T]);ax.axhline(5.8,color='#BCC4C9',lw=.6);ax.grid(axis='x',color='#E5E8EA',lw=.5);ax.set_axisbelow(True)
ax=axes[0]
for yy,x,z in zip(y,a.similarity_median,b.similarity_median):ax.plot([x,z],[yy,yy],color='#AAB3B8',lw=1,zorder=2)
ax.scatter(a.similarity_median,y,s=28,color=blue,label='Source exact',zorder=3)
ax.scatter(b.similarity_median,y,s=30,color=orange,marker='D',label='Year exact',zorder=3)
ax.set_xlim(.2,.8);ax.set_xticks([.2,.4,.6,.8]);ax.set_xlabel('Median nearest-training Tanimoto')
ax.legend(loc='lower left',bbox_to_anchor=(-.03,1.03),frameon=False,ncol=1,fontsize=8.5,handletextpad=.35)
ax.set_title('(a) Structural proximity',loc='left',y=1.22,fontweight='bold',fontsize=10)
ax=axes[1];ax.axvline(0,color='#68767E',lw=.8)
for yy,x,z in zip(y,p.rf_mse_log_change,p.baseline_log_change):ax.plot([x,z],[yy,yy],color='#CAD0D4',lw=1,zorder=2)
ax.scatter(p.rf_mse_log_change,y-.17,s=26,color=orange,marker='s',label='RF error',zorder=3)
ax.scatter(p.baseline_log_change,y,s=28,color=blue,label='Training-mean error',zorder=3)
ax.scatter(p.rf_relative_error_log_change,y+.17,s=25,color=ink,marker='D',label='Relative error (RF / mean)',zorder=3)
ax.set_xlim(-.52,1.32);ax.set_xticks([-.4,0,.4,.8,1.2]);ax.set_xlabel('Natural log ratio: year exact / source exact')
ax.legend(loc='lower left',bbox_to_anchor=(-.03,1.03),frameon=False,fontsize=8.5,handletextpad=.4)
ax.set_title('(b) Error components',loc='left',y=1.22,fontweight='bold',fontsize=10)
fig.subplots_adjust(left=.105,right=.987,bottom=.16,top=.73,wspace=.29)
out=R/'02_editable_source/figures/Figure_6.pdf';fig.savefig(out,metadata={'Title':'Frozen-prediction holdout diagnostics','Creator':'Matplotlib from frozen numerical records'})
fig.savefig(R/'analysis/Figure_6.png',dpi=180);plt.close(fig)
mapping=pd.DataFrame({'target':T,'source_median_similarity':a.similarity_median.values,'year_median_similarity':b.similarity_median.values,'rf_mse_log_change':p.rf_mse_log_change.values,'baseline_log_change':p.baseline_log_change.values,'relative_error_log_change':p.rf_relative_error_log_change.values})
mapping.to_csv(R/'analysis/Figure_6_data.csv',index=False)
summary=dict(contexts=44,prediction_files=220,test_context_appearances=int(c.test_n.sum()),targets=11,year_exact=dict(lower_median_similarity=int(p.similarity_change.lt(0).sum()),median_similarity_change=float(p.similarity_change.median()),smaller_training_sets=int(p.train_n_ratio.lt(1).sum()),median_training_size_ratio=float(p.train_n_ratio.median()),lower_rf_skill=int(p.rf_skill_change.lt(0).sum()),median_rf_skill_change=float(p.rf_skill_change.median()),rf_mse_increased=int(p.rf_mse_log_change.gt(0).sum()),baseline_mse_increased=int(p.baseline_log_change.gt(0).sum()),ptgs2=dict(rf_mse_ratio=float(np.exp(p.loc['ptgs2','rf_mse_log_change'])),baseline_mse_ratio=float(np.exp(p.loc['ptgs2','baseline_log_change'])),skill_change=float(p.loc['ptgs2','rf_skill_change']))))
(R/'RESULT_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))
