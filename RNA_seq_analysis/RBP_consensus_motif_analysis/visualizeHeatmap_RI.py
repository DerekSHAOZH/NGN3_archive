import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib import gridspec


threshold=0.001
days=["day"+str(i) for i in range(4)]+["day5"]
dire="/project/Neurodifferentiation_System/Analysis_NGN3/AS_rMAPS/Results/"

m={"dn":"retained","up":"skipped"}
for direction in m.keys():
    for i in range(len(days)-1):
        control=days[i]
        for j in range(i+1, len(days)):
            condition=days[j]
            
            input_fn = dire+"%s_background_%s/RI/pVal.%s.vs.bg.RNAmap.txt"%(condition,control,direction)
            df = pd.read_csv(input_fn, sep = "\t")
            df=df[(df.drop(columns='RBP').min(1)<threshold)]
            if df.shape[0]:
                df['RBP']=df["RBP"].str.split('.',expand=True)[0]
                df=df.sort_values("RBP").groupby('RBP').min()
                vals=-1*np.log10(df)
                
                r=max(vals.shape[0]*0.5,3)
                fig = plt.figure(figsize=(7, r)) 
                plt.suptitle("RBP motif enriched for introns %s at %s vs %s" %(m[direction],condition,control),fontsize=16)
                gs = gridspec.GridSpec(2,1, height_ratios=[1,vals.shape[0]])
                
                ax0 = plt.subplot(gs[0])
                y = [0,1,1,0]
                x = [0.0,0.0,1.0,0.75]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="grey"))
                x = [1.0,1.25,2.0,2.0]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="grey"))
                ax0.plot([2,3],[0.5,0.5],color="lightgreen",linewidth=2)
                x = [3.0,3.0,4.0,3.75]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="grey"))
                x = [4.0,4.25,5.0,5.0]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="grey"))
                ax0.set_xlim([0,5])
                ax0.axis('off')
                
                ax1 = plt.subplot(gs[1])
                heatmap=ax1.imshow(vals,aspect='auto',cmap='Reds',vmin=2,vmax=5)
                ax1.set_ylim(len(vals.index)-0.5, -0.5)
                ax1.set_yticks(np.arange(len(vals.index)))
                ax1.set_yticklabels(vals.index)
                ax1.set_xticklabels([])
                ax1.set_xticks([])
                fig.tight_layout()
    
                plt.savefig(dire+"Plots/RI_%s_%s_vs_%s.pdf" %(m[direction],condition,control))