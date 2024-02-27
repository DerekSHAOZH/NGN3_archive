import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as patches
from matplotlib import gridspec


threshold=0.001
days=["day"+str(i) for i in range(4)]+["day5"]
dire="/project/Neurodifferentiation_System/Analysis_NGN3/AS_rMAPS/Results/"

m={"dn":"longer" ,"up":"shorter"}
for direction in m.keys():
    for i in range(len(days)-1):
        control=days[i]
        for j in range(i+1, len(days)):
            condition=days[j]
            
            input_fn = dire+"%s_background_%s/A5SS/pVal.%s.vs.bg.RNAmap.txt"%(condition,control,direction)
            df = pd.read_csv(input_fn, sep = "\t")
            df=df[(df.drop(columns='RBP').min(1)<threshold)]
            if df.shape[0]:
                df['RBP']=df["RBP"].str.split('.',expand=True)[0]
                df=df.sort_values("RBP").groupby('RBP').min()
                vals=-1*np.log10(df)
    
                fig = plt.figure(figsize=(8, vals.shape[0]*0.3)) 
                plt.suptitle("RBP motif enriched for %s exons at %s vs %s via A5SS" %(m[direction],condition,control),fontsize=16)
                gs = gridspec.GridSpec(2,1, height_ratios=[1,df.shape[0]])
                
                ax0 = plt.subplot(gs[0])
                y = [0,1,1,0]
                x = [0.0,0.0,1.0,0.75]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="grey"))
                x = [1.0,1.25,2.0,2.0]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="grey"))
                x = [2.0,2.0,3.0,2.75]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="lightgreen"))
                x = [3.0,3.25,4.0,4.0]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="lightgreen"))
                ax0.plot([4,4.9],[0.5,0.5],color="black",linewidth=2)
                ax0.plot([4.8,5.0],[0.25,0.75],color="black",linewidth=2)
                ax0.plot([5.0,5.2],[0.25,0.75],color="black",linewidth=2)
                ax0.plot([5.1,6],[0.5,0.5],color="black",linewidth=2)
                x = [6.0,6.0,7.0,6.75]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="grey"))
                x = [7.0,7.25,8.0,8.0]
                ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="grey"))
                ax0.set_xlim([0,8])
                ax0.axis('off')
                
                ax1 = plt.subplot(gs[1])
                heatmap=ax1.imshow(vals,aspect='auto',cmap='Reds',vmin=2,vmax=5)
                ax1.set_ylim(len(vals.index)-0.5, -0.5)
                ax1.set_yticks(np.arange(len(vals.index)))
                ax1.set_yticklabels(vals.index)
                ax1.set_xticklabels([])
                ax1.set_xticks([])
                fig.tight_layout()
    
                plt.savefig(dire+"Plots/A5SS_%s_%s_vs_%s.pdf" %(m[direction],condition,control))