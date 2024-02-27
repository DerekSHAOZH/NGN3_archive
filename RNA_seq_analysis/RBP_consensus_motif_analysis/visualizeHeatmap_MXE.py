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
            
            input_fn = dire+"%s_background_%s/MXE/pVal.%s.vs.bg.RNAmap.txt"%(condition,control,direction)
            df = pd.read_csv(input_fn, sep = "\t")
            df=df[(df.drop(columns='RBP').min(1)<threshold)]
            df['RBP']=df["RBP"].str.split('.',expand=True)[0]
            df=df.sort_values("RBP").groupby('RBP').min()
            vals=-1*np.log10(df)

            fig = plt.figure(figsize=(12, df.shape[0]*0.3)) 
            plt.suptitle("RBP motif enriched for first exons %s at %s vs %s" %(m[direction],condition,control),fontsize=16)
            gs = gridspec.GridSpec(2,1, height_ratios=[1,df.shape[0]])
            
            ax0 = plt.subplot(gs[0])
            ax0.add_patch(Rectangle((0,0),1,1,edgecolor="black",facecolor="grey"))
            ax0.plot([1,1.9],[0.5,0.5],color="black",linewidth=2)
            ax0.plot([1.8,2.0],[0.25,0.75],color="black",linewidth=2)
            ax0.plot([2.0,2.2],[0.25,0.75],color="black",linewidth=2)
            ax0.plot([2.1,3],[0.5,0.5],color="black",linewidth=2)
            y = [0,1,1,0]
            x = [3.0,3.0,4.0,3.75]
            ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="lightgreen"))
            x = [4.0,4.25,5.0,5.0]
            ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="lightgreen"))
            ax0.plot([5,5.9],[0.5,0.5],color="black",linewidth=2)
            ax0.plot([5.8,6.0],[0.25,0.75],color="black",linewidth=2)
            ax0.plot([6.0,6.2],[0.25,0.75],color="black",linewidth=2)
            ax0.plot([6.1,7],[0.5,0.5],color="black",linewidth=2)
            
            x = [7.0,7.0,8.0,7.75]
            ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="darkgreen"))
            x = [8.0,8.25,9.0,9.0]
            ax0.add_patch(patches.Polygon(xy=list(zip(x,y)),edgecolor="black",facecolor="darkgreen"))
            ax0.plot([9,9.9],[0.5,0.5],color="black",linewidth=2)
            ax0.plot([9.8,10.0],[0.25,0.75],color="black",linewidth=2)
            ax0.plot([10.0,10.2],[0.25,0.75],color="black",linewidth=2)
            ax0.plot([10.1,11],[0.5,0.5],color="black",linewidth=2)
            
            ax0.add_patch(Rectangle((11.0,0),1,1,edgecolor="black",facecolor="grey"))
            ax0.set_xlim([0,12])
            ax0.axis('off')
            
            ax1 = plt.subplot(gs[1])
            heatmap=ax1.imshow(vals,aspect='auto',cmap='Reds',vmin=2,vmax=5)
            ax1.set_ylim(len(df.index)-0.5, -0.5)
            ax1.set_yticks(np.arange(len(df.index)))
            ax1.set_yticklabels(df.index)
            ax1.set_xticklabels([])
            ax1.set_xticks([])
            fig.tight_layout()

            plt.savefig(dire+"Plots/MXE_%s_%s_vs_%s.pdf" %(m[direction],condition,control))