#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages

dire='/project/Neurodifferentiation_System/owlmayerTemporary/derek/plotExpression/'
iso_fn='/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Tables/RSEM_isoform_TPM_shortreads_byDay.csv'
gene_fn='/project/Neurodifferentiation_System/Analysis_NGN3/RSEM/Tables/RSEM_gene_TPM_shortreads_byDay.csv'
# c=["#25224c","#c7890f","#b4242f","#444036"]
c=['#ffb480', '#f8f38d', '#08cad1', '#c780e8']

# In[2]:


iso_df=pd.read_csv(iso_fn)
gene_df=pd.read_csv(gene_fn)

days=[0,1,2,3,5]
#sem: standard error of mean
for d in days:
    temp=iso_df.filter(like='D'+str(d))
    iso_df['mean_'+str(d)]=temp.mean(1)
    iso_df['sem_'+str(d)]=temp.std(1)/np.sqrt(temp.shape[1])
    
    temp=gene_df.filter(like='D'+str(d))
    gene_df['mean_'+str(d)]=temp.mean(1)
    gene_df['sem_'+str(d)]=temp.std(1)/np.sqrt(temp.shape[1])
    
    
    
iso_df['max']=iso_df.filter(like='D').max(axis=1)


# In[ ]:


#expression_df is expected to have columns "mean", "sem", and "max"
#gene_list is expected to have no dot
def plotIsoformExpression(output_filename, gene_df, iso_df, gene_list, days, colors, includeTotal = True, includeIsoform = True):
    with PdfPages(output_filename) as pdf:
        #sort by column lables so that keep mean & sem values in order
        gene_df=gene_df.sort_index(axis=1)
        iso_df=iso_df.sort_index(axis=1)
        for g in gene_list:
            gene_data=gene_df[gene_df['gene_id']==g]
            iso_data=iso_df[iso_df['gene_id']==g]
            iso_data=iso_data.sort_values("max",ascending=False)
            if gene_data.shape[0]:
                fig, ax = plt.subplots(figsize=(6,4))
                cur_axes = plt.gca()
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)

                if includeTotal:
                    #total expression
                    plt.errorbar(days,gene_data.iloc[0].filter(like='mean_').values,
                            yerr=gene_data.iloc[0].filter(like='sem_').values,
                            color='#22A39F',linewidth=2, 
                            # label="Total"
                            label="RNA-seq"
                            )

                
                if includeIsoform:
                    #max. number of isoform plotted = 4
                    for i in range(min(iso_data.shape[0],len(colors))):
                        plt.errorbar(days,iso_data.iloc[i].filter(like='mean_').values,
                                yerr=iso_data.iloc[i].filter(like='sem_').values,
                                color=colors[i],linewidth=2,label=iso_data.iloc[i]['transcript_id'])
                                      
                        
                plt.xticks(days,[str(d) for d in days])
                plt.xlabel('Day')
                plt.ylabel('Expression [TPM]')
                plt.title(iso_data.iloc[0]['gene_symbol'])
                a=plt.axis()
                yl=plt.yticks()
                plt.xlim([-0.3,5.3])
                plt.ylim([max(-5,-0.01*a[3]),a[3]])  
                yl=plt.yticks()
                ax.spines['left'].set_bounds(max(-5,-0.01*a[3]),yl[0][yl[0]<a[3]][-1])
                ax.spines['bottom'].set_bounds(-0.3,5)
                plt.legend()
                plt.tight_layout()
                pdf.savefig()


# In[ ]:


#Plot for isoform switches of specified genes
plt.ioff()

####### TO-EDIT #######
#input any gene symbols of interest
# gene_symbol = ["RTCA","FNBP1", "ATPAF1", "EMC3", "GPANK1"]
gene_symbol = ["ACTB","RPL41", "PTMA", "RPS15", "PPIA"]
# gene_symbol = ["DNMT3B","FLNA", "HNRNPD"]
# gene_symbol=list(pd.read_csv('/project/Neurodifferentiation_System/owlmayerTemporary/derek/plotExpression/List/SE_day3_vs_day2.tsv',header=None)[0])
#######################



#convert gene symbol to gene id
genes = []
for x in gene_symbol:
    if len(gene_df.loc[gene_df["gene_symbol"] == x]["gene_id"].unique()):
        genes.append(gene_df.loc[gene_df["gene_symbol"] == x]["gene_id"].unique()[0])


# genes = [gene_df.loc[gene_df["gene_symbol"] == x]["gene_id"].unique()[0] for x in gene_symbol]



####### TO-EDIT #######
#specify whether to include Total, Isoform, or both
#specify output filename
plotIsoformExpression(dire+'Plots/'+"DTE_but_no_DTU.RNAseq"+'.pdf', gene_df, iso_df, genes, days, c, includeTotal = True, includeIsoform = True)
#######################
