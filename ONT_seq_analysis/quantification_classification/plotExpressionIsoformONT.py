#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages

dire='/project/Neurodifferentiation_System/owlmayerTemporary/derek/plotExpression/'
iso_fn='/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/Results/Quantification/counts_TPM.txt'
# c=["#25224c","#c7890f","#b4242f","#444036"]
c=['#42d6a4','#ff6961', '#59adf6', '#9d94ff']

# In[2]:


iso_df=pd.read_csv(iso_fn)

days=[0,3,5]



#expression_df is expected to have columns "mean", "sem", and "max"
#gene_list is expected to have no dot
def plotIsoformExpressionONT(output_filename, iso_df, gene_list, days, colors, includeTotal = True, includeIsoform = True):
    with PdfPages(output_filename) as pdf:
        #sort by column lables so that keep mean & sem values in order
        iso_df=iso_df.sort_index(axis=1)
        for g in gene_list:
            df_temp = iso_df[iso_df['gene_name'] == g]
            if df_temp.shape[0]:
                gene_data = df_temp.filter(like='day').sum().to_frame().transpose()

                #sem: standard error of mean
                for d in days:
                    temp=df_temp.filter(like='day'+str(d))
                    df_temp['mean_'+str(d)]=temp.mean(1)
                    df_temp['sem_'+str(d)]=temp.std(1)/np.sqrt(temp.shape[1])

                    temp=gene_data.filter(like='day'+str(d))
                    gene_data['mean_'+str(d)]=temp.mean(1)
                    gene_data['sem_'+str(d)]=temp.std(1)/np.sqrt(temp.shape[1])

                df_temp['max']=df_temp.filter(like='day').max(axis=1)

                iso_data = df_temp
                iso_data=iso_data.sort_values("max",ascending=False)

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
                            label="ONT-seq"
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
                plt.title(iso_data.iloc[0]['gene_name'])
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

            else:
                print(f'{g} NOT detected in ONT')
# In[ ]:


#Plot for isoform switches of specified genes
plt.ioff()

####### TO-EDIT #######
#input any gene symbols of interest
genes = ['YDJC']
# genes=list(pd.read_csv('/project/Neurodifferentiation_System/owlmayerTemporary/derek/plotExpression/List/SE_day1_vs_day0.tsv',header=None)[0])
#######################



#convert gene symbol to gene id
# genes = [gene_df.loc[gene_df["gene_symbol"] == x]["gene_id"].unique()[0] for x in gene_symbol]



####### TO-EDIT #######
#specify whether to include Total, Isoform, or both
#specify output filename
plotIsoformExpressionONT(dire+'Plots/'+"YDJC.ONTseq"+'.pdf', iso_df, genes, days, c, includeTotal = True, includeIsoform = True)
#######################
