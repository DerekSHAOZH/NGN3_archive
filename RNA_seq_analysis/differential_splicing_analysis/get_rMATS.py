#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import os

pval=0.05
direout='/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitches_RBPs/rMATS/'
direg='/project/Neurodifferentiation_System/Analysis_NGN3/IsoSwitchAnalyser/GenesWithSwitches/'
direas='/project/Neurodifferentiation_System/Analysis_NGN3/AS_rMATS/Results/'

days=[0,1,2,3,5]

for i in range(len(days)-1):
    dayi=days[i]
    for j in range(i+1,len(days)):
        dayj=days[j]
        os.makedirs(direout+'day'+str(dayj)+'_background_day'+str(dayi),exist_ok=True)
        genes=list(pd.read_csv(direg+'GenesWithIsoSwitch_day'+str(dayi)+'_vs_day'+str(dayj)+'.txt',
                          header=None)[0])
        l=[]
        for fn in os.listdir(direas+'day'+str(dayj)+'_background_day'+str(dayi)+'/'):
            if fn.endswith('.MATS.JC.txt'):
                ase=fn.split('.')[0]
                df=pd.read_csv(direas+'day'+str(dayj)+'_background_day'+str(dayi)+'/'+fn,sep='\t')
                l=l+list(df[df['GeneID'].isin(genes)]['GeneID'])
                with open(direout+'day'+str(dayj)+'_background_day'+str(dayi)+'/'+ase+'_stats.txt','w') as f:
                    f.write('Number of genes with isoform switches: '+str(len(genes))+'\n')
                    f.write('Genes with detected AS: '+str(len(set(df[df['GeneID'].isin(genes)]['GeneID'])))+'\n')
                    f.write('Number of events within the genes: '+str(sum(df['GeneID'].isin(genes)))+'\n')
                df=df[(df['FDR']>pval) | (df['GeneID'].isin(genes))]
                df.to_csv(direout+'day'+str(dayj)+'_background_day'+str(dayi)+'/'+ase+'.MATS.JC.txt',sep='\t',index=False)
                df=df[df['GeneID'].isin(genes)]
                df.to_csv(direout+'day'+str(dayj)+'_background_day'+str(dayi)+'/OnlyIsoformSwitching_'+ase+'.MATS.JC.txt',sep='\t',index=False)
        with open(direout+'day'+str(dayj)+'_background_day'+str(dayi)+'/stats.txt','w') as f:
            f.write('Number of genes with isoform switches: '+str(len(genes))+'\n')
            f.write('Genes with detected AS: '+str(len(set(l)))+'\n')
                

