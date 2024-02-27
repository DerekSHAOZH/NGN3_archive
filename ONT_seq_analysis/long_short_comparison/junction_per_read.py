#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pysam
from collections import Counter
import pandas as pd
import numpy as np
import os
import pyranges as pr

###Load protein-coding multi-exonic genes
gtf_fn = '/project/Neurodifferentiation_System/GeneralResources/gencode.v32.primary_assembly.annotation.sorted.gtf'
gtf = pr.read_gtf(gtf_fn, as_df = True)

curr_gene = ''
goi = list()
curr_num_exon = 0
for i in range(gtf.shape[0]):
    if gtf.iloc[i,2] == 'gene' and gtf.iloc[i,9] == 'protein_coding':
        if curr_num_exon > 1:
            goi.append(curr_gene)
        curr_gene = gtf.iloc[i,8]
        curr_num_exon = 0
    if gtf.iloc[i,2] == 'exon' and gtf.iloc[i,8] == curr_gene:
        curr_num_exon += 1


genic_gtf = gtf[gtf['Feature'] == 'gene']
genic_gtf = genic_gtf[genic_gtf['gene_id'].isin(goi)]
genic_gtf = genic_gtf.iloc[:, 0:14]


###Get bam file names
ont_alignment_dir = '/project/Neurodifferentiation_System/Analysis_NGN3_ONT/Nanopore/IGV/ONT_genome_alignment/'
ont_alignment = list()

for file in os.listdir(ont_alignment_dir):
    if file.endswith(".bam"):
        ont_alignment.append(os.path.join(ont_alignment_dir, file))
ont_alignment.sort()

rna_alignment_dir = '/project/owlmayerTemporary/Derek/RNA_alignGenome/alignment/'
rna_alignment = list()
for root, dirs, files in os.walk(rna_alignment_dir):
    for file in files:
        if file.endswith(".bam"):
             rna_alignment.append(os.path.join(root, file))

rna_alignment.sort()
rna_alignment.pop(0)

alignment = ont_alignment + rna_alignment

###Define the function to compute the numebr of exon-exon junctions covered by one read based on CIGAR string
def compute_num_junction_per_read(cigar_string):
    
    #m: number of exon-exon splice junctions
    m = cigar_string.count('N')
        
    return m


###Calculate for each bam file
max_junction = 11 #can be changed
nums = list()
df_dic = {}
# whole_cigar = ''
for fn in alignment:
    bamfile = pysam.AlignmentFile(fn, "rb") 
    for i in range(genic_gtf.shape[0]):
    # for read in bamfile.fetch('chr3', 149964904, 149970895): #PFN2, at most 2 exon-exon junctions
    # for read in bamfile.fetch('chr14', 69767112, 69772005): #SRSF5, at most 8 exon-exon junctions
        for read in bamfile.fetch(genic_gtf.iloc[i,0], genic_gtf.iloc[i,3], genic_gtf.iloc[i,4]):  #whole genome
        #     print(read.mapq)
        #     print(read)
        #     break
    #         name = read.query_name
    #         names.append(name)
    #         mq = read.mapq
    #         mqs.append(mq)
            num = compute_num_junction_per_read(read.cigarstring)
            nums.append(num)
    #         cigar = read.cigarstring
    #         whole_cigar = whole_cigar + cigar
    #     print(Counter(nums))
    tmp_counter = Counter(nums)
    # print(tmp_counter)
    tmp_num = list()
    for i in range(max_junction):   
        tmp_num.append(tmp_counter[i])
    df_dic[fn] = tmp_num    
    nums = list()

df = pd.DataFrame(df_dic).transpose()
per_df = df.div(df.sum(axis=1), axis=0)

###Calculate summary statistics for two types of sequencing
ont_per_df = per_df.iloc[0:9, ].transpose()
ont_per_df['mean']=ont_per_df.mean(1)
ont_per_df['sem']=ont_per_df.std(1)/np.sqrt(ont_per_df.shape[1])

rna_per_df = per_df.iloc[9:per_df.shape[0], ].transpose()
rna_per_df['mean']=rna_per_df.mean(1)
rna_per_df['sem']=rna_per_df.std(1)/np.sqrt(rna_per_df.shape[1])

summary_dic = {'group': ['long read']*max_junction + ['short_read']*max_junction, 
              'num': list(range(max_junction))*2,
              'mean': [0]*2*max_junction, 
              'sem': [0]*2*max_junction}
summary = pd.DataFrame(summary_dic)

summary.iloc[0:max_junction, 2:summary.shape[1]] = ont_per_df.loc[:, ['mean', 'sem']]
summary.iloc[max_junction:max_junction*2, 2:summary.shape[1]] = rna_per_df.loc[:, ['mean', 'sem']]


output_fn = '/project/Neurodifferentiation_System/owlmayerTemporary/derek/ONT/Results/perc_read_cover_junction.gene.csv'
summary.to_csv(output_fn)