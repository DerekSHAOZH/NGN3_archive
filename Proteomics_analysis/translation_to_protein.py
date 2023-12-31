# Translation of transcript sequences to amino acid sequences

from Bio import SeqIO
import os
import pandas as pd
import sys

codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }

# Kozak scoring
def score(seq,start):
    kozak = {
        "A":[0.25,0.61,0.27,0.15,1.00,0.00,0.00,0.23],
        "C":[0.53,0.02,0.49,0.55,0.00,0.00,0.00,0.16],
        "G":[0.15,0.36,0.13,0.21,0.00,0.00,1.00,0.46],
        "T":[0.07,0.01,0.11,0.09,0.00,1.00,0.00,0.15]
    }

    score = 1.0
    for i in range(start,len(seq)):
        score *= kozak[seq[i]][i]
    return score

# Translation
def translate(seq, i):
    translating = True
    aa = ""

    while(translating):
        if (len(seq) < 3):
            translating = False
            #aa += "*"
        else:
            codon = seq[0:3]
            if (codon_table[codon] == "_"):
                translating = False
            else:
                aa += codon_table[codon]
            seq = seq[3:]
            i += 3
    return aa,i

# Translation scoring function
def translate_aa_seq_length(seq):
    i = 0
    uorf_seq = "M"
    uorf_sc = 0
    maybe = False
    uorf_st = 0
    uorf_end = 0
    translating = True
    uorfs_seq_frames = {0:"M",1:"M",2:"M"}

    while (translating):
        if (i > len(seq)):
            break
        elif (seq[i:i+3] == "ATG"):
            
            aa,end = translate(seq[i:], i)
            sc = score(seq[i-4:i+4],0)

            if (aa in uorfs_seq_frames[i % 3]):
                pass
            elif (len(aa) < 15 * 3):
                pass
            elif (len(aa) > 15 * 10):
                if (sc > uorf_sc):
                    return aa,i,end
                if (len(aa) > 15 * 20):
                    return aa,i,end
                if (maybe):
                    return uorf_seq,uorf_st,uorf_end
                else:
                    return aa,i,end
            elif (sc > uorf_sc):
                uorf_seq = aa
                uorf_sc = sc
                uorf_st = i
                uorf_end = end
                maybe = True

            if (aa not in uorfs_seq_frames[i % 3]):
                uorfs_seq_frames[i % 3] = aa

        i += 1
    return uorf_seq,uorf_st,uorf_end

def find_all_aa_seqs(seq):
    longest_aa_seq,start,end = translate_aa_seq_length(seq)
    return longest_aa_seq,start,end


transcripts_filename = sys.argv[1] #"/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS/nanopore.combined.transcript.fa"
outfile= sys.argv[2] #"/project/owlmayerTemporary/Bressin/Paper_Repository/Neurodifferentiation/Result/IsoformMS/nanopore.combined.translation.fa"
transcripts = SeqIO.index(transcripts_filename, "fasta")
ids = []
seqs = []
starts = []
ends = []
# Translate each transcript seperatly and place in their own fasta file
#transcript_stats = []
for transcript in transcripts:
    seq = str(transcripts[transcript].seq).strip()
    transcript_id = str(transcripts[transcript].id).strip()
    protein,start,end = find_all_aa_seqs(seq)

    seqs.append(protein)
    starts.append(start)
    ends.append(end)
    ids.append(transcript_id)


translation_file = open(outfile, "w+")
    
for i in range(1,len(seqs)):
    translation_file.write(">" + ids[i] + "\t" + str(starts[i]) + "\t" + str(ends[i]) + "\n" + seqs[i] + "\n")

translation_file.close()

