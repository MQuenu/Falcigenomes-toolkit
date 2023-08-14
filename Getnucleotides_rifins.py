#!/usr/bin/env python

import sys
from Bio import SeqIO
import gffutils
import pandas as pd

### functions

def write_fasta_from_dict(sequence_dict, output_file):
    with open(output_file, 'w') as f:
        for sequence_id, sequence in sequence_dict.items():
            f.write(">{0}\n{1}\n".format(sequence_id, sequence))

### Import data 

genome_fasta = "scafs.fasta"
genome_annotation = "scaffold.out.gff3"
table_rifins_stevors = "STRIDE/annotated.proteins.txt"

### Extract the protein IDs of rifins and stevor from the stride output, put them in two lists

df_Rifins_Stevor = pd.read_csv(table_rifins_stevors, delimiter='\t', header = None)

df_Rifins = df_Rifins_Stevor[df_Rifins_Stevor[1].str.match(r'RIFIN-[ab]', case=False)]
list_Rifins_with_suffixes = df_Rifins[0].to_list()
list_Rifins = []
for record in list_Rifins_with_suffixes:
    simplified_record = record.rsplit('.', 1)[0]
    list_Rifins.append(simplified_record)

df_Stevor= df_Rifins_Stevor[df_Rifins_Stevor[1].str.match(r'STEVOR', case=False)]
list_Stevor_with_suffixes = df_Stevor[0].to_list()
list_Stevor = []
for record in list_Stevor_with_suffixes:
    simplified_record = record.rsplit('.', 1)[0]
    list_Stevor.append(simplified_record)

### Create a gffutils database and convert the fasta file to a dictionary

db = gffutils.create_db(genome_annotation,"gffutils_db", force = True, merge_strategy="warning", disable_infer_transcripts = True)
genome_sequences = {}
with open(genome_fasta, 'r') as fasta_file:
    for record in SeqIO.parse(fasta_file, 'fasta'):
        genome_sequences[record.id] = str(record.seq)

### Extract the rifins 

### here I have to check if the gene is indeed in the gffutils db. 
# It seems STRIDE can calls rifins that are pseudogenes, and pseudogenes from companion are not imported the gffutils db

all_ids_db = []  
for feature in db.all_features():
    if 'ID' in feature.attributes:
        all_ids_db.append(feature['ID'][0])

dict_rifins = {}
count_rifins = 0

for gene_id in list_Rifins:
    if gene_id in all_ids_db:
        gene = db[gene_id]
        seq_id = gene.seqid
        start = gene.start - 1  
        end = gene.end
        sequence = genome_sequences[seq_id][start:end]
        dict_rifins[gene_id] = sequence
        count_rifins +=1
    else:
        next

print("Rifins extraction done, %i rifin genes have been found" % count_rifins)
write_fasta_from_dict(dict_rifins, "rifins.fasta")

### Extract the STEVORs

dict_stevors = {}
count_stevor = 0

for gene_id in list_Stevor:
    if gene_id in all_ids_db:
        gene = db[gene_id]
        seq_id = gene.seqid
        start = gene.start - 1  
        end = gene.end
        sequence = genome_sequences[seq_id][start:end]
        dict_stevors[gene_id] = sequence
        count_stevor += 1
    else:
        next

print("Stevors extraction done, %i stevors genes have been found" % count_stevor)
write_fasta_from_dict(dict_stevors, "stevors.fasta")