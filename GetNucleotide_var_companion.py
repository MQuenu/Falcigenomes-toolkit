#!/usr/bin/env python

### Import libraries

import sys
import re
import gffutils
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


### import data

Pfgenome_annotation = sys.argv[1]
Pfgenome_fasta = sys.argv[2]
PfGenome_db = "Pfgenome.db"
db = gffutils.create_db(Pfgenome_annotation, PfGenome_db, force=True, keep_order=True, merge_strategy="merge", sort_attribute_values=True)
features = db.features_of_type("polypeptide")
Genome_dict = SeqIO.to_dict(SeqIO.parse(Pfgenome_fasta,"fasta"))

### Usage

usage = 'Usage: ' + sys.argv[0] + ' <gff file> <fasta file>'

if len(sys.argv) != 3:
    print("Error in input data")
    print(usage)
    exit()

### Write a dunction to get a sequence from a particular ID
    ### Inputs: feature_id: ID of the sequence to be retrieved (will be looped through)
             #  g-dict: a genome dictionary parsed with SeqIO 
             #  dbase : the gffutils database

def get_feature_sequence(feature_id, g_dict, dbase):
    feature = dbase[feature_id]
    seq_object = g_dict[feature.seqid]
    seq = seq_object[feature.start-1:feature.end].seq
    if feature.strand == "-":
        seq = seq.reverse_complement()
    seq_record = SeqRecord(seq=seq, id=feature.id, description="")
    return seq_record

### Define and look for the regular expression pattern (either var or rifin / setvor) for the product of interest

product_pattern = re.compile(r"PfEMP1", re.IGNORECASE)
matching_genes = []
for feature in features:
    if "product" in feature.attributes and product_pattern.search(feature.attributes["product"][0]):
        matching_genes.append(feature.id)

### Create a seq_record object containing the sequences

seq_records = []
Nb_complete_genes = 0
for sequence_id in matching_genes:
    seq_record = get_feature_sequence(sequence_id, Genome_dict, db)
    if len(seq_record.seq) > 2500:
        Nb_complete_genes += 1
        seq_records.append(seq_record)

SeqIO.write(seq_records, "J0272815_1701_var.fa", "fasta")

print("There were %i polypeptides matching the expression 'PfEMP1' in the gff file" % (Nb_complete_genes))