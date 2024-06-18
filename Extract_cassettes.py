#!/usr/bin/env python3
import csv
from Bio import SeqIO

### Input the domains.tsv file and extract cassette coordinates

Domains_input_table = "domains_cassette.tsv"

with open(Domains_input_table, 'r') as input_tsv:
    dictionary_start_seq = {}
    dictionary_end_seq = {}
    for line in input_tsv:
        columns = line.split('\t')
        ID = columns[0]
        Domain = columns[1]
        if Domain == 'DBLa':
            seq_start_p = int(columns[3])
            dictionary_start_seq[ID] = [seq_start_p]
        elif Domain == 'CIDRa':
            seq_end_p = int(columns[4])
            dictionary_end_seq[ID] = [seq_end_p]

merged_dictionary = {}
for key in set(dictionary_start_seq.keys()).union(dictionary_end_seq.keys()):
    start_seq = dictionary_start_seq.get(key, [None])
    end_seq = dictionary_end_seq.get(key, [None])
    merged_dictionary[key] = start_seq + end_seq

## output a .tsv file with the start and end locations of the cassette (optional here)

output_tsv_file = "domain_locations_cassette.tsv"

with open(output_tsv_file, 'w', newline='') as output_tsv:
    tsv_writer = csv.writer(output_tsv, delimiter='\t')
    tsv_writer.writerow(['ID', 'Start Sequence', 'End Sequence'])
    for key, values in merged_dictionary.items():
        tsv_writer.writerow([key] + values)

### Extract the sequences from the protein sequences

fasta_input = 'Var_proteins_all_samples.fa'
var_fasta = SeqIO.parse(fasta_input,'fasta')

list_dbla = []

for record in var_fasta:
    if record.id in merged_dictionary:
        original_sequence = record.seq
        start, stop = merged_dictionary[record.id]
        if start == "" or start is None:
            continue
        if stop == "" or stop is None:
            continue
        new_sequence = original_sequence[int(start-1):int(stop)] 
        record.seq = new_sequence
        list_dbla.append(record)

SeqIO.write(list_dbla,'Cassettes_DBLa_CIDRa_sequences.fa','fasta')
