#/bin/python

import sys
from Bio import SeqIO
import re

#### functions ###

def extract_chr_number(input_string):
    pattern = r'[^_]+_(.*)$'
    match = re.search(pattern, input_string)
    if match:
        extracted_string = match.group(1)
    return(extracted_string)

#### import files ###

Usage = "Usage: " + sys.argv[0] + "fasta_file" + "chromosome to extract (format: 01 for chromosome 1)"

if len(sys.argv) != 3:
    print(Usage)
    exit()

Input_fasta = sys.argv[1]
Chromosome_target = str(sys.argv[2])

#### parse through the fasta headers ####

fasta = SeqIO.parse(Input_fasta, "fasta")

target_sequence = []

for record in fasta:
    header = record.id
    chr_number = extract_chr_number(header)
    print(chr_number)
    if chr_number == Chromosome_target:
        target_sequence.append(record)

if target_sequence == []:
    print('error in the extraction process, please make sure the chromosome headers are of the forme : NameSample_02')
    exit()

output_name = "%s_chromosome%s.fasta" % (Input_fasta, Chromosome_target)

SeqIO.write(target_sequence, output_name, "fasta")

