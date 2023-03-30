import sys
import re
from Bio import SeqIO

### Input the data

PfAAfasta = sys.argv[1]

Protein_sequences = SeqIO.parse(PfAAfasta,"fasta")

### Usage

usage = 'Usage: ' + sys.argv[0] + ' <protein fasta file>'

if len(sys.argv) != 2:
    print("Error in input data")
    print(usage)
    exit()

### loop through the record for 'PfEMP1' ID using regular expression

var_pattern = re.compile(r"PfEMP1", re.IGNORECASE)
Nb_protein = 0
Var_sequences = []

for record in Protein_sequences:
    header_match = var_pattern.search(record.description)
    if header_match:
        Nb_protein += 1
        Var_sequences.append(record)

### Output

SeqIO.write(Var_sequences,"Var_potential_proteins.fa","fasta")
print("Done ! %i protein sequences were found to have an id that matched with PfEMP1" % Nb_protein)