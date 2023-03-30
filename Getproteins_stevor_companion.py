import sys
import re
from Bio import SeqIO

PfAAfasta = sys.argv[1]

Protein_sequences = SeqIO.parse(PfAAfasta,"fasta")

var_pattern = re.compile(r"stevor", re.IGNORECASE)
Nb_protein = 0
Stevor_sequences = []

for record in Protein_sequences:
    header_match = var_pattern.search(record.description)
    if header_match:
        Nb_protein += 1
        Stevor_sequences.append(record)

SeqIO.write(Stevor_sequences, "J0272815_1701_var.fa", "fasta")
print("Done ! %i protein sequences were found to have an id that matched with stevor" % Nb_protein)