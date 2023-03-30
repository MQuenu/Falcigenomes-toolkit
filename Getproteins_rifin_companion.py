import sys
import re
from Bio import SeqIO

PfAAfasta = sys.argv[1]

Protein_sequences = SeqIO.parse(PfAAfasta,"fasta")

var_pattern = re.compile(r"RIFIN", re.IGNORECASE)
Nb_protein = 0
Rif_sequences = []

for record in Protein_sequences:
    header_match = var_pattern.search(record.description)
    if header_match:
        Nb_protein += 1
        Rif_sequences.append(record)

SeqIO.write(Rif_sequences, "Rifin_proteins.fa", "fasta")
print("Done ! %i protein sequences were found to have an id that matched with RIFIN" % Nb_protein)