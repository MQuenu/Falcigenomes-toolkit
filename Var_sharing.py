import subprocess
import sys
from Bio import SeqIO


### Import the two fasta files, usage

fasta_file_1 = sys.argv[1]
fasta_file_2 = sys.argv[2]

usage = 'Usage: ' + sys.argv[0] + ' <var_nucleotide_fasta_1>  <var_nucleotide_fasta_2>'

if len(sys.argv) != 3:
    print(usage)
    exit()

### Check for duplicated genes within each of the fasta files 

print("Checking for duplicated var genes within the first fasta file...")
cmd = ['cd-hit', '-i', fasta_file_1, '-o', 'fasta_1_clustered.fasta', '-c', '1']
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
if stderr:
    print(stderr.decode('utf-8'))
else:
    print('clustering step done !')

clustered_sequences_fasta_1 = SeqIO.parse('fasta_1_clustered.fasta.clstr', 'fasta')
Nb_duplicated_var_genes_1 = 0
for record in clustered_sequences_fasta_1:
    if "100.00%" in record.seq:
        Nb_duplicated_var_genes_1 += 1
print('there were %i duplicated var genes in the first fasta file' % (Nb_duplicated_var_genes_1))


print("Checking for duplicated var genes within the second fasta file...")
cmd = ['cd-hit', '-i', fasta_file_2, '-o', 'fasta_2_clustered.fasta', '-c', '1']
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
if stderr:
    print(stderr.decode('utf-8'))
else:
    print('clustering step done !')

clustered_sequences_fasta_2 = SeqIO.parse('fasta_2_clustered.fasta.clstr', 'fasta')
Nb_duplicated_var_genes_2 = 0
for record in clustered_sequences_fasta_1:
    if "100.00%" in record.seq:
        Nb_duplicated_var_genes_1 += 1
print('there were %i duplicated var genes in the second fasta file' % (Nb_duplicated_var_genes_2))

### Concatenate the fasta files from the cd-hit clustering steps

Nb_unique_var_genes_fasta1 = 0
Nb_unique_var_genes_fasta2 = 0

with open("concatenated_fastas.fa", "w") as outfile:
    for record in SeqIO.parse('fasta_1_clustered.fasta', "fasta"):
        Nb_unique_var_genes_fasta1 += 1
        SeqIO.write(record, outfile, "fasta")
    for record in SeqIO.parse('fasta_2_clustered.fasta', "fasta"):
        Nb_unique_var_genes_fasta2 += 1
        SeqIO.write(record, outfile, "fasta")

### Subprocess CD-hit to cluster sequences

print("Now cheking for shared var genes between the two fasta files...")
cmd = ['cd-hit', '-i', 'concatenated_fastas.fa', '-o', 'clustered_sequences.fasta', '-c', '0.95']
p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
stdout, stderr = p.communicate()
if stderr:
    print(stderr.decode('utf-8'))
else:
    print('CD-HIT completed successfully !')

### Parse the cd-hit output

clustered_sequences = SeqIO.parse("clustered_sequences.fasta.clstr", "fasta")

nb_shared_var_genes = 0

for record in clustered_sequences:
    if "100.00%" in record.seq:
        nb_shared_var_genes += 1

### Printing results

print("Number of unique var genes in the first fasta file: %i" % Nb_unique_var_genes_fasta1)
print("Number of unique var genes in the second fasta file: %i" % Nb_unique_var_genes_fasta2)
print("Number of shared var genes: %s" % nb_shared_var_genes)
print("details of the clustered sequences can be found in the file 'clustered_sequences.fasta.clstr'")