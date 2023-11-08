#!/bin/bash

### Extract the var genes annotations from the complete gff file

grep -wFf var_IDs.txt scaffold.out.gff3 > filtered_var_genes_step1.gff
grep "gene" filtered_var_genes_step1.gff > filtered_var_genes.gff
rm filtered_var_genes_step1.gff

## Create separate files for the different orientations

awk '$7 == "+"' filtered_var_genes.gff > forward_var.gff
awk '$7 == "-"' filtered_var_genes.gff > reverse_var.gff

### Index the genome

samtools faidx scafs.fasta

### Get the upstream region locations in a gff file (and tidy up the file)

bedtools flank  -i forward_var.gff -g scafs.fasta.fai -l 600 -r 0 > upstream_forward_raw.gff
bedtools flank  -i reverse_var.gff -g scafs.fasta.fai -l 0 -r 600 > upstream_reverse_raw.gff
awk 'BEGIN {OFS = FS = "\t"} {temp = $3; $3 = $9; $9 = temp; print}' upstream_forward_raw.gff > upstream_fvar.gff
awk 'BEGIN {OFS = FS = "\t"} {temp = $3; $3 = $9; $9 = temp; print}' upstream_reverse_raw.gff > upstream_rvar.gff
rm upstream_*raw.gff

### Get the nucleotide 

bedtools getfasta -fi scafs.fasta -bed upstream_fvar.gff -nameOnly -fullHeader > upstream_fvar.fa
bedtools getfasta -fi scafs.fasta -bed upstream_rvar.gff -nameOnly -fullHeader > upstream_rvar.fa

### reverse complement the upstream regions on the reverse strand

python Reverse_complement.py

### concatenate all records in a single fasta file

cat upstream_fvar.fa upstream_rvar.fa > upstream_var.fa

### Put all the files in a new directory for tidyness

mkdir upstreams_var/
mv upstream*fa upstreams_var/
mv upstream*gff upstreams_var/
mv reverse_var.gff upstreams_var/
mv forward_var.gff upstreams_var/