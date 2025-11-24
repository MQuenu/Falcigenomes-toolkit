# Falcigenomes-toolkit
A series of python and R scripts I developped to facilitate extraction of information from a set of falciparum genomes. Most scripts work from the annotation outputs of the  Companion pipeline (https://companion.gla.ac.uk/) and the outputs of the STRIDE pipeline (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04515-8).

1. GetNucleotide_var_companion.py

Extract nucleotide sequences of var genes from a genome assembly. gff and fasta provided by the companion pipeline.

Usage:
```bash
python GetNucleotide_var_companion.py <gff file> <fasta file>
```

2. Getnucleotides_rifins.py

Extract nucleotide sequences of rif and stevor genes from a set pf genome assembly and STRIDE outputs. gff and fasta provided by the companion pipeline, the file annotated.proteins.txt is the output of STRIDE. Change paths within script.

Usage:
```bash
python GetNucleotide_var_companion.py
```

3. LongReads_variant_calling.sh

A pipeline to perform core genome P. falciparum SNPs and small indels variant calling from raw pacbio hifi raw reads data. Works with nanopore data. It needs a reference 3D7 genome assembly for mapping.

Usage:
```bash
bash LongReads_variant_calling.sh <bam_reads_file>
```
