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

4. Var_sharing.py

Given two fasta files containing var nucleotide sequences from two genomes, will output the number of common shared var genes between genomes. Var genes can here be substituted by rif of stevor to get measures of rif/stevor sharing. The script relies on CD-hit (https://academic.oup.com/bioinformatics/article/22/13/1658/194225) for clustering of var genes.

Usage:
```bash
python Var_sharing.py <fasta_var1> <fasta_var_2>
```

5. extract_var_ups.sh

Given an fasta file scaffold.fa, an annotation file scafs.gff and a list of var genes IDs var_IDs.txt, the script will extract the 600 upstream sequence of each var gene and output it in a fasta file upstream_var.fa

Usage:
```bash
bash extract_var_ups.sh
```

6. Getproteins_var_companion.py, Getproteins_rifin_companion.py and Getproteins_rifin_companion.py (non-recommended)

Given a proteins.fa companion output file, will extract protein sequences of var, rifin, stevor (depending on script). Method less sensible than STRIDE for rifins and stevors.

Usage:
```bash
bash Getproteins_var_companion.py <protein fasta file>
```
