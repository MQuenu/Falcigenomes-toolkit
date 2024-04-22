#!/bin/bash 

## input the data and create working directory

input_file=$1
sample_name=$(basename "$input_file" .bam)

if [ -e $sample_name ]; then
  echo "Sample allready processed, exiting the variant calling script"
  exit
else
  echo "processing $sample_name"
fi

mkdir "$sample_name"
cp 3D7_files/* "$sample_name"
cp bam_reads_files/$sample_name.bam "$sample_name" 
cd "$sample_name"

## convert the bam file to fastq 

samtools fastq $sample_name.bam > $sample_name.fastq
rm DC01m1.bam

## map the reads to the reference genome, convert and sort to .bam

minimap2 -a 3D7_reference.fasta $sample_name.fastq > alignment_$sample_name.sam
samtools view -bS alignment_$sample_name.sam > alignment_$sample_name.bam
rm alignment_$sample_name.sam
samtools sort alignment_$sample_name.bam -o sorted_$sample_name.bam
samtools index sorted_$sample_name.bam
rm alignment_$sample_name.bam
samtools addreplacerg -r '@RG\tID:samplename\tSM:samplename' sorted_$sample_name.bam -o $sample_name.bam ## I had to add this line otherwise gatk was bugging
rm sorted_$sample_name.bam
samtools index $sample_name.bam

## do the SNP calling using GATK

gatk HaplotypeCaller \
    -R 3D7_reference.fasta \
    -I $sample_name.bam \
    -O $sample_name.vcf 

cd ../  