#!/bin/bash

# yeast pooled segregant code analysis

# this script is not intended to be run as is
# it provides a synopsis of the parts you will need to build your pipeline scripts
# please mine this script for useful bits as suits your needs
# this code is not complete and you will need to troubleshoot to get yours running

# the expected input is:
#   - Illumina paired end reads, typically 2 x 150
#   - two input samples, typically:
#       - one wild-type pool
#       - one mutant pool

# the overall logic is that each haploid pool is treated as if it were diploid
# since the reads could potentially show two values at every genome position
# the two "pseduo-diploids" are compared to find genotype differences
# where the mutation you seek will be:
#   - "homozgyous", i.e., pure mutant in the mutant pool
#   - "homozygous", i.e., pure wild-type in the wild-type pool

#------------------------------------------------------------
# STEP 1 - align sample reads to the yeast with bwa mem
# prerequisites: install bwa and create the required genome index
# do this for each pool separately
#------------------------------------------------------------
SAMPLE=xyz
GENOME_FASTA=/path/to/sacCer3.fa
FASTQ_GLOB1=/path/to/sample/input/*.R1.fastq.gz
FASTQ_GLOB2=/path/to/sample/input/*.R2.fastq.gz
BAM_FILE=/path/to/sample/output/$SAMPLE.bam
#------------------------------------------------------------
bwa mem -R '@RG\tID:'$SAMPLE'\tSM:'$SAMPLE \
$GENOME_FASTA \
<(zcat $FASTQ_GLOB1) \
<(zcat $FASTQ_GLOB2) |
samtools view -hSb - |
samtools sort -o > $BAM_FILE

#------------------------------------------------------------
# STEP 2 - run samtools mpileup, merging the reads from two samples, e.g., wild-type and mutant pools
# prerequisites: install samtools and bcftools; ensure your bam files have RG (read group) tags
# note: you can also do this with a mutant pool and a parent heterozyous diploid, etc.
#------------------------------------------------------------
SAMPLE1=abc
SAMPLE2=xyz
BAM_FILE1=/path/to/sample/output/$SAMPLE1.bam
BAM_FILE2=/path/to/sample/output/$SAMPLE2.bam
PLOIDY=2 
BCF_FILE=/path/to/sample/output/$SAMPLE1.$SAMPLE2.bcf
#------------------------------------------------------------
samtools mpileup -DSu \
-f $GENOME_FASTA \
-d 250 \
-L 250 \
$BAM_FILE1 $BAM_FILE2 | 
bcftools view -bvcg \
-T pair \
-s <(echo -e "$SAMPLE1\t$PLOIDY\n$SAMPLE2\t$PLOIDY") - > $BCF_FILE
bcftools index $BCF_FILE

#------------------------------------------------------------
# STEP 3 - find genotype differences between two samples, e.g., the wild-type and mutant pools
# prerequisites: make a copy of this Perl script
#   - https://github.com/wilsontelab/pooled-segregant-public/find.pl
#------------------------------------------------------------
FIND_SCRIPT=/path/to/find.pl
DIFFERENCES_TSV=/path/to/sample/output/$SAMPLE1.$SAMPLE2.differences.tsv
#------------------------------------------------------------
bcftools view $BCF_FILE |
perl $FIND_SCRIPT > $DIFFERENCES_TSV
