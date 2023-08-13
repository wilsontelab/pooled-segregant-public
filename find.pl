#!/usr/bin/perl
use strict;
use warnings;

# IMPORTANT NOTE about CLR!!
# CLR is "Phred log ratio of genotype likelihoods
#         with and without the trio/pair constraint"
# Note: _genotype_ likelihood, not of allele frequencies
# Thus, if genotypes were called the same, CLR is NOT REPORTED!
# In other words, CLR is null (not present) for same-genotype variants.
# If insufficient, consider using FQ (see below)
# However, for pooled backcrosses, seek genotype 0/0 vs. 1/1,
# so realistically the interesting mutants should always have CLR reported.

# collect the CLR scores for all variants
my ($nCLR, $nNoCLR, %CLRs, %vars) = (0, 0);
while(my $var = <STDIN>){
    chomp $var;
    $var =~ m|^\#| and next;
    if($var =~ m|CLR=(.+?)\;|){
        my $CLR = $1;
        my @f = split("\t", $var);
        my $s = $f[1]-101;
        my $e = $f[1]+100;
        $CLRs{$CLR}++;
        my $reg = "$f[0]:$s-$e";
        my $wtGen  = "_".(split(":", $f[$#f-1]))[0];
        my $mutGen = "_".(split(":", $f[$#f  ]))[0];
        push @{$vars{$CLR}}, join("\t", $reg, $f[0], $f[1], $f[3], $f[4], $wtGen, $mutGen);
        $nCLR++;
    } else {
        $nNoCLR++;
    }
}
print STDERR "$nCLR\tCLR found\n$nNoCLR\tCLR not found\n\n";

# print the histogram and sorted difference lists
print join("\t", qw(CLR REGION CHROM POS REF ALT WT_GEN MUT_GEN)), "\n";
foreach my $CLR(sort {$b<=>$a} keys %CLRs){
    print STDERR "$CLR\t$CLRs{$CLR}\n";
    foreach my $var(@{$vars{$CLR}}){
        print "$CLR\t$var\n";
    }
}


##fileformat=VCFv4.1
##samtoolsVersion=0.1.18 (r982:295)
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">

##INFO=<ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same">

##INFO=<ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)">
##INFO=<ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count(no HWE assumption)">
##INFO=<ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies">
##INFO=<ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3">

##INFO=<ID=CLR,Number=1,Type=Integer,Description="Log ratio of genotype likelihoods with and without the constraint">

##INFO=<ID=UGT,Number=1,Type=String,Description="The most probable unconstrained genotype configuration in the trio">
##INFO=<ID=CGT,Number=1,Type=String,Description="The most probable constrained genotype configurationin the trio">
##INFO=<ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=PC2,Number=2,Type=Integer,Description="Phred probability of the nonRef allele frequency in group1 samples being larger (,smaller) than in group2.">
##INFO=<ID=PCHI2,Number=1,Type=Float,Description="Posterior weighted chi^2 P-value for testing the association between group1 and group2 samples.">
##INFO=<ID=QCHI2,Number=1,Type=Integer,Description="Phred scaled PCHI2.">
##INFO=<ID=PR,Number=1,Type=Integer,Description="# permutations yielding a smaller PCHI2.">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  vac11_WT        vac11_mut

# with CLR
#chrI    13807   .       T       G       3.14    .       DP=287;VDB=0.0361;AF1=0.3333;CLR=3;AC1=1;DP4=118,99,30,40;MQ=51;FQ=4.13;PV4=0.1,8.7e-115,1,1  GT:PL:DP:SP:GQ  0/1:32,0,255:150:18:28  0/0:0,3,255:137:2:8
#
# without CLR
#chrI    11321   .       G       A       268     .       DP=186;VDB=0.0386;AF1=0.5;AC1=2;DP4=39,48,43,51;MQ=60;FQ=271;PV4=1,5.2e-100,1,1       GT:PL:DP:SP:GQ  0/1:143,0,255:101:1:99  0/1:160,0,255:80:1:99

