## Yeast pooled segregant analysis

This repository has code to perform yeast pooled segregant analysis after
[Birkeland et al.](https://academic.oup.com/genetics/article/186/4/1127/6063702?login=false)
but using modern code tools.

Unfortunately, the code is not supported and not completely formed as a pipeline to run as is, but the key bits are present.

### Usage

Examine the following script
and use its parts to finish assembly of your pipeline scripts.
- pooled-segregant.sh

As noted in the file, the script is not intended to run as is, 
you should mine it for parts, which you will have to troubleshoot.
To run all parts, you will need to copy and use this perl script:
- find.pl

The final result is a TSV (tab separated values) file that you 
can mine for interesting mutations using command line tools or Excel.
A basic strategy is to load your bam files into the 
[Integrated Genomics Viewer (IGV)](https://software.broadinstitute.org/software/igv/)
and go to regions with high CLR values indicative of genotype differences. 
