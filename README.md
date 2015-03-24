# pe-sam-edit-v2.pl version 1.0 
***
## Description
From the sam files from pair-end reads  
- aligined reads are chosen (the outputted file: xxx.aln)  
- sorted based based on the chromosome and the position where reads are aligned to (the outputted file: xxx.sort.aln)   
- the number of reads in each bin (the indicated size) across the genome is counted (the outputted file: xxx.w[bin size].count.csv) These proess sequencially are excuted.  
***
##Comannd
perl ~/tools/perl/pe-sam-edit-v2.pl -i [***.sam] --window [int] --ref [string] (--trim5 [int] --strand --end [0, 1 or 2] --min [int] --max [int])  
***
##Options
--input:        the inputting sam file (required)  
--window:  the window size for bin-counting (required)  
--ref:          the reference fasta file (required)  
--trim5:     the number of trimed base of  5' side during arigment by bowtie  
--strand:   make two (f & r) files to sort reads depending on the orientation  
--end:       1:outputting the position of the edge R1, 2:R2, 0:center (default)  
--min, --max:   For size selection, minimum and maximum length of inserts (seq of interest)  