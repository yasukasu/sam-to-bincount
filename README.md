# pe-sam-bincount.pl version 1.0
## Description
From the sam files from pair-end reads  
\- aligined reads are chosen (the outputted file: xxx.aln)  
\- sorted based based on the chromosome and the position where reads are aligned to (the outputted file: xxx.sort.aln)   
\- the number of reads in each bin (the indicated size) across the genome is counted (the outputted file: xxx.w[bin size].count.csv)  
These processes are sequencially excuted.
## Comand
perl pe-sam-to-bincount.pl -i [xxx.sam] --window [int] --ref [xxx.fasta]  --trim5 [int] --strand --end [0, 1 or 2] --min [int] --max [int]
## Options
--input:        the inputting sam file (required)  
--window:  the window size for bin-counting (required)  
--ref:          the reference fasta file (required)  
--trim5:     the number of trimed base of  5' side during arigment by bowtie  
--strand:   make two (f & r) files to sort reads depending on the orientation  
--end:       1:outputting the position of the edge R1, 2:R2, 0:center (default)  
--nei:        outputting postion shfited outward from the end by the indicated number (only in --end 1 or 2 mode)  
--min, --max:   For size selection, minimum and maximum length of inserts (seq of interest)
## Example
perl pe-sam-to-bincount.pl -i chip-seq.sam --window 300 --ref pombe.fasta  --trim5 2 --strand --end  1  --max 500