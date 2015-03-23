### pe-sam-edit-v2.pl version 1.0 ### 
perl ~/tools/perl/pe-sam-edit-v2.pl -i [***.sam] --window [int] --ref [string] --trim5 [int] --strand --end [0, 1 or 2] --min [int] --max [int]
--input     the inputting sam file (required)
--window  the window size for bin-counting (required)
--ref          the reference fasta file (required)
--trim5     the number of trimed base of  5' side during arigment by bowtie
--strand   make two (f & r) files to sort reads depending on the orientation
--end       1:outputting the position of the edge R1, 2:R2, 0:center (default)
--min, --max   For size selection, minimum and maximum length of inserts (seq of interest)