# /usr/bin/perl
use strict;
use Getopt::Long;

### initail parameters

my $trim5 = 0; 
my $maxc = 0;
my $minc = 0;
my $m1 = 0;
my $m2 = 0;
my $m3 = 0;
my $m4 = 0;
my $ss =0;
my $help=0;
my $end=0;
my $nei=0;
my $val=0;

my $file_aln = "output.aln";
my $filenames ="";
 
### processing  value entered in the command-line
my %opts=();
GetOptions(\%opts, "input=s", "window=i", "ref=s", "output=s" =>\$file_aln,  "trim5=i" => \$trim5,  "strand" => \$ss,  "end=i" => \$end, "min=i" => \$minc, "max=i" => \$maxc, "nei=i" => \$nei, "help" => \$help) or exit 1;

if($help){&help(); exit(0)}

foreach my $field ("input", "window", "ref"){
    if ( ! exists $opts{$field} ){
	&help();
	die "\nThe option '--$field' is required.\n"; 
    }
}

if($end>2){die "The option '--$end' must be '0', '1' or '2' (center/R1/R2).\n"; }

my $file_in = $opts{input};
my $win_s = $opts{window};
my $file_ref = $opts{ref};

unless(-e $file_in){print "$file_in does not exist!\n"; exit(1);}
unless(-e $file_ref){print "$file_ref does not exist!\n"; exit(1);}


print "### sam to bin-counting ###\n";
if($trim5){ print "5'-trimming before aligment:\t", $trim5, " bp\n";}
if($minc) { print "Minmum length of inserted seq:\t",$minc, " bp\n";}
if($maxc){ print "Maxmum length of inserted seq\t",$maxc, " bp\n";}
print "bin_size for bin-counting:\t",  $win_s, " bp\n";
print "referece fasta file:\t", $file_ref, "\n";
print "\n";

### Conversion of sam(PE) to aln file

# Setting files
#input

my $file_aln1 = "";
my $file_aln2 = "";


if($file_in =~ /.sam$/) {
    open IN, $file_in or die "$file_in could  - not be opened!:$!", "\n";
}
else {die "$file_in is not sam file.", "\n";}

my $mod ="";

if($minc && $maxc){ $mod = ".". $minc. "-". $maxc;}


if($end){$mod = $mod. ".e" . $end;}
else{$mod  = $mod.".cen";}

if(!$ss){
    my $ext = $mod . ".aln";
    
    $file_aln = $file_in;
    $file_aln =~ s/.sam$/$ext/;
    
    open OUT, ">$file_aln" or die "$file_aln could not be opened!:$!", "\n";

    $filenames = $file_aln;
}
else{
    my $ext1 = $mod . ".f.aln";
    my $ext2 = $mod . ".r.aln";     

    $file_aln1 = $file_in;
    $file_aln2 = $file_in;
    $file_aln1 =~ s/.sam$/$ext1/;
    $file_aln2 =~ s/.sam$/$ext2/;
  
    open OUT1, ">$file_aln1" or die "$file_aln could not be opened!:$!", "\n";
    open OUT2, ">$file_aln2" or die "$file_aln could not be opened!:$!", "\n";

    $filenames = $file_aln1. " & ". $file_aln2;
}

my $line = 0;

my $var1; #reads mapped in proper pair
my $var2; #fisrt read unmapped
my $var3; #secondary read unmapped

my $slt; #selected reads based on min/max mode

my $rn=0; # read length
my $d = 0; #length of cluster
my $pos_middle = 0;
 
print "## Extracting aligned PE-reads...\n";

while(<IN>){
	chop;
	if(substr($_,0,1) eq "@"){next;}

	$line++;

	(my $qname, my $flag, my $chro, my $pos, my $mapq, my $ciger, my $rname, my $rpos, my $len, my $seq, my @others) = split /\t+/, $_;

	if($flag & 64){ #selecting R1 elemtent

	    if($flag & 2){$var1++;} # the number of reads properly aligned according to the aligner
	    if($flag & 4){$var2++;} # R1 is unmapped
	    if($flag & 8){$var3++;} # R2 is unmapped
	    
	    $rn = length($seq);
	    $d = abs($len)+2*$trim5; # length before triming

	    if($minc || $maxc){
		if( $d<$minc && $minc){next;}
		if( $d>$maxc && $maxc){next;}
		$slt++;
	    }

	    if(!$ss){
		if(($flag & 2) && ($flag & 32)){   # selecting read in which both ends are properly aligned.  

		    if($end==1){$val = $pos-$trim5-$nei;} #left-side
		    elsif($end==2){$val = $rpos+$rn+$trim5-1+$nei;}#right-side
		    else{$val =int(($pos+$rpos)/2);}
		    print OUT $chro, "\t", $val, "\n";
		}
		elsif(($flag & 2) && ($flag & 16)){ 

		    if($end==1){$val = $pos+$rn+$trim5-1+$nei;}#right-side
		    elsif($end==2){$val = $rpos-$trim5-$nei;}
		    else{$val =int(($pos+$rpos)/2);}
		    print OUT $chro, "\t", $val, "\n";
		}
	    }
	    else{     
		if(($flag & 2) && ($flag & 32)){   # selecting read in which both ends are properly aligned.  

		    if($end==1){$val = $pos-$trim5-$nei;}
		    elsif($end==2){$val = $rpos+$rn+$trim5-1+$nei;}
		    else{$val = int(($pos+$rpos)/2);}
		    print OUT1 $chro, "\t", $val, "\n";
		}
		elsif(($flag & 2) && ($flag & 16)){ 

		    if($end==1){$val = $pos+$rn+$trim5-1+$nei;}
		    elsif($end==2){$val = $rpos-$trim5-$nei;}
		    else{$val = int(($pos+$rpos)/2);}
		    print OUT2 $chro, "\t", $val, "\n";
		}
	    }
		     
	    if($d < 200){$m1++;}
	    elsif($d <400){$m2++;}
	    elsif($d <600){$m3++;}
	    else {$m4++;}

	    # check d value	print $d, "\n";
	}
}

close(IN);
if(!$ss){close(OUT);}
else{close(OUT1); close(OUT2);}

print "produced file: ", $filenames,  "\nNumber of PE-clusters:\t", $line/2, "\n";
print "Number of aligned pairs:\t", $var1, "\n";
if($minc || $maxc){
    print "Selected pairs (", $minc, "-";
    if(!$maxc){ print " bp)"; }
    else { print $maxc, " bp)";}
    print "\t". $slt, "\n";
}
print "Length of inserted seq: <200 N=", $m1, "; 201-400 N=", $m2, "; 401-600 N=", $m3, "; >600 N=", $m4, "\n";print "Number of first reads that was not aligned :\t", $var2, "\n";
print "Number of second reads that was not aligned :\t", $var3, "\n";

###






### Sorting data in aln files (making .sort.aln files)

print "\n## Sorting aligned PE-reads...\n";

if(!$ss){&aln_sort($file_aln, $file_ref );}
else{
    print("f: watson str\n");
    &aln_sort($file_aln1, $file_ref); 
    print("r: Crick str\n");
    &aln_sort($file_aln2, $file_ref);
}

### Bin-counting

print "\n## Counting the number of reads in each $win_s bin...\n";
 
if(!$ss){
    (my $file_sorted  = $file_aln) =~ s/.aln$/.sort.aln/;
    &bin_count($file_sorted, $win_s, $file_ref, $trim5);
}
else{
    (my $file_sorted1  = $file_aln1) =~ s/.aln$/.sort.aln/;
    (my $file_sorted2  = $file_aln2) =~ s/.aln$/.sort.aln/;
    print("f: watson str\n");
    &bin_count($file_sorted1, $win_s, $file_ref, $trim5);
    print("r: Crick str\n");
    &bin_count($file_sorted2, $win_s, $file_ref, $trim5);
}



######### subrutines ##############

sub bin_count{

# setting options

    my $file_in = @_[0];
    my $win_swn = @_[1];  # window (bin) size (-w)
    my $file_ref = @_[2]; 
    my $trim5 = @_[3];
    my $file_out ="";

# file mode setting

    if($file_in !~ /.sort.aln$/) {die "$file_in is not a \"sorted\" file.", "\n";}
    else{ ($file_out = $file_in) =~ s/.sort.aln$/-w$win_s.count.csv/;}
  

    open IN, $file_in or die "$file_in could not be opened!:$!", "\n";
    open OUT, ">$file_out" or die "$file_out could not be opened!:$!", "\n";


# setting parmeters


    my $win_num = 0;   # the number of windows on the chromosome 
    my $chro_num = 0;

    my ($id_ref, $n_seq_ref) = extract_FASTA_info($file_ref);

    my @chro_name =@{$id_ref};
    my @chro_s = @{$n_seq_ref};

    my $count = 0; 
    my $chro_count = 0;
    my $chro_non_al = 0;
    my $win_x = 0; # position of the shiting window




# outputting row names 
    print OUT "chro,pos,count\n";

    while(<IN>){
	chop;
	
	(my $chro, my $pos) = split /\t+/, $_;


	if($chro eq $chro_name[$chro_num]){

	    # when pos is out of chromosome-range...
	    if($pos < 1 || $chro_s[$chro_num] < $pos){$chro_non_al++; next;}

  
          # when the read is within the bin... 
	    if($win_x < $pos && $pos <= ($win_x+$win_s)){$count++;}
	    
	    # when the read is out of the bin... 
	    elsif(($win_x+$win_s)< $pos){

		&output_reads_bin($chro_name[$chro_num], $win_x, $win_s, $count);

		# shifting bin...
		$win_x += $win_s;
		$win_num++; 
		$chro_count +=  $count;
		$count = 0;

		# outputting empty bins (if appliable)
		($win_x, $win_num) = output_empty_bin($win_x, $pos-$win_s, $chro_name[$chro_num], $win_s, $win_num);		



		#when the read mached in the bin...
		$count++;
	    }		
	    else {
		close(IN);
		close(OUT);
		die "Error: Unexpected value during shifing bins.\n";
	    }	        
	    
	}

        # when encountering the read on the next chromosome... 
	elsif($chro eq $chro_name[$chro_num+1]) {

	    # outputting the bin being dealt with
	    &output_reads_bin($chro_name[$chro_num], $win_x, $win_s, $count);

	    # shifting bin...
	    $win_x += $win_s;
	    $win_num++;
	    $chro_count +=  $count;
	    $count = 0;

	    # outputting the empty bins on the rest of chromosomal region 
	    ($win_x, $win_num) = output_empty_bin($win_x, $chro_s[$chro_num], $chro_name[$chro_num], $win_s, $win_num);


            #chromsome report		
	    &chro_report($chro_name[$chro_num], $chro_s[$chro_num],  $chro_count, $chro_non_al, $win_num, $win_x, $win_s); 
  	    
	    # moving to the next chromosome...
	    $chro_count = 0;
	    $chro_non_al = 0;
	    $chro_num++;
	    $win_x = 0;
 	    $win_num =0;
	    
	    # outputting empty bins
	    ($win_x, $win_num) = output_empty_bin($win_x, $pos-$win_s, $chro_name[$chro_num], $win_s, $win_num);
	    
            #when the read matched in the bin...
	    $count++;  
	}

	else {
	    close(IN);
	    close(OUT);
	    $" = ", ";
	    die "Error: Seq (chromosome): identity is not recognised. (name: $chro, ref seq contains '@chro_name')\n";
	}
    }
    



#outputting the last bin of the last hromosome
    &output_reads_bin($chro_name[$chro_num], $win_x, $win_s, $count);

# shifting bin...
    $win_x += $win_s;
    $win_num ++;
    $chro_count +=  $count;
    $count = 0;

# outputting the empty bins on the rest of chromosomal region
($win_x, $win_num) = output_empty_bin($win_x, $chro_s[$chro_num], $chro_name[$chro_num], $win_s, $win_num); 

#chromsome report
&chro_report($chro_name[$chro_num], $chro_s[$chro_num],  $chro_count, $chro_non_al, $win_num, $win_x, $win_s); 


close(IN);
close(OUT);

}



#outputting bins with reads
sub output_reads_bin{

    my $chro = @_[0];
    my $win_pos = @_[1];
    my $win_size = @_[2];
    my $count = @_[3];

    
    if($win_size>1){ print OUT $chro, ",", $win_pos+$win_size/2, ",",  $count, "\n";}
    else{ print OUT $chro, ",", $win_pos+1, ",",  $count, "\n";}
}	    
	  

# outputting the empty bins  
sub output_empty_bin{

    my $win_pos = @_[0];
    my $X2 = @_[1];
    my $chro_name = @_[2];
    my $win_size = @_[3];
    my $win_num =@_[4];



    while ($win_pos < $X2){
	&output_reads_bin($chro_name, $win_pos, $win_size, 0);
	$win_pos += $win_size;
	$win_num++;
    }
    
    return ($win_pos, $win_num);
}

# the report of the chromosome before
sub chro_report {

    my $chro_name = @_[0];
    my $chro_size = @_[1];
    my $chro_count =@_[2];
    my $chro_non_al =@_[3];
    my $win_num = @_[4];
    my $win_pos = @_[5];
    my $win_size = @_[6];


    print $chro_name, "\tsize:", $chro_size, " bp\ttotal reads:", $chro_count, "\tout of chro-range:", $chro_non_al, "\tbins:", $win_num, "\treads per window:", int($chro_count/$win_num), "\n";   
    if (($win_pos -$win_size) > $chro_size) {print "## warning! Some reads is out of the ", $chro_name, "-size. \n";}
} 

#obtaining IDs & lengths of sequences in fasta file 

sub extract_FASTA_info{

    my $file_seq =@_[0];
    my @id =();
    my @n_seq=();

    my $i=-1;

    open SEQ, $file_seq or die "Ref:$file_seq could not be opened!\n";

    while(<SEQ>){
	chomp;

	my $c = substr($_, 0, 1);
	
	if ($c eq ">"){
 
	    $i++;	    

	    $id[$i] =  (split /\s+/, substr($_, 1, length($_)))[0]; 
	    
	    next; 
	}
	
	if($id[$i] eq ""){ die("Error during processing the fasta file.\n");}
        $n_seq[$i] = $n_seq[$i]+length($_);


    }
    close(SEQ);
    return(\@id, \@n_seq);
}

# Opening aln files and sorting data inside

sub aln_sort{

    my $file_in = @_[0];
    my $file_ref = @_[1];    
    

    (my $file_out = $file_in) =~ s/.aln$/.sort.aln/;
    my ($chro_list_ref, $n_seq_ref) = extract_FASTA_info($file_ref);  

    open IN, $file_in  or die "Can't open $file_in.";
    my %group=();
    my $n=0;
    my %chk=();
    my $rstr = randstr(6);


    while (my $line = <IN>) {
	chomp($line);
	(my $chr, my $pos) = split(/\t/, $line); # chr name will be the key
	
	push @{$group{$chr}}, int($pos);

	if (@{$group{$chr}} == 1000) {   



	    my $tmpfile = $rstr.".".$chr."."."tmp";

	    if (-e $tmpfile &&  $chk{$chr}!=1){unlink $tmpfile;  print("The $tmpfile file is cleared.\n")}
	    $chk{$chr}=1;
	    
	    open OUT, ">>$tmpfile" or die "Can't open $tmpfile: $!";
	    print OUT join("\n", @{$group{$chr}})."\n";     
	    close OUT;
	    @{$group{$chr}} = ();          
	}
	$n++;
    }
    
    print("Total line: $n\n");
    

    close IN;
    open OUT, ">$file_out" or die "Can't open outfile: $!";

    my $tmp="";
    foreach my $chr (@{$chro_list_ref}) {  
	if (-e "$rstr.$chr.tmp") {              # if there are temp file
	    open IN, "$rstr.$chr.tmp" or die "Can't open $rstr.$chr.tmp: $!";
	    
	    while(my $tmp = <IN>){
		chomp($tmp);
		push @{$group{$chr}}, $tmp+0;     # inputting arrray data
	    }
	    close IN;
	    unlink "$rstr.$chr.tmp" or die "Can't delete $rstr.$chr.tmp: $!";    # deleting temp files
	}

	my $n_chr =@{$group{$chr}};
	print "$chr: $n_chr\n";

	foreach (sort {$a <=> $b} @{$group{$chr}}) {
	    print OUT "$chr\t$_\n"; #out
	}

	@{$group{$chr}} = ();
   }

close OUT;
}

sub randstr {
    my $length = $_[0];

    my @char_tmp=();

    # storing  characters
    # (a-z,A-Z, 0-9)
    push @char_tmp, ('a'..'z');
    push @char_tmp, ('A'..'Z');
    push @char_tmp, (0..9);

    # generating a indcated number of radom string.
    my $rand_str_tmp = '';
    my $i;
    for ($i=1; $i<=$length; $i++) {
    	$rand_str_tmp .= $char_tmp[int(rand($#char_tmp+1))];
    }

    return $rand_str_tmp;
}


sub help{
	print "###  Sam-to-Bincount version 1.04 ###\n"; 
	print "perl sam-to-bincount.pl -i [***.sam] --window [int] --ref [string] --trim5 [int] --strand --end [0, 1 or 2] --nei [int] --min [int] --max [int]\n\n";
	print "--input;\tthe inputting sam file (required)\n";
	print "--window;\tthe window size for bin-counting (required)\n";
	print "--ref;\t\tthe location of the reference fasta file (required)\n";
	print "--trim5;\tthe number of trimed base of  5' side during arigment by bowtie\n";
	print "--strand;\tmake two (f & r) files to sort reads depending on the orientation\n";
	print "--end;\t\t1:outputting the position of edge R1, 2:R2; 0:center (default)\n";
	print "--nei;\t\toutputting postion shfited outward from the end by the indicated number (only in --end 1 or 2 mode)\n";
	print "--min, --max;\tFor size selection, minimum and maximum length of inserts (seq of interest)\n";

}
