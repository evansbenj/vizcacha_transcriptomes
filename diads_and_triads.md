# This is the info on diad and triad analysis

On info, I am working here: 
```
/home/ben/2014_Tympanoctomys_transcriptomes/diads_and_triads
```

Here is the script I used to compile the blast results (Makes_diads_and_triads_min200.pl):

```perl
#!/usr/bin/perl
use warnings;
use strict;

# This program reads in two blast outputs.  One is trop unigene blasted 
# against laevis unigene blast db and saving top two hits.  The other is 
# laevis blasted against trop saving only top hit.  We need to find trop ids
# that have reciprocal best blast hit to one or two laevis unigene ids.


# on my computer I made a blast db like this:
# /usr/local/ncbi/blast/bin/makeblastdb -in /projects/Tympanoctomys/guineapig_1700_unique -dbtype nucl -out /projects/Tympanoctomys/guineapig_1700__unique_blastable
# /usr/local/ncbi/blast/bin/makeblastdb -in tuco_tuco_transcriptome -dbtype nucl -out tuco_tuco_transcriptome_blastable
# /usr/local/ncbi/blast/bin/makeblastdb -in /projects/Tympanoctomys/Tympanoctymys_contigs_and_singlets_tgicl_default -dbtype nucl -out /projects/Tympanoctomys/Tympanoctymys_contigs_and_singlets_tgicl_default_blastable

# To get the blast output that is parsed by this program, I used these commands
# /usr/local/ncbi/blast/bin/blastn -query Tympanoctymys_contigs_and_singlets_tgicl_default -db guineapig_1700__unique_blastable -outfmt 6 -out Tympanoctymys_contigs_and_singlets_tgicl_default_to_guineapig -evalue 1e-20 -task megablast -max_target_seqs 1
# /usr/local/ncbi/blast/bin/blastn -query guineapig_1700_unique -db Tympanoctymys_contigs_and_singlets_tgicl_default_blastable -outfmt 6 -out guineapig_to_Tympanoctymys_contigs_and_singlets_tgicl_default -evalue 1e-20 -task megablast -max_target_seqs 2

# /usr/local/ncbi/blast/bin/blastn -query Tympanoctymys_contigs_and_singlets_tgicl_default -db /projects/Tympanoctomys/tuco_tuco/tuco_tuco_transcriptome_blastable -outfmt 6 -out Tympanoctymys_contigs_and_singlets_tgicl_default_to_tuco -evalue 1e-20 -task megablast -max_target_seqs 1
# /usr/local/ncbi/blast/bin/blastn -query /projects/Tympanoctomys/tuco_tuco/tuco_tuco_transcriptome -db /projects/Tympanoctomys/Tympanoctymys_contigs_and_singlets_tgicl_default_blastable -outfmt 6 -out tuco_tuco_to_Tympanoctymys_contigs_and_singlets_tgicl_default -evalue 1e-20 -task megablast -max_target_seqs 2


# To check guineapig do this
# /usr/local/ncbi/blast/bin/blastn -query guineapig_Nov_18_2014_est.fa_unique -db guineapig_Nov_18_2014_est__unique_blastable -outfmt 6 -out Tympanoctymys_contigs_and_singlets_tgicl_default_to_guineapig -evalue 1e-20 -task megablast -max_target_seqs 3


# Add the additional criterion that the length of the blast hit must be at least $min_length in order for it to qualify.

#my $outputfile = "octomys_tympa_dyads_and_triads_reverse_minlength200";
my $outputfile = "octomys_tympa_dyads_and_triads_reverse_minlength200";

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile  $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";

#my $outputfile2 = "octomys_IDs_with_no_reciprocal_best_blast_hit_unigene_minlength200";
my $outputfile2 = "tympa_IDs_with_no_reciprocal_best_blast_hit_unigene_reverse";

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2  $!\n\n";
	exit;
}
print "Creating output file: $outputfile2\n";



# first open laevis blast results
#open (DATA, "Tymp_cdhit100_to_Octomys_cdhit100.out") or die "Failed to open Blast results";
open (DATA, "Oct_cdhit100_to_Tymp_cdhit100_save1.out") or die "Failed to open Blast results";


my @temp;
my %laevis_blast_results;
my %diad_triad_hash;
my %trop_IDs_with_no_reciprocal_best_blast_hit;
my $key;
my $min_length=200;

while ( my $line = <DATA>) {
	@temp = split("\t",$line);
	my $a=$temp[1];
	# assign the laevis hash irrespective of hit size
	$laevis_blast_results{$temp[0]} = $temp[1];
}	

close DATA;

#open (DATA2, "Oct_cdhit100_to_Tymp_cdhit100.out") or die "Failed to open Blast results";
open (DATA2, "Tymp_cdhit100_to_Octomys_cdhit100_save2.out") or die "Failed to open Blast results";

	while ( my $line2 = <DATA2>) {
		@temp = split("\t",$line2);
		if ((defined($laevis_blast_results{$temp[1]}))&&($laevis_blast_results{$temp[1]} eq $temp[0])){
			# this laevis seq matched the trop seq
			# now check if it is in a triad that has already been established
			if(defined ($diad_triad_hash{$temp[0]}[0]) ){
				# this is a check for duplicated entries:
				if($diad_triad_hash{$temp[0]}[0] ne $temp[1]){
					if($temp[3] >= $min_length){
						$diad_triad_hash{$temp[0]}[1]=$temp[1];
					}	
				}
			}
			else{
				if($temp[3] >= $min_length){
					$diad_triad_hash{$temp[0]}[0]=$temp[1];
				}
			}
		}
		
		elsif ((defined($laevis_blast_results{$temp[1]}))&&($laevis_blast_results{$temp[1]} ne $temp[0])){
			$trop_IDs_with_no_reciprocal_best_blast_hit{$temp[0]} = 1;
		}
		else{
			print "There was not a passable blast result for Tymp seq ",$temp[1],"\n";
		}	
	}

close DATA2;

foreach $key (sort(keys %diad_triad_hash)) {
	if(defined($diad_triad_hash{$key}[0])){
		print OUTFILE $key,"\t",$diad_triad_hash{$key}[0];
		if(defined $diad_triad_hash{$key}[1]){
			print OUTFILE "\t",$diad_triad_hash{$key}[1];
		}
		print OUTFILE "\n";
	}	
}

foreach $key (sort(keys %trop_IDs_with_no_reciprocal_best_blast_hit)) {
	print OUTFILE2 $key,"\n";
}
```
and to align the blast bits (Aligns_XLunigene_STensembl_triads_and_diads_onlyBLASTbits_baseml.pl):
```
#!/usr/bin/perl
use warnings;
use strict;

# This program will open a list of 
# triads and diads based on XL unigene and ST ensembl
# It will then open some files that have the sequences for
# each species.  

# It will then make fasta files that include the entire trop seq
# and the bits of the XL seq(s) that BLAST identified as homologous

# I'm doing this because using the entire XL sequences generated
# tons of bogus alignments

# Then the program will output triad and diad files in fasta
# format.  These will be aligned with mafft.  This will then be 
# converted with a program called Elconcatenero.py.  



my @temp;
my %bits;
my $key;
my $switch=0;
my $temp;
my $position;
my $string;
my $w;
my $y;
my $temparoo;
my @liner;

my @files;
my $status;
my $ST_name;
my $XL_name_0;
my $XL_name_1;
my @prefix;
my $temperoo;
my $previousline;
my %laevis_bits;

# read in the diad and triad info into a hash; 
# we can then use this information to retrieve
# the sequences

	
my %trop_key_to_laevis_seqs;

# first open laevis blast results
open (DATA, "octomys_tympa_dyads_and_triads_minlength200") or die "Failed to open diad triad file results";

while ( my $line = <DATA>) {
	@temp = split(/\s+/,$line);
	$trop_key_to_laevis_seqs{$temp[0]}[0]= $temp[1];
	if(defined($temp[2])){
		$trop_key_to_laevis_seqs{$temp[0]}[1]= $temp[2];
	}
}	
close DATA;

foreach $key (sort(keys %trop_key_to_laevis_seqs)) {
	print "fkey",$key,"\n";
	print "fofo",$trop_key_to_laevis_seqs{$key}[0],"\n";
}


# now I need to open up the trop to XL blast results and read in the coordinates of the XL match
open (DATA44, "Oct_cdhit100_to_Tymp_cdhit100.out") or die "Failed to open diad triad file results";


while ( my $line = <DATA44>) {
	@temp = split(/\s+/,$line);
	if(defined($laevis_bits{$temp[1]})){
		# we've seen this one before so expand the margins
		if($temp[8] < $temp[9]){
			if($temp[8] < $laevis_bits{$temp[1]}[0]){
				if(defined($trop_key_to_laevis_seqs{$temp[0]}[0])){
					if($trop_key_to_laevis_seqs{$temp[0]}[0] eq $temp[1]){
						$laevis_bits{$temp[1]}[0]=$temp[8];
					}
				}		
				if(defined($trop_key_to_laevis_seqs{$temp[0]}[1])){
					if($trop_key_to_laevis_seqs{$temp[0]}[1] eq $temp[1]){
						$laevis_bits{$temp[1]}[0]=$temp[8];
					}
				}
			}
			else{
				if(defined($trop_key_to_laevis_seqs{$temp[0]}[0])){
					if($trop_key_to_laevis_seqs{$temp[0]}[0] eq $temp[1]){
						$laevis_bits{$temp[1]}[1]=$temp[9];
					}
				}		
				if(defined($trop_key_to_laevis_seqs{$temp[0]}[1])){
					if($trop_key_to_laevis_seqs{$temp[0]}[1] eq $temp[1]){
						$laevis_bits{$temp[1]}[1]=$temp[9];
					}
				}	
			}
		}
		else{
			if($temp[9] < $laevis_bits{$temp[1]}[0]){
				if(defined($trop_key_to_laevis_seqs{$temp[0]}[0])){
					if($trop_key_to_laevis_seqs{$temp[0]}[0] eq $temp[1]){
						$laevis_bits{$temp[1]}[0]=$temp[9];
					}
				}		
				if(defined($trop_key_to_laevis_seqs{$temp[0]}[1])){
					if($trop_key_to_laevis_seqs{$temp[0]}[1] eq $temp[1]){
						$laevis_bits{$temp[1]}[0]=$temp[9];
					}
				}	
			}
			else{
				if(defined($trop_key_to_laevis_seqs{$temp[0]}[0])){
					if($trop_key_to_laevis_seqs{$temp[0]}[0] eq $temp[1]){
						$laevis_bits{$temp[1]}[1]=$temp[8];
					}
				}		
				if(defined($trop_key_to_laevis_seqs{$temp[0]}[1])){
					if($trop_key_to_laevis_seqs{$temp[0]}[1] eq $temp[1]){
						$laevis_bits{$temp[1]}[1]=$temp[8];
					}
				}	
			}		
		}
	}
	else{
		if($temp[8] < $temp[9]){
			if(defined($trop_key_to_laevis_seqs{$temp[0]}[0])){
				if($trop_key_to_laevis_seqs{$temp[0]}[0] eq $temp[1]){
					$laevis_bits{$temp[1]}[0]=$temp[8];
					$laevis_bits{$temp[1]}[1]=$temp[9];
				}
			}	
			if(defined($trop_key_to_laevis_seqs{$temp[0]}[1])){
				if($trop_key_to_laevis_seqs{$temp[0]}[1] eq $temp[1]){
					$laevis_bits{$temp[1]}[0]=$temp[8];
					$laevis_bits{$temp[1]}[1]=$temp[9];
				}
			}
		}
		else{
			if(defined($trop_key_to_laevis_seqs{$temp[0]}[0])){
				if($trop_key_to_laevis_seqs{$temp[0]}[0] eq $temp[1]){
					$laevis_bits{$temp[1]}[1]=$temp[8];
					$laevis_bits{$temp[1]}[0]=$temp[9];
				}
			}	
			if(defined($trop_key_to_laevis_seqs{$temp[0]}[1])){
				if($trop_key_to_laevis_seqs{$temp[0]}[1] eq $temp[1]){
					$laevis_bits{$temp[1]}[1]=$temp[8];
					$laevis_bits{$temp[1]}[0]=$temp[9];
				}
			}	
		}
	}	
}	
close DATA;

my $outputfile;



foreach $key (sort(keys %trop_key_to_laevis_seqs)) {
	#print "mkey",$key,"\n";
	if(exists($trop_key_to_laevis_seqs{$key}[0])){
		print $trop_key_to_laevis_seqs{$key}[0],"\n";	
	}
	else{
		delete $trop_key_to_laevis_seqs{$key},"\n";
	}
}


foreach $key (sort(keys %trop_key_to_laevis_seqs)) {
	print "key",$key,"h",$trop_key_to_laevis_seqs{ $key }[0],"\n";
	my $outputfile = "./alignments_minlength200_blastbits/triads_and_diads_blastbit_".$key.".fa";
	unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile  $!\n\n";
	exit;
	}
	
	print "Creating output file: $outputfile\n";

	# open ST data
	open (DATA2, "Octomys_all_transcriptomes_assembled_together_unique.fasta") or die "Failed to open trop data results";
	$switch=0;
	while ( my $line = <DATA2>) {
		chomp($line);
		@liner=split(/\s+/,$line);
		if($switch == 1){
			if($liner[0] =~ />/){
				print OUTFILE $previousline,"\n";
				$switch = 2;
				last;
			}
			else{
				print OUTFILE $previousline,"\n";
				$previousline = $line;
			}
		}
		elsif($liner[0] =~/$key/){
			$previousline = $liner[0];
			$switch=1;
		}
	}
	close DATA2;


	my $temp_laevis1="";


	# now get the XL seq(s)
	open (DATA3, "Tympa_all_transcriptomes_assembled_together_unique.fasta") or die "Failed to open laevis data first time";
	$switch=0;
	$temp_laevis1="";
	print "looking for ",$trop_key_to_laevis_seqs{ $key }[0],"\n";
	while ( my $line = <DATA3>) {
		chomp($line);
		@liner=split(/\s+/,$line);
		if($switch == 1){
			if($liner[0] =~ />/){
				# ok we are at the next entry so stop
				$switch = 2;
				last;
			}
			else{
				$temp_laevis1=$temp_laevis1.$line;
			}
		}
		elsif($liner[0] eq ">".$trop_key_to_laevis_seqs{ $key }[0] ){
			$XL_name_0 = substr($liner[0],1);
			print "XL name ",$XL_name_0,"\n";
			$switch=1; 
			print "found it\n";
			print OUTFILE ">",$XL_name_0,"\n";
		}
	}
	close DATA3;
	# now print only the bit of this seq that we want
	print OUTFILE substr($temp_laevis1,$laevis_bits{$XL_name_0}[0],($laevis_bits{$XL_name_0}[1]-$laevis_bits{$XL_name_0}[0])),"\n";

	if(exists($trop_key_to_laevis_seqs{ $key } [1])){
	print "this is a triad\n";
		print "looking for ",$trop_key_to_laevis_seqs{ $key }[1],"\n";
		open (DATA4, "Tympa_all_transcriptomes_assembled_together_unique.fasta") or die "Failed to open laevis data first time";
		$switch=0;
		$temp_laevis1="";
		while ( my $line = <DATA4>) {
			chomp($line);
			@liner=split(/\s+/,$line);
			if($switch == 1){
				if($liner[0] =~ />/){
					# ok we are at the next entry so stop
					$switch = 2;
					last;
				}
				else{
					$temp_laevis1=$temp_laevis1.$line;
				}
			}
			elsif($liner[0] eq ">".$trop_key_to_laevis_seqs{ $key } [1]){
				$previousline= $liner[0];
				$XL_name_1 = substr($liner[0],1);
				print "XL name ",$XL_name_1,"\n";
				$switch=1;
				print "found it\n";
				print OUTFILE ">",$XL_name_1,"\n";
			}
		}
		close DATA4;
		print OUTFILE substr($temp_laevis1,$laevis_bits{$XL_name_1}[0],($laevis_bits{$XL_name_1}[1]-$laevis_bits{$XL_name_1}[0])),"\n";

	}
	close OUTFILE;
	print "Doing mafft alignment ";
	print "with outputfile ",$outputfile," yeah\n";
	# replace bars with backslash bar so that the maaft command works
	$outputfile =~ s/\|/\\\|/g;
	# replace forwardslash with backslash forwardslaps so that the maaft command work
	$outputfile =~ s/\//\\\//g;
	print "with outputfile ",$outputfile," yeah\n";
	$status = system("mafft --adjustdirectionaccurately --genafpair $outputfile > $outputfile.aligned");
	# now get rid of the backslash bars so that the base ml command works
	$outputfile =~ s/\\//g;

# I am commenting out the baseml stuff because it doesn't work with the ens72 and/or XL_Taejoon names

#	if(exists($trop_key_to_laevis_seqs{ $key } [1])){
#	# this is a triad so make paml files (data and tree files, and then run baseml
#		# make datafile
#		print "making paml outputfile: ",$outputfile,"\n";
#		@prefix= split(/\./,$outputfile);
#		$temperoo=$prefix[0];
#		print $prefix[0], " hi \n";
#		# replace bars with backslash bar so that the ElConcatenero command works
#		$outputfile =~ s/\|/\\\|/g;
#
#		$status = system("python ElConcatenero.py -c -if fasta -of phylip -in $outputfile.aligned -o $temperoo.phy");
#		# now get rid of the backslash bars 
#		$outputfile =~ s/\\//g;
#		
#		# search for exclamation points inserted by macseand replace with dash for gap
#		# also get rid of number signs and bars in the species names
#		$status = system("sed -i \"\" 's/!/-/g' $temperoo.phy");
#		$status = system("sed -i \"\" 's/|/_/g' $temperoo.phy");
#		$status = system("sed -i \"\" 's/#/_/g' $temperoo.phy");
#		
#		# make the paml control file
#		my $outputfile4 = $temperoo."_baseml.ctl";
#		unless (open(OUTFILE4, ">$outputfile4"))  {
#			print "I can\'t write to $outputfile4   $!\n\n";
#			exit;
#		}
#		print OUTFILE4 "      seqfile = ./$temperoo.phy\n";
#		print OUTFILE4 "     treefile = ./$temperoo.tre\n";
#		print OUTFILE4 "      outfile = ./$temperoo.mlb\n";	
#
#		print OUTFILE4 "        noisy = 0  * 0,1,2,3,9: how much rubbish on the screen\n";
#		print OUTFILE4 "      verbose = 1  * 0: concise; 1: detailed, 2: too much\n";
#		print OUTFILE4 "      runmode = 0  * 0: user tree;  1: semi-automatic;  2: automatic\n";
#		print OUTFILE4 "                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise\n";
#
#
#		print OUTFILE4 "        model = 7   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85\n";
#		print OUTFILE4 "                            * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu\n";
#
#		print OUTFILE4 "        clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis\n";
#		print OUTFILE4 "    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated\n";
#		print OUTFILE4 "        kappa = 5  * initial or fixed kappa\n";
#		print OUTFILE4 "    fix_alpha = 0  * 0: estimate alpha; 1: fix alpha at value below \n";
#		print OUTFILE4 "        alpha = 0.5  * initial or fixed alpha, 0:infinity (constant rate)\n";
#
#		print OUTFILE4 "       Malpha = 0  * different alphas for genes\n";
#		print OUTFILE4 "        ncatG = 5  * # of categories in the dG, AdG, or nparK models of rates\n";
#
#		print OUTFILE4 "        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK \n";
#		print OUTFILE4 "        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2\n";
#		print OUTFILE4 "        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates\n";
#		print OUTFILE4 "        RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states\n";
#
#		print OUTFILE4 "   Small_Diff = 1e-7\n";
#		print OUTFILE4 "    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?\n";
#		print OUTFILE4 "       method = 0   * 0: simultaneous; 1: one branch at a time\n";
#		
#		close OUTFILE4;
#		
#		# make tree file
#		my $outputfile5 = $temperoo.".tre";
#		unless (open(OUTFILE5, ">$outputfile5"))  {
#			print "I can\'t write to $outputfile5   $!\n\n";
#			exit;
#		}
#		print OUTFILE5 "3  1\n";
#		print OUTFILE5 "\(",$ST_name,",",$XL_name_0," #1,",$XL_name_1," #2\);\n\n";
#	
#		$status = system("./baseml_4.8_longname $temperoo\_baseml.ctl");
#	}
}
```
and to get patristic distances (gets_patristic_distances_from_fasta_minlength200.pl):
```perl
#!/usr/bin/perl 

use strict;

######
#
#	This program reads in a bunch of aligned fasta files
#	and calculates the pairwise divergence between one or both
# 	laevis paralogs.  The distances from baseml seem unreliable
#	so hopefully this will be better.
#
#	We include only characters that are not a gap between the two
#	XL seqs and not an "N" in either sequence
#
#	We also only report alignments that have at least $min_length
#	positions that are aligned, ungapped, and not Ns
#
#####

my @files;
my @line;
my $length;
my $switch=0;

my $outputfile = "tymp_tet_oct_dip_triad_and_diad_divergences_minlength200";

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile   $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";

print OUTFILE "Oct_ID\tTym_or_D\tbp\tdivergence_from_Oct\tparalog_bp\ttriad_paralog_divergence\n";


@files = glob("./alignments_minlength200_blastbits/*blastbit*.fa.aligned");
$switch=0;
my $counter;
my $trop;
my @trop;
my $laevis1;
my $laevis2;
my @laevis1;
my @laevis2;
my $z;
my $laevis1_divergence=0;
my $laevis2_divergence=0;
my $laevis1_laevis2_divergence=0;
my $laevis1_bp=0;
my $laevis2_bp=0;
my $laevis1_laevis2_bp=0;
my $tropseq;
my $min_length=200;

foreach (@files) {
	open (DATAINPUT, "$_") or die "Failed to open trop data results";
	$counter=0;
	$trop=();
	$laevis1=();
	$laevis2=();
	$laevis1_divergence=0;
	$laevis2_divergence=0;
	$laevis1_laevis2_divergence=0;
	$laevis1_bp=0;
	$laevis2_bp=0;
	$laevis1_laevis2_bp=0;
	while ( my $line = <DATAINPUT>) {
		#@line = split (/\s+/,$line);
		chomp($line);
		if((substr($line,0,1) ne ">")&&($counter == 1)){
			$trop = $trop.$line;
		}
		elsif((substr($line,0,1) ne ">")&&($counter != 1)){
			if($counter == 2){
				$laevis1 = $laevis1.$line;
			}
			elsif($counter == 3){
				$laevis2 = $laevis2.$line;
			}
		}
		elsif(substr($line,0,1) eq ">"){
			$counter+=1;
			if($counter == 1){
				$tropseq = substr($line,1);
				print $tropseq,"\n";
			}
		}
	}#end while   
	#   close DATAINPUT;
	close DATAINPUT;
		
	# now compare them
	@trop=();
	@laevis1=();
	@laevis2=();
	@trop = split ("",$trop);
	@laevis1 = split ("",$laevis1);
	if(length($laevis2) > 0){
		@laevis2 = split ("",$laevis2);
	}	
	if(length($laevis2) > 0){	
		# this means it is a triad 
		# calculate st divergence for first paralog
		for ($z=0; $z<= $#laevis1; $z += 1){
			# first calculate divergence between laevis 1 and trop
			if(($trop[$z] ne "-")&&($laevis1[$z] ne "-")&&(uc($trop[$z]) ne "N")&&(uc($laevis1[$z]) ne "N")){
				if(uc($trop[$z]) ne uc($laevis1[$z])){
					$laevis1_divergence+=1;
					#print uc($laevis1[$z]),"h",uc($laevis2[$z]),"\n";
				}
				$laevis1_bp+=1;
			}
		}
		# calculate st divergence for second paralog
		for ($z=0; $z<= $#laevis1; $z += 1){
			# first calculate divergence between laevis 1 and trop
			if(($trop[$z] ne "-")&&($laevis2[$z] ne "-")&&(uc($trop[$z]) ne "N")&&(uc($laevis2[$z]) ne "N")){
				if(uc($trop[$z]) ne uc($laevis2[$z])){
					$laevis2_divergence+=1;
					#print uc($laevis1[$z]),"h",uc($laevis2[$z]),"\n";
				}
				$laevis2_bp+=1;
			}
		}
		# calculate divergence between paralogs
		for ($z=0; $z<= $#laevis1; $z += 1){
			# first calculate divergence between laevis 1 and trop
			if(($laevis1[$z] ne "-")&&($laevis2[$z] ne "-")&&(uc($laevis1[$z]) ne "N")&&(uc($laevis2[$z]) ne "N")){
				if(uc($laevis1[$z]) ne uc($laevis2[$z])){
					$laevis1_laevis2_divergence+=1;
					#print uc($laevis1[$z]),"h",uc($laevis2[$z]),"\n";
				}
				$laevis1_laevis2_bp+=1;
			}
		}
		#print average of both distances
		if(($laevis1_bp>0)&&($laevis2_bp>0)){
			print OUTFILE $tropseq,"\ttriad\t",($laevis1_bp+$laevis2_bp)/2,"\t",(($laevis1_divergence/$laevis1_bp)+($laevis2_divergence/$laevis2_bp))/2,"\t";
		}
		if(($laevis1_bp>0)&&($laevis2_bp>0)&&($laevis1_laevis2_bp>100)){
			print OUTFILE "\t",$laevis1_laevis2_bp,"\t",$laevis1_laevis2_divergence/$laevis1_laevis2_bp,"\n";     
		}
		elsif(($laevis1_bp>0)&&($laevis2_bp>0)&&($laevis1_laevis2_bp<100)){
			print OUTFILE "\tNAN\n";
		} 
	}
	else{
		# this means it is a diad 
		# calculate divergence for first paralog only
		for ($z=0; $z<= $#laevis1; $z += 1){
			# first calculate divergence between laevis 1 and trop
			if(($trop[$z] ne "-")&&($laevis1[$z] ne "-")&&(uc($trop[$z]) ne "N")&&(uc($laevis1[$z]) ne "N")){
				if(uc($trop[$z]) ne uc($laevis1[$z])){
					$laevis1_divergence+=1;
				}
				$laevis1_bp+=1;
			}
		}
		#print first divergence only
		print OUTFILE $tropseq,"\tdiad\t",$laevis1_bp,"\t",$laevis1_divergence/$laevis1_bp,"\n";     
	}


}

#  close the output file1
close OUTFILE;   
print "Closing output file: $outputfile\n";


exit;
```
