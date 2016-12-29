# Genome assembly for Tymp and Oct using Abyss

(On iqaluk) I am working in this directory:
```
/work/ben/2017_Tymp_Oct_HiSeqX
```
First unload this module:
```
module unload intel mkl openmpi
```
Then load other modules:
```
module load gcc/4.9.2
module load openmpi/gcc-4.9.2/std/1.8.7
module load boost/gcc492-openmpi187/1.59.0
export PATH=/work/ben/abyss/bin:$PATH
```
This worked:
```
module load gcc/4.9.2
module load openmpi/gcc492-std/1.8.7
module load boost/gcc492-openmpi187std/1.59.0
export PATH=/work/ben/abyss/bin:$PATH
```

Assemble genome of each species using only paired end reads. 
```
abyss-pe np=8 name=AO245 k=64 in='AO245_R1_trim_paired.fastq.gz AO245_R2_trim_paired.fastq.gz' 
abyss-pe np=8 name=AO248 k=64 in='AO248_R1_trim_paired.fastq.gz AO248_R2_trim_paired.fastq.gz' 
```

Assembly of Octomys genome is completed.  Here are the stats:
```
The minimum coverage of single-end contigs is 0.472727.
The minimum coverage of merged contigs is 2.83636.
Consider increasing the coverage threshold parameter, c, to 2.83636.
ln -sf AO248-8.fa AO248-scaffolds.fa
PathOverlap --overlap   -k64 --dot AO248-7.dot AO248-7.path >AO248-8.dot
ln -sf AO248-8.dot AO248-scaffolds.dot
abyss-fac  AO248-unitigs.fa AO248-contigs.fa AO248-scaffolds.fa |tee AO248-stats.tab
n       n:500   L50     min     N80     N50     N20     E-size  max     sum    name
6294392 962523  180653  500     1677    3812    7731    5100    68869   2.404e9AO248-unitigs.fa
5585105 819908  139635  500     2035    4872    10309   6740    73917   2.44e9 AO248-contigs.fa
5515621 783585  129154  500     2143    5223    11254   7336    74984   2.448e9AO248-scaffolds.fa
ln -sf AO248-stats.tab AO248-stats
tr '\t' , <AO248-stats.tab >AO248-stats.csv
abyss-tabtomd AO248-stats.tab >AO248-stats.md
```

directory:
`/work/ben/2017_Tymp_Oct_HiSeqX/AO248_quaked_data`
assembly is either `AO248-contigs.fa` or `AO248-unitigs.fa`


# Assembling repeat elements using Repark from wgs reads

From within this directory: `/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS`
```
../RepARK.pl -l AO245_R1_trim_paired.cor.fastq.gz -l AO245_R2_trim_paired.cor.fastq.gz -l AO245_R1_trim_paired.cor_single.fastq.gz -l AO245_R1_trim_single.cor.fastq.gz -l AO245_R2_trim_paired.cor_single.fastq.gz -l AO245_R2_trim_single.cor.fastq.gz -k 31 -o repArc_kmer_31
```
and this directory `/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS`
```
../RepARK.pl -l AO248_R1_trim_paired.cor.fastq.gz -l AO248_R2_trim_paired.cor.fastq.gz -l AO248_R1_trim_paired.cor_single.fastq.gz -l AO248_R2_trim_paired.cor_single.fastq.gz -l AO248_R1_trim_single.cor.fastq.gz -l AO248_R2_trim_single.cor.fastq.gz -k 31 -o repArc_AO248_kmer_31
```
I also did this with a k-mer size of 45. I used this one liner to output the highest coverage contig in the velvet_repeat_lib directory:
```
grep -o -P '(?<=cov_).*(?=)' contigs.fa | sort -rn | head -n 1
```
This gets the headers of the entries with the highest coverage:
```
grep -o -P '(?<=cov_).*(?=)' contigs.fa | sort -rn | head -n 10
```

print out a fasta file
```
awk -v seq="NODE_9_length_29_cov_3.000000" -v RS='>' '$1 == seq {print RS $0}' contigs.fa
```

Make a file with all of the contig lengths from RepARK:
```
grep -o -P '(?<=length_).*(?=_cov)' contigs.fa > length_of_all_highabundancekmers
```
and one with the coverage of each contig:
```
grep -o -P '(?<=cov_).*(?=)' contigs.fa > coverage_of_all_highabundancekmers
```

# Make a blast database out of the genome assembly
```
/usr/local/blast/2.3.0/bin/makeblastdb -in /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248-scaffolds.fa -dbtype nucl -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248-scaffolds.fa_blastable
```
and to blast the repeats from Octomys against the genome assembly
```
/usr/local/blast/2.3.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa -db /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248-scaffolds.fa_blastable -outfmt 6 -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/repArc_AO248_kmer_31_to_AO248-scaffolds.fa_blastable -evalue 1e-20 -task megablast 
```
In the above command, could add `-max_target_seqs 1` if we want to restrict the output to some number of entries.

# Also make a blast database out of both transcriptome assemblies

The assemblies are here
```
/home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together.fasta
/home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together.fasta
```
Here's the commands to make the blast dbs out of the transcriptome assemblies:
```
/usr/local/blast/2.3.0/bin/makeblastdb -in /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together.fasta -dbtype nucl -out /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together.fasta_blastable
/usr/local/blast/2.3.0/bin/makeblastdb -in /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together.fasta -dbtype nucl -out /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together.fasta_blastable
```

and to blast the repeats from Octomys against the Octomys and tympa transcriptome assemblies
```
/usr/local/blast/2.3.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa -db /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together.fasta_blastable -outfmt 6 -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa_blasted_to_octomys_transcriptiome -evalue 1e-20 -task megablast 

/usr/local/blast/2.3.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/repArc_kmer_31/velvet_repeat_lib/contigs.fa -db /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together.fasta_blastable -outfmt 6 -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/repArc_kmer_31/velvet_repeat_lib/contigs.fa_blasted_to_tympa_transcriptiome -evalue 1e-20 -task megablast 
```

# Examining the incidence of repetitive elements in the RNAseq Assembly

I wrote a script to quantify repeats in UTRs and CDS for each species using gff files male from TransDecoder and the blast results generated above.  The paths for each of these files are hardcoded in the script; just change the comments for Oct and tymp.

The script is in this directory:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa_blasted_to_octomys_transcriptiome
```
and is called 'quantifies_repeats_in_RNAseq.pl':

```
#!/usr/bin/perl
use warnings;
use strict;

# This program reads in a gff file for a transciptome assembly and a 
# blast output from high abundence kmer contigs to a RNAseq transcriptome assembly

# it then quantifies the proportion of coding regions and non-coding regions that 
# are repeats

# here is the format of the blast output:

# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore, where
#  1.	 qseqid	 query (e.g., gene) sequence id
#  2.	 sseqid	 subject (e.g., reference genome) sequence id
#  3.	 pident	 percentage of identical matches
#  4.	 length	 alignment length
#  5.	 mismatch	 number of mismatches
#  6.	 gapopen	 number of gap openings
#  7.	 qstart	 start of alignment in query
#  8.	 qend	 end of alignment in query
#  9.	 sstart	 start of alignment in subject
#  10.	 send	 end of alignment in subject
#  11.	 evalue	 expect value
#  12.	 bitscore	 bit score

#my $outputfile = "Repeats_in_RNAseq_assembly_Octomys.out";
#unless (open(OUTFILE, ">$outputfile"))  {
#	print "I can\'t write to $outputfile  $!\n\n";
#	exit;
#}
#print "Creating output file: $outputfile\n";

# open the gff file results
# for octomys
open (DATAgff, "/home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta.transdecoder_dir/longest_orfs.gff3") or die "Failed to open laevis Blast results";
# for tympa
#open (DATAgff, "/home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta.transdecoder_dir/longest_orfs.gff3") or die "Failed to open laevis Blast results";

my @temp;
my %gff_genelength;
my %gff_CDSlength;
my %gff_UTRlength;
my %gff_CDS;

while ( my $line = <DATAgff>) {
	@temp = split("\t",$line);
	if($#temp>1){
		if($temp[2] eq 'gene'){
			unless(defined ($gff_genelength{$temp[0]})){ # just go with the first entry
				$gff_genelength{$temp[0]} = $temp[4];
			}	
		}
		elsif($temp[2] eq 'CDS'){
			unless(defined ($gff_CDS{$temp[0]})){ # just go with the first entry
				$gff_CDS{$temp[0]}[0]=$temp[3];
				$gff_CDS{$temp[0]}[1]=$temp[4];
				$gff_CDSlength{$temp[0]} = abs($temp[4]-$temp[3]);
				$gff_UTRlength{$temp[0]} = $gff_genelength{$temp[0]}-abs($temp[4]-$temp[3]);
			}	
		}
	}
}	

close DATAgff;

# OK now all the gff data are loaded into hashes
my $begin;
my $end;
my $CDS_repeats=0;
my $UTR_repeats=0;
# open blast results
# for octomys
open (DATAblast, "/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa_blasted_to_octomys_transcriptiome") or die "Failed to open laevis Blast results";
# for tympa
#open (DATAblast, "/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/repArc_kmer_31/velvet_repeat_lib/contigs.fa_blasted_to_tympa_transcriptiom") or die "Failed to open laevis Blast results";

my %blast_results;

while ( my $line = <DATAblast>) {
	@temp = split("\t",$line);
	if($temp[8] < $temp[9]){
		$begin=$temp[8];
		$end=$temp[9];
	}
	else{
		$begin=$temp[9];
		$end=$temp[8];
	}
	# RNAseq assembly id is $temp[1]
	# check what portion of the gene this hits. 4 scenarios
	# 1. all in the CDS
	# 2. on one end of the CDS plus UTR
	# 3. on the other end of the CDS plus UTR
	print $temp[1],"n";
	if(defined($gff_CDS{$temp[1]}[0])){
		if((($begin < $gff_CDS{$temp[1]}[0])&&($end < $gff_CDS{$temp[1]}[1]))||(($begin > $gff_CDS{$temp[1]}[0])&&($end > $gff_CDS{$temp[1]}[1]))){
			# the repeat is either higher or lower than the CDS; it is entirely in one or the other UTR
			$UTR_repeats += ($end-$begin);
		}
		elsif(($begin > $gff_CDS{$temp[1]}[0])&&($end < $gff_CDS{$temp[1]}[1])){
			# the repeat within the CDS; it is entirely in one or the other UTR
			$CDS_repeats += ($end-$begin);
		}
		elsif(($begin < $gff_CDS{$temp[1]}[0])&&($end > $gff_CDS{$temp[1]}[0])&&($end < $gff_CDS{$temp[1]}[1])){
			# the repeat spans one UTR and the CDS
			$UTR_repeats += ($gff_CDS{$temp[1]}[0]-$begin);
			$CDS_repeats += ($end-$gff_CDS{$temp[1]}[0]);
		}
		elsif(($begin > $gff_CDS{$temp[1]}[0])&&($begin < $gff_CDS{$temp[1]}[1])&&($end > $gff_CDS{$temp[1]}[1])){
			# the repeat spans the CDS and one UTR
			$CDS_repeats += ($gff_CDS{$temp[1]}[1]-$begin);		
			$UTR_repeats += ($end-$gff_CDS{$temp[1]}[1]);
		}
		elsif(($begin < $gff_CDS{$temp[1]}[0])&&($end > $gff_CDS{$temp[1]}[1])){
			# the repeat includes the CDS and portions of both UTRs
			$UTR_repeats += ($gff_CDS{$temp[1]}[0]-$begin)+($end-$gff_CDS{$temp[1]}[1]);		
			$CDS_repeats += ($gff_CDS{$temp[1]}[1]-$gff_CDS{$temp[1]}[0]);
		}
	}
}	

close DATAblast;

my $total_length;
my $counter=0;
# calculate the mean and total length
foreach my $key (sort(keys %gff_genelength)) {
	$total_length+=$gff_genelength{$key};
	$counter+=1;
}

print "The total number of base pairs in the transcriptome is ",$total_length,"\n";
print "The average transcriptome length is ",$total_length/$counter,"\n";

my $CDScounter=0;
my $CDStotal_length=0;
# calculate the mean and total length of CDS
foreach my $key (sort(keys %gff_CDSlength)) {
	$CDStotal_length+=$gff_CDSlength{$key};
	$CDScounter+=1;
}

my $UTRcounter=0;
my $UTRtotal_length=0;
# calculate the mean and total length of CDS
foreach my $key (sort(keys %gff_UTRlength)) {
	$UTRtotal_length+=$gff_UTRlength{$key};
	$UTRcounter+=1;
}

print "The total number of CDS base pairs in the transcriptome is ",$CDStotal_length,"\n";
print "The average CDS length is ",$CDStotal_length/$CDScounter,"\n";
print "The proportion of CDS bp that is repeats is ",$CDS_repeats/$CDStotal_length,"\n";

print "The total number of UTR base pairs in the transcriptome is ",$UTRtotal_length,"\n";
print "The average UTR length is ",$UTRtotal_length/$UTRcounter,"\n";
print "The proportion of UTR bp that is repeats is ",$UTR_repeats/$UTRtotal_length,"\n";
```

# Examining the incidence of repetitive elements in the Genome Assembly

Here is the repeat element assembly:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa
```
Here is the results of the blast output
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/repArc_AO248_kmer_31_to_AO248-scaffolds.fa_blastable
```

So, how many times does the most repetitive element show up in the assembly?  First find the coverage of the 10 most repetive elements like this:
```
grep -o -P '(?<=cov_).*(?=)' /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa | sort -rn | head -n 10
```
```
82.937500
82.903229
77.638885
74.660004
73.687500
73.628571
73.393936
72.500000
71.366669
70.621620
```

Now, get the header of the highest coverage one.
```
grep '82.937500' /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/repArc_kmer_31/velvet_repeat_lib/contigs.fa
```
```
>NODE_251811_length_32_cov_82.937500
```
print out a fasta file
```
awk -v seq="NODE_251811_length_32_cov_82.937500" -v RS='>' '$1 == seq {print RS $0}' /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/repArc_kmer_31/velvet_repeat_lib/contigs.fa
```
