# Readmapping to genome assemblies

Prepare genomes for mapping using bwa, samtools and piccard as previously.  To map paired end reads, I did this in this directory ofr octomys:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly
```

prepare the genomez
```
/usr/local/bin/bwa index -a bwtsw AO248_newtrim_scaffolds.fa
/usr/local/bin/samtools faidx AO248_newtrim_scaffolds.fa
~/jre1.8.0_111/bin/java -Xmx2g -jar ~/picard-tools-1.131/picard.jar CreateSequenceDictionary REFERENCE=AO248_newtrim_scaffolds.fa OUTPUT=AO248_newtrim_scaffolds.dict
```

```
/usr/local/bin/bwa index -a bwtsw AO245_newtrim_scaffolds.fa
/usr/local/bin/samtools faidx AO245_newtrim_scaffolds.fa
~/jre1.8.0_111/bin/java -Xmx2g -jar ~/picard-tools-1.131/picard.jar CreateSequenceDictionary REFERENCE=AO245_newtrim_scaffolds.fa OUTPUT=AO245_newtrim_scaffolds.dict
```

mapping reads

```
/usr/local/bin/bwa mem -M -t 16 AO248_newtrim_scaffolds.fa ../AO248_R1_newtrim_paired.cor.fastq.gz ../AO248_R2_newtrim_paired.fastq.gz > octomys_WGS_to_newgenome_aln.sam
```


```
/usr/local/bin/bwa mem -M -t 16 AO245_newtrim_scaffolds.fa ../AO245_R1_newtrim_paired.cor.fastq.gz ../AO245_R2_newtrim_paired.fastq.gz > tympa_WGS_to_newgenome_aln.sam
```

This caused problems with the octomys file due to differently named paired end reads.  I suspect this is related to quake.  I corrected these reads like this:

```
../../bbmap/bbmap/repair.sh in1=../AO248_R1_newtrim_paired.cor.fastq.gz in2=../AO248_R2_newtrim_paired.cor.fastq.gz out1=fixed1.fq out2=fixed2.fq outsingle=single.fq
```
and then re-ran the ampping using the `fixed1.fq` and `fixed2.fq` files. This worked even though none of the paired reads did not have a mate.


make bam files

```
/usr/local/bin/samtools view -bt AO248_newtrim_scaffolds.fa -o octomys_WGS_to_newgenome_aln.bam octomys_WGS_to_newgenome_aln.sam
```

```
/usr/local/bin/samtools view -bt AO248_newtrim_scaffolds.fa -o tympa_WGS_to_newgenome_aln.bam tympa_WGS_to_newgenome_aln.sam
```

delete sam files

```
rm -f octomys_WGS_to_newgenome_aln.sam
```
```
rm -f tympa_WGS_to_newgenome_aln.sam
```

sort bam files
```
/usr/local/bin/samtools sort octomys_WGS_to_newgenome_aln.bam -o octomys_WGS_to_newgenome_aln_sorted
```
```
/usr/local/bin/samtools sort tympa_WGS_to_newgenome_aln.bam -o tympa_WGS_to_newgenome_aln_sorted
```

add readgroups:

```
java -jar ~/picard-tools-1.131/picard.jar AddOrReplaceReadGroups \
      I=octomys_WGS_to_newgenome_aln_sorted.bam \
      O=octomys_WGS_to_newgenome_aln_sorted_rg.bam \
      RGID=octWGS \
      RGLB=octWGS \
      RGPL=illumina \
      RGPU=octWGS \
      RGSM=octWGS
```

```
java -jar ~/picard-tools-1.131/picard.jar AddOrReplaceReadGroups \
      I=tympa_WGS_to_newgenome_aln_sorted.bam \
      O=tympa_WGS_to_newgenome_aln_sorted_rg.bam \
      RGID=tympaWGS \
      RGLB=tympaWGS \
      RGPL=illumina \
      RGPU=tympaWGS \
      RGSM=tympaWGS
```
readgrooup did not work, but this doesn't matter for genotyping without BQSR with samtools.

make a bai file
```
/usr/local/bin/samtools index octomys_WGS_to_newgenome_aln_sorted.bam
```
```
/usr/local/bin/samtools index tympa_WGS_to_newgenome_aln_sorted.bam
```
Mark duplicates

```
java -Xmx5G -jar ~/picard-tools-1.131/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=octomys_WGS_to_newgenome_aln_sorted.bam OUTPUT=octomys_WGS_to_newgenome_aln_sorted_dedup.bam METRICS_FILE=octomys_WGS_to_newgenome_aln_sorted_dedup_metrics.txt
```
```
java -Xmx5G -jar ~/picard-tools-1.131/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=tympa_WGS_to_newgenome_aln_sorted.bam OUTPUT=tympa_WGS_to_newgenome_aln_sorted_dedup.bam METRICS_FILE=tympa_WGS_to_newgenome_aln_sorted_dedup_metrics.txt
```

picard is not working because of java issues so I used samtools:
```
samtools rmdup octomys_WGS_to_newgenome_aln_sorted.bam octomys_WGS_to_newgenome_aln_sorted_dedup.bam
```


Use samtools and bcftools to call genotypes and filter


```
~/samtools_2016/bin/samtools mpileup -d8000 -ugf AO248_newtrim_scaffolds.fa -t DP,AD octomys_WGS_to_newgenome_aln_sorted_dedup.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o oct_WGS_to_newgenome_aln_sorted_dedup.bam.vcf.gz
```

```
~/samtools_2016/bin/samtools mpileup -d8000 -ugf AO245_newtrim_scaffolds.fa -t DP,AD tymp_WGS_to_newgenome_aln_sorted_dedup.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o tymp_WGS_to_newgenome_aln_sorted_dedup.bam.vcf.gz
```

# Coverage

```
samtools depth  oct_WGS_to_newgenome_aln_sorted_dedup.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'

```
```
samtools depth  tymp_WGS_to_newgenome_aln_sorted_dedup.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
```

otherstuff
```

/usr/local/bin/samtools index tymp_WGS_to_newgenome_aln_sorted_rg.bam

java -Xmx5G -jar ~/picard-tools-1.131/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=tymp_WGS_to_newgenome_aln_sorted_rg.bam OUTPUT=tymp_WGS_to_newgenome_aln_sorted_rg_dedup.bam METRICS_FILE=tymp_WGS_to_newgenome_aln_sorted_rg_dedup_metrics.txt

~/samtools_2016/bin/samtools mpileup -d8000 -ugf AO245_newtrim_scaffolds.fa -t DP,AD tymp_WGS_to_newgenome_aln_sorted_rg_dedup.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o tymp_WGS_to_newgenome_aln_sorted_rg_dedup.bam.vcf.gz
```

# Genotyping
```
samtools mpileup -d8000 -ugf ../xenXL_MT.fasta -t DP,AD BMNH1947_2_24_78_stampy_to_XLmtDNA_sorted.bam | bcftools call -V indels --format-fields GQ -m -g 1 -O z -o BMNH1947_2_24_78_stampy_to_XLmtDNA.bam.g.vcf.gz
```
```
samtools mpileup -d8000 -ugf ../xenXL_MT.fasta -t DP,AD BMNH1947_2_24_79_stampy_to_XLmtDNA_sorted.bam | bcftools call -V indels --format-fields GQ -m -g 1 -O z -o BMNH1947_2_24_79_stampy_to_XLmtDNA.bam.g.vcf.gz
```
```
samtools mpileup -d8000 -ugf ../xenXL_MT.fasta -t DP,AD 16294_stampy_to_XLmtDNA_sorted.bam | bcftools call -V indels --format-fields GQ -m -g 1 -O z -o 16294_stampy_to_XLmtDNA.bam.g.vcf.gz
```
```
tabix -p vcf myvcf.vcf.gz
```
```
java -Xmx8G -cp /mnt/expressions/ben_evans/bin/GenomeAnalysisTK-nightly-2017-10-07-g1994025/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R ../xenXL_MT.fasta -V BMNH1947_2_24_78_stampy_to_XLmtDNA.bam.g.vcf.gz -V BMNH1947_2_24_79_stampy_to_XLmtDNA.bam.g.vcf.gz -V 16294_stampy_to_XLmtDNA.bam.g.vcf.gz -out ancientfrogz.g.vcf.gz -assumeSorted
```
extract allele depth column:
```
vcftools --gzvcf Merged.vcf.gz --extract-FORMAT-info AD
```
to get only the transcripts that are sex-biased:
```
grep -f SL_trans_IDs_sig_Sex_biased Merged.vcf.gz_out.AD.FORMAT > SL_allelic_depth_sig_sex_biased.AD.FORMAT
```
where `SL_trans_IDs_sig_Sex_biased` is from 
```
XB_adult_liver_SL_transcripts_sig_Sex_biased <- borTad_laevisGenome_deseq2_tmm_combine_st46_chr8L %>% filter((start< 54000000)&((logFC< -2)|(logFC> 2))&(padj<0.01))
library(data.table) # install if not installed already
fwrite(list(XB_adult_liver_SL_transcripts$trans_id), file = "SL_trans_IDs")
```
I wrote this script to try to parse the results, but the finds were not impressive because there were very few SNPs
```
#!/usr/bin/env perl
use strict;
#use warnings;
#no warnings qw(numeric);

#  This program reads in allele depths of male and female individuals
# For the current example, I have a mother and a father in the first two
# columns, and then 3 daughters and three sons.  

# I want to identify transcripts in the sex-linked region that have SNPs
# that permit me to quantify W- and Z- expression

#mother		dad		daughters		sons
#X,0		X,Y		X,0 or Y,0		X,Y or X,0 or 0,Y	Z-linked SNP from father 
#X,Y		X,0		X,0     		X,Y					Z-linked SNP from mother on biparental (W&Z) transcript 
#0,Y		X,0		X,0     		X,Y					Z-linked SNP from mother on Z-specific transcript
#X,0		0,0		X,0				0,0					W-linked SNP from mother on W-specific transcript 
#X,Y		X,0		X,Y				X,0					W-linked SNP from mother on biparental (W&Z) transcript 

#X,0		X,Y		X,Y				X,Y					paternal SNP on autosomes
#X,Y		X,0		X,Y				X,Y					maternal SNP on autosomes

# Some assumptions: if no genotype is present in dad and sons, we cant
# tell whether alleles expressed in daughters are W- or Z-



# I outputted a list of transIDs in the sex linked region like this:
# library(data.table) # install if not installed already
# fwrite(list(XB_adult_liver_SL_transcripts$trans_id), file = "SL_trans_IDs")

# then I can use this to grep the allelic depth for these specific transcripts
# grep -f SL_trans_IDs_XB_adult_liver Merged.vcf.gz_out.AD.FORMAT > SL_allelic_depth


my $inputfile = $ARGV[0]; # file name
my $outputfile = $ARGV[1]; # output file name


unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}


unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

my @temp = ();
my @mom = ();
my @dad = ();
my @girl1 = ();
my @girl2 = ();
my @girl3 = ();
my @boy1 = ();
my @boy2 = ();
my @boy3 = ();
my $sum_girls_allele1;
my $sum_girls_allele2;
my $sum_boys_allele1;
my $sum_boys_allele2;
my $num_girl_genotypes = 0;
my $num_boy_genotypes = 0;
my %transcripts = ();

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'CHROM'){
		@mom=split(',',$temp[2]);
		@dad=split(',',$temp[3]);
		@girl1=split(',',$temp[4]);
		@girl2=split(',',$temp[5]);
		@girl3=split(',',$temp[6]);
		@boy1=split(',',$temp[7]);
		@boy2=split(',',$temp[8]);
		@boy3=split(',',$temp[9]);
		# sum the alleles for available genotypes
		$sum_girls_allele1=0;
		$sum_girls_allele2=0;
		$sum_boys_allele1=0;
		$sum_boys_allele2=0;
		$num_girl_genotypes=0;
		$num_boy_genotypes=0;
		if($girl1[0] ne "."){
			$sum_girls_allele1+=$girl1[0];
			$sum_girls_allele2+=$girl1[1];
			$num_girl_genotypes+=1;
		}
		if($girl2[0] ne "."){
			$sum_girls_allele1+=$girl2[0];
			$sum_girls_allele2+=$girl2[1];
			$num_girl_genotypes+=1;
		}
		if($girl3[0] ne "."){
			$sum_girls_allele1+=$girl3[0];
			$sum_girls_allele2+=$girl3[1];
			$num_girl_genotypes+=1;
		}
		if($boy1[0] ne "."){
			$sum_boys_allele1+=$boy1[0];
			$sum_boys_allele2+=$boy1[1];
			$num_boy_genotypes+=1;
		}
		if($boy1[0] ne "."){
			$sum_boys_allele1+=$boy2[0];
			$sum_boys_allele2+=$boy2[1];
			$num_boy_genotypes+=1;
		}
		if($boy1[0] ne "."){
			$sum_boys_allele1+=$boy3[0];
			$sum_boys_allele2+=$boy3[1];
			$num_boy_genotypes+=1;
		}
		if((($mom[0] ne ".")&&($mom[0] != 0)&&($mom[1] != 0)&&($dad[0] != ".")&&
			($dad[0] != 0)&&($dad[1] != 0))||($mom[0] ne ".")&&($dad[0] ne ".")){
			# there is not a sex-specific SNP or no genotype in both parents
		#	print "no sex-specific SNP or no genotype in both parents $line\n";
		}
		elsif(($#mom > 1)||($#dad > 1)||($#girl1 > 1)||($#girl2 > 1)||
			($#girl3 > 1)||($#boy1 > 1)||($#boy2 > 1)||($#boy3 > 1)){
		#		print "someone has a triallelic SNP $line\n";
		}
		elsif(($mom[0] ne ".")&&($mom[0] != 0)&&($mom[1] != 0)){
			# there is a SNP in mom 
			# this is the first class of SNPs we care about (maternal W- or Z- SNPs)
			# we can use this to identity W- and Z- expression in daughters
			if(($dad[0] == 0)||($dad[0] ne ".")){
				# the dad is homoz for the second allele or a missing genotype
				# check if all the daughters are missing the first allele
				if(($sum_girls_allele1 == 0)&&($num_girl_genotypes > 0)){
					# allele 1 is a maternal Z-SNP that was not inherited by any of the girls
					print "allele 1 M Z-SNP $line\n";
					# later I can print out the expression of each allele in girls and boys
				}
				elsif(($sum_girls_allele1 != 0)&&($num_girl_genotypes > 0)){
					# allele 1 may be a maternal W-SNP that was only inherited by the girls
					if(($sum_boys_allele1 == 0)&&($num_boy_genotypes>0)){
						print "allele 1 M W-SNP $line\n";
						# later I can print out the expression of each allele in girls and boys
					}
					elsif(($sum_boys_allele1 != 0)&&($num_boy_genotypes > 0)&&
						(($sum_girls_allele2 != 0)||($num_girl_genotypes == 0))){
						# if we see expression of both alleles in both sexes, this must be a double SNP
						print "Double SNP $line\n";
						# later I can print out the expression of each allele in girls and boys
					}
					elsif(($sum_boys_allele1 > 0)&&($num_boy_genotypes>0)&&
						(($sum_girls_allele2 == 0)||($num_girl_genotypes == 0))){
						# if we see expression of a M SNP in boys and not girls, this must be a M Z-SNP
						print "allele 1 M Z-SNP $line\n";
						# later I can print out the expression of each allele in girls and boys
					}
					elsif(($sum_boys_allele1 == 0)&&($num_boy_genotypes == 0)){
						# if we see expression in boys, this must be a Z-SNP
						print "M-SNP both alleles girl-specific $line\n";
						# later I can print out the expression of each allele in girls and boys
					}
				}
				elsif(($num_girl_genotypes == 0)&&($num_boy_genotypes == 0)){
				#	print "there were missing genotypes in all offspring of at least one sex $line\n";
				}
			}
			elsif(($dad[1] == 0)||($dad[0] ne ".")){
				# the dad is homoz for the first allele or a missing genotype
				# check if all the daughters are missing the second allele
				if(($sum_girls_allele2 == 0)&&($num_girl_genotypes > 0)){
					# allele 2 is a maternal Z-SNP that was not inherited by any of the girls
					print "allele 2 M Z-SNP $line\n";
					# later I can print out the expression of each allele in girls and boys
				}
				elsif(($sum_girls_allele2 != 0)&&($num_girl_genotypes > 0)){
					# allele 1 may be a maternal W-SNP that was only inherited by the girls
					# allele 1 may be a maternal W-SNP that was only inherited by the girls
					if(($sum_boys_allele2 == 0)&&($num_boy_genotypes>0)){
						print "allele 2 M W-SNP $line\n";
						# later I can print out the expression of each allele in girls and boys
					}
					elsif(($sum_boys_allele2 != 0)&&($num_boy_genotypes > 0)&&
						(($sum_girls_allele1 != 0)||($num_girl_genotypes == 0))){
						# if we see expression of both alleles in both sexes, this must be a double SNP
						print "Double SNP $line\n";
						# later I can print out the expression of each allele in girls and boys
					}
					elsif(($sum_boys_allele2 > 0)&&($num_boy_genotypes>0)&&
						(($sum_girls_allele1 == 0)||($num_girl_genotypes == 0))){
						# if we see expression of a M SNP in boys and not girls, this must be a M Z-SNP
						print "allele 2 M Z-SNP $line\n";
						# later I can print out the expression of each allele in girls and boys
					}
					elsif(($sum_boys_allele2 == 0)&&($num_boy_genotypes == 0)){
						# if we see expression in boys, this must be a Z-SNP
						print "M-SNP both alleles girl-specific $line\n";
						# later I can print out the expression of each allele in girls and boys
					}
				}
				elsif(($num_girl_genotypes == 0)&&($num_boy_genotypes == 0)){
				#	print "there were missing genotypes in all offspring of at least one sex $line\n";
				}
			}	
		}		
		elsif(($dad[0] ne ".")&&($dad[0] != 0)&&($dad[1] != 0)){ 
			# there is a SNP in dad and the mom has a homozygous or missing genotype
			# this isn't very useful unless (which is unlikely) the mother is fixed
			# for a different SNP.  In this case one could distinguish between W and Z
			# in daughters.  Otherwise, pretty useless.
		}		
		elsif(($mom[0] ne ".")&&($mom[0] != 0)&&($mom[1] != 0)&&
			($dad[0] ne ".")&&($dad[0] != 0)&&($dad[1] != 0)){ 
				print "Double SNP $line";
		}
	}
}

close OUTFILE1;
```
