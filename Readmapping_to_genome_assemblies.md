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

Use samtools and bcftools to call genotypes and filter


```
~/samtools_2016/bin/samtools mpileup -d8000 -ugf AO248_newtrim_scaffolds.fa -t DP,AD octomys_WGS_to_newgenome_aln_sorted_dedup.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o oct_WGS_to_newgenome_aln_sorted_dedup.bam.vcf.gz
```

```
~/samtools_2016/bin/samtools mpileup -d8000 -ugf AO245_newtrim_scaffolds.fa -t DP,AD tymp_WGS_to_newgenome_aln_sorted_dedup.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o tymp_WGS_to_newgenome_aln_sorted_dedup.bam.vcf.gz
```

othetstuff
```

/usr/local/bin/samtools index tymp_WGS_to_newgenome_aln_sorted_rg.bam

java -Xmx5G -jar ~/picard-tools-1.131/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=tymp_WGS_to_newgenome_aln_sorted_rg.bam OUTPUT=tymp_WGS_to_newgenome_aln_sorted_rg_dedup.bam METRICS_FILE=tymp_WGS_to_newgenome_aln_sorted_rg_dedup_metrics.txt

~/samtools_2016/bin/samtools mpileup -d8000 -ugf AO245_newtrim_scaffolds.fa -t DP,AD tymp_WGS_to_newgenome_aln_sorted_rg_dedup.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o tymp_WGS_to_newgenome_aln_sorted_rg_dedup.bam.vcf.gz
```
