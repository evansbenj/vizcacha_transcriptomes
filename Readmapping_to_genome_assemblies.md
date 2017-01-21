# Readmapping to genome assemblies

Prepare genomes for mapping using bwa, samtools and piccard as previously.  To map paired end reads, I did this in this directory ofr octomys:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly
```
```
/usr/local/bin/bwa mem -M -t 16 AO248-scaffolds.fa ../AO248_R1_trim_paired.cor.fastq.gz ../AO248_R2_trim_paired.fastq.gz > octomys_WGS_to_genome_aln.sam
```

and later bits for tympa:

```
java -jar ~/picard-tools-1.131/picard.jar AddOrReplaceReadGroups \
      I=tymp_WGS_to_genome_aln_sorted.bam \
      O=tymp_WGS_to_genome_aln_sorted_rg.bam \
      RGID=tympaWGS \
      RGLB=tympaWGS \
      RGPL=illumina \
      RGPU=tympaWGS \
      RGSM=tympaWGS


/usr/local/bin/samtools index tymp_WGS_to_genome_aln_sorted_rg.bam

java -Xmx5G -jar ~/picard-tools-1.131/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=tymp_WGS_to_genome_aln_sorted_rg.bam OUTPUT=tymp_WGS_to_genome_aln_sorted_rg_dedup.bam METRICS_FILE=tymp_WGS_to_genome_aln_sorted_rg_dedup_metrics.txt

~/samtools_2016/bin/samtools mpileup -d8000 -ugf AO245-scaffolds.fa -t DP,AD tymp_WGS_to_genome_aln_sorted_rg_dedup.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o tymp_WGS_to_genome_aln_sorted_rg_dedup.bam.vcf.gz
```
