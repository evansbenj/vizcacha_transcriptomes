# Mapping reads to the Trinity assembly

The goal here is to map reads and call genotypes in Tympa and Octomys to see if there are more heterozygous sites in Tympa due to chimerial contigs comprised of homeologous reads. I also want to map the tympa reads to octomys and vice versa.

To prepare the genome on info (I am using the unique reads recovered from cdhit), for tympa:
```
/usr/local/bin/bwa index -a bwtsw Tympa_all_transcriptomes_assembled_together_unique.fasta
/usr/local/bin/samtools faidx Tympa_all_transcriptomes_assembled_together_unique.fasta
java -jar ~/picard-tools-1.131/picard.jar CreateSequenceDictionary REFERENCE=Tympa_all_transcriptomes_assembled_together_unique.fasta OUTPUT=Tympa_all_transcriptomes_assembled_together_unique.dict
```

mapping tympareads to tympaRNAsequnique:
```
/usr/local/bin/bwa mem -M -t 16 /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_all_R1_trim_paired.fastq.gz /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_all_R2_trim_paired.fastq.gz > /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln.sam
```
To prepare the genome on info (I am using the unique reads recovered from cdhit), for octomys:
```
/usr/local/bin/bwa index -a bwtsw Octomys_all_transcriptomes_assembled_together_unique.fasta
/usr/local/bin/samtools faidx Octomys_all_transcriptomes_assembled_together_unique.fasta
java -jar ~/picard-tools-1.131/picard.jar CreateSequenceDictionary REFERENCE=Octomys_all_transcriptomes_assembled_together_unique.fasta OUTPUT=Octomys_all_transcriptomes_assembled_together_unique.dict
```
mapping ocromysreads to octoRNAsequnique
```
/usr/local/bin/bwa mem -M -t 16 /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/AO248_all_R1_trim_paired.fastq.gz /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/AO248_all_R2_trim_paired.fastq.gz > /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln.sam
```

mapping tympareads to octoRNAsequnique
```
/usr/local/bin/bwa mem -M -t 16 /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_all_R1_trim_paired.fastq.gz /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_all_R2_trim_paired.fastq.gz > /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympareads_mappedto_octomysassembly.sam
```
mapping ocromysreads to tympaRNAsequnique
```
/usr/local/bin/bwa mem -M -t 16 /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/AO248_all_R1_trim_paired.fastq.gz /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/AO248_all_R2_trim_paired.fastq.gz > /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octreads_mappedto_tympsassembly.sam
```


make bam files
```
/usr/local/bin/samtools view -bt /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta -o /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln.bam /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln.sam

/usr/local/bin/samtools view -bt /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta -o /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln.bam /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln.sam

/usr/local/bin/samtools view -bt /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta -o /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octreads_mappedto_tympsassembly.bam /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octreads_mappedto_tympsassembly.sam

/usr/local/bin/samtools view -bt /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta -o /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympareads_mappedto_octomysassembly.bam /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympareads_mappedto_octomysassembly.sam
```

delete the sam files

```
rm /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octreads_mappedto_tympsassembly.sam

rm /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympareads_mappedto_octomysassembly.sam


```
sort the bam files

```
/usr/local/bin/samtools sort /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln.bam /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted

/usr/local/bin/samtools sort /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln.bam /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted

/usr/local/bin/samtools sort /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympareads_mappedto_octomysassembly.bam /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympareads_mappedto_octomysassembly_sorted

/usr/local/bin/samtools sort /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octreads_mappedto_tympsassembly.bam /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octreads_mappedto_tympsassembly_sorted

```

add a readgroup with picard
```
java -jar ~/picard-tools-1.131/picard.jar AddOrReplaceReadGroups \
      I=/home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted.bam \
      O=/home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted_rg.bam \
      RGID=tympa_combined \
      RGLB=tympa_combined \
      RGPL=illumina \
      RGPU=tympa_combined \
      RGSM=tympa_combined
      
java -jar ~/picard-tools-1.131/picard.jar AddOrReplaceReadGroups \
      I=/home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted.bam \
      O=/home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted_rg.bam \
      RGID=4 \
      RGLB=oct_combined \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=oct_combined      
      
      
java -jar ~/picard-tools-1.131/picard.jar AddOrReplaceReadGroups \
      I=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octreads_mappedto_tympsassembly_sorted.bam \
      O=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octreads_mappedto_tympsassembly_sorted_rg.bam \
      RGID=4 \
      RGLB=oct_to_tympassemb \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=oct_to_tympassemb 

java -jar ~/picard-tools-1.131/picard.jar AddOrReplaceReadGroups \
      I=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympareads_mappedto_octomysassembly_sorted.bam \
      O=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympareads_mappedto_octomysassembly_sorted_rg.bam \
      RGID=4 \
      RGLB=tymp_to_octassemb \
      RGPL=illumina \
      RGPU=unit1 \
      RGSM=tymp_to_octassemb 


```

make a bai file

```
/usr/local/bin/samtools index /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted_rg.bam

/usr/local/bin/samtools index /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted_rg.bam

/usr/local/bin/samtools index /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympareads_mappedto_octomysassembly_sorted_rg.bam

/usr/local/bin/samtools index /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octreads_mappedto_tympsassembly_sorted_rg.bam
```

Moving to this directory: `/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/`
Mark duplicates

```
java -Xmx5G -jar ~/picard-tools-1.131/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=/home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted_rg.bam OUTPUT=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octomys_aln_sorted_rg_dedup.bam METRICS_FILE=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octomys_aln_sorted_rg_dedup_metrics.txt

java -Xmx5G -jar ~/picard-tools-1.131/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=/home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted_rg.bam OUTPUT=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympa_aln_sorted_rg_dedup.bam METRICS_FILE=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympa_aln_sorted_rg_dedup_metrics.txt

java -Xmx5G -jar ~/picard-tools-1.131/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympareads_mappedto_octomysassembly_sorted_rg.bam OUTPUT=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympareads_mappedto_octomysassembly_sorted_rg_dedup.bam METRICS_FILE=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/tympareads_mappedto_octomysassembly_aln_sorted_rg_dedup_metrics.txt


java -Xmx5G -jar ~/picard-tools-1.131/picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 INPUT=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octreads_mappedto_tympsassembly_sorted_rg.bam OUTPUT=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octreads_mappedto_tympsassembly_sorted_rg_dedup.bam METRICS_FILE=/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_RNAseq/octreads_mappedto_tympsassembly_aln_sorted_rg_dedup_metrics.txt

```
change names
```
mv tympa_aln_sorted_rg_dedup.bam tympareads_mappedto_tympaassembly_sorted_rg_dedup.bam
mv octomys_aln_sorted_rg_dedup.bam octreads_mappedto_octomysassembly_sorted_rg_dedup.bam
```

Use unified genotyper to call bases

```
java -Xmx4G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta -I /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted_rg.bam -out_mode EMIT_ALL_CONFIDENT_SITES -o /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_allconfident.vcf


java -Xmx4G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta -I /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted_rg.bam -out_mode EMIT_ALL_CONFIDENT_SITES -o /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_allconfident.vcf
```
Indel realignment not needed because this is now done with HaplotypeCaller
Actually use haplotypecaller and then genotypegvcfs to call bases

```
java -Xmx4G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta -I /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted_rg.bam -out_mode EMIT_ALL_CONFIDENT_SITES -o /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_allconfident.vcf


java -Xmx4G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta -I /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted_rg.bam -out_mode EMIT_ALL_CONFIDENT_SITES -o /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_allconfident.vcf
```

Actually use samtools and bcftools to call genotypes and filter
```
~/samtools_2016/bin/samtools mpileup -d8000 -ugf /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta -t DP,AD octreads_mappedto_octomysassembly_sorted_rg_dedup.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o octreads_mappedto_octomysassembly_sorted_rg_dedup.bam.vcf.gz

~/samtools_2016/bin/samtools mpileup -d8000 -ugf /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta -t DP,AD tympareads_mappedto_tympaassembly_sorted_rg_dedup.bam  | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o tympareads_mappedto_tympaassembly_sorted_rg_dedup.bam.vcf.gz

~/samtools_2016/bin/samtools mpileup -d8000 -ugf /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta -t DP,AD octreads_mappedto_tympsassembly_sorted_rg_dedup.bam  | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o octreads_mappedto_tympsassembly_sorted_rg_dedup.bam.vcf.gz

~/samtools_2016/bin/samtools mpileup -d8000 -ugf /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta -t DP,AD tympareads_mappedto_octomysassembly_sorted_rg_dedup.bam | ~/samtools_2016/bcftools-1.3.1/bcftools call -V indels --format-fields GQ -m -O z | ~/samtools_2016/bcftools-1.3.1/bcftools filter -e 'FORMAT/GT = "." || FORMAT/DP < 10 || FORMAT/GQ < 20 || FORMAT/GQ = "."' -O z -o tympareads_mappedto_octomysassembly_sorted_rg_dedup.bam.vcf.gz

```

Make tab delimited files
```
~/tabix-0.2.6/bgzip tympa_allconfident.vcf
~/tabix-0.2.6/tabix -p vcf tympa_allconfident.vcf.gz
zcat tympa_allconfident.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > tympa_allconfident.vcf.gz.tab
```
and
```
~/tabix-0.2.6/bgzip octomys_allconfident.vcf
~/tabix-0.2.6/tabix -p vcf octomys_allconfident.vcf.gz
zcat octomys_allconfident.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > octomys_allconfident.vcf.gz.tab
```
