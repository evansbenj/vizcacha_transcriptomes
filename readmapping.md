# Mapping reads to the Trinity assembly

The goal here is to map reads and call genotypes in Tympa and Octomys to see if there are more heterozygous sites in Tympa due to chimerial contigs comprised of homeologous reads.

To prepare the genome on info (I am using the unique reads recovered from cdhit), for tympa:
```
/usr/local/bin/bwa index -a bwtsw Tympa_all_transcriptomes_assembled_together_unique.fasta
/usr/local/bin/samtools faidx Tympa_all_transcriptomes_assembled_together_unique.fasta
java -jar ~/picard-tools-1.131/picard.jar CreateSequenceDictionary REFERENCE=Tympa_all_transcriptomes_assembled_together_unique.fasta OUTPUT=Tympa_all_transcriptomes_assembled_together_unique.dict
```
```
/usr/local/bin/bwa mem -M -t 16 /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_all_R1_trim_paired.fastq.gz /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_all_R2_trim_paired.fastq.gz > /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln.sam
```
To prepare the genome on info (I am using the unique reads recovered from cdhit), for octomys:
```
/usr/local/bin/bwa index -a bwtsw Octomys_all_transcriptomes_assembled_together_unique.fasta
/usr/local/bin/samtools faidx Octomys_all_transcriptomes_assembled_together_unique.fasta
java -jar ~/picard-tools-1.131/picard.jar CreateSequenceDictionary REFERENCE=Octomys_all_transcriptomes_assembled_together_unique.fasta OUTPUT=Octomys_all_transcriptomes_assembled_together_unique.dict
```
```
/usr/local/bin/bwa mem -M -t 16 /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/AO248_all_R1_trim_paired.fastq.gz /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/AO248_all_R2_trim_paired.fastq.gz > /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln.sam
```
make bam files
```
/usr/local/bin/samtools view -bt /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta -o /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln.bam /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln.sam

/usr/local/bin/samtools view -bt /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta -o /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln.bam /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln.sam

```
sort the bam files

```
/usr/local/bin/samtools sort /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln.bam /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted

/usr/local/bin/samtools sort /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln.bam /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted
```
make a bai file

```
/usr/local/bin/samtools index /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted.bam

/usr/local/bin/samtools index /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted.bam
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
      
      
```


Use unified genotyper to call bases

```
java -Xmx4G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta -I /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_aln_sorted_rg.bam -out_mode EMIT_ALL_CONFIDENT_SITES -o /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/tympa_allconfident.vcf


java -Xmx4G -jar /usr/local/gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta -I /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_aln_sorted_rg.bam -out_mode EMIT_ALL_CONFIDENT_SITES -o /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/octomys_allconfident.vcf
```
