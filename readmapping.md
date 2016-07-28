# Mapping reads to the Trinity assembly

The goal here is to map reads and call genotypes in Tympa and Octomys to see if there are more heterozygous sites in Tympa due to chimerial contigs comprised of homeologous reads.

To prepare the genome on info (I am using the unique reads recovered from cdhit)
```
/usr/local/bin/bwa index -a bwtsw Tympa_all_transcriptomes_assembled_together_unique.fasta
/usr/local/bin/samtools faidx Tympa_all_transcriptomes_assembled_together_unique.fasta
java -jar ~/picard-tools-1.131/picard.jar CreateSequenceDictionary REFERENCE=Tympa_all_transcriptomes_assembled_together_unique.fasta OUTPUT=Tympa_all_transcriptomes_assembled_together_unique.dict
```
``
/apps/bwa/0.7.12/bwa aln reference_genome individual_1.fastq > individual_1.sai
/apps/bwa/0.7.12/bwa samse reference_genome.fa individual_1.sai individual_1.fastq > individual_1.sam
/apps/bwa/0.7.12/bwa samse -r "@RG\tID:FLOWCELL1.LANE6\tSM:Individual_1\tPL:illumina" reference.fa Individual_1.sai Individual_1.fastq > Individual_1.sam
/apps/samtools/0.1.19/samtools view -bt reference_genome -o Individual_1.bam Individual_1.sam
/apps/samtools/0.1.19/samtools sort Individual_1.bam Individual_1_sorted
```
