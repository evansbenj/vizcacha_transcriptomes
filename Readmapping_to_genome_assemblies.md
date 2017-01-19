# Readmapping to genome assemblies

Prepare genomes for mapping using bwa, samtools and piccard as previously.  To map paired end reads, I did this in this directory ofr octomys:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly
```
```
/usr/local/bin/bwa mem -M -t 16 AO248-scaffolds.fa ../AO248_R1_trim_paired.cor.fastq.gz ../AO248_R2_trim_paired.fastq.gz > octomys_WGS_to_genome_aln.sam
```
