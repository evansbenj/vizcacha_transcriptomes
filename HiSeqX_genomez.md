# HiSeqX Genomez

OK, we now have HiSeqX data from Tympa AO245 and Octomys AO248. For Tympa, there are 66 Gb of zipped data for its ~8.4 pg genome and for Octomys there are 32 Gb of zipped data for its ~4pg genome. So there should be similar coverage for each one.

# FastQC

This hopefully will help identify adapter sequences (might also pick up some TE seqs)


Wilson says that these adaptor sequences were used:

AO245 - GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG (TruSeq Adapter, Index 22)
AO248 - GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG (TruSeq Adapter, Index 25)
Both - AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT (TruSeq Universal Adapter)

So I made this adapter file using these plus the overrepresented seqs from fastqc:
>seq1/1
ATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATG
>seq2/1
ATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGT
>seq3/1
ATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATG
>seq4/1
ATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGT
>TruSeqAdapterIndex22/1
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG
>TruSeqAdapterIndex25/1
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG
>TruSeqUniversalAdapter/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT



# Trimming

Trim these data like this:
```
java -jar /home/ben/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog AO245.txt AO245_S7_L007_R1_001.fastq.gz AO245_S7_L007_R2_001.fastq.gz AO245_R1_trim_paired.fastq.gz AO245_R1_trim_single.fastq.gz AO245_R2_trim_paired.fastq.gz AO245_R2_trim_single.fastq.gz ILLUMINACLIP:/home/ben/Trimmomatic-0.36/adapters/HiSeqX_overrep.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36

java -jar /home/ben/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog AO248.txt AO248_S8_L007_R1_001.fastq.gz AO248_S8_L007_R2_001.fastq.gz AO248_R1_trim_paired.fastq.gz AO248_R1_trim_single.fastq.gz AO248_R2_trim_paired.fastq.gz AO248_R2_trim_single.fastq.gz ILLUMINACLIP:/home/ben/Trimmomatic-0.36/adapters/HiSeqX_overrep.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
```
