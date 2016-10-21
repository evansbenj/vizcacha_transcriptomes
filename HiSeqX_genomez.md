# HiSeqX Genomez

OK, we now have HiSeqX data from Tympa AO245 and Octomys AO248. For Tympa, there are 66 Gb of zipped data for its ~8.4 pg genome and for Octomys there are 32 Gb of zipped data for its ~4pg genome. So there should be similar coverage for each one.

# FastQC

This hopefully will help identify adapter sequences (might also pick up some TE seqs)


# Trimming

Trim these data like this:

java -jar /home/ben/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog AO245.txt AO245_S7_L007_R1_001.fastq.gz AO245_S7_L007_R2_001.fastq.gz AO245_R1_trim_paired.fastq.gz AO245_R1_trim_single.fastq.gz AO245_R2_trim_paired.fastq.gz AO245_R2_trim_single.fastq.gz ILLUMINACLIP:/home/ben/Trimmomatic-0.36/adapters/TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
