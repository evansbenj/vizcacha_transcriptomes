# HiSeqX Genomez

OK, we now have HiSeqX data from Tympa AO245 and Octomys AO248. For Tympa, there are 66 Gb of zipped data for its ~8.4 pg genome and for Octomys there are 32 Gb of zipped data for its ~4pg genome. So there should be similar coverage for each one.

# FastQC
This hopefully will help identify adapter sequences (might also pick up some TE seqs)


Wilson says that these adaptor sequences were used:

```
AO245 - GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG (TruSeq Adapter, Index 22)
AO248 - GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG (TruSeq Adapter, Index 25)
Both - AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT (TruSeq Universal Adapter)
```

So I made this adapter file using these plus the overrepresented seqs from fastqc:
```
>seq1/1
ATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATG
>seq2/1
ATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGT
>seq3/1
ATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATG
>TruSeqAdapterIndex22/1
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG
>TruSeqAdapterIndex25/1
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG
>TruSeqUniversalAdapter/1
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
```


# Update

New trimmomatic adapter file (`HiSeqX_overrep_new.fa`) includes the rev comp of each sequence and removes the /1 which specifies that reads only be checked in the forward read

```
>seq1
ATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATG
>seq2
ATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGT
>seq3
ATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATG
>TruSeqAdapterIndex22
GATCGGAAGAGCACACGTCTGAACTCCAGTCACGAGTGGATATCTCGTATGCCGTCTTCTGCTTG
>TruSeqAdapterIndex25
GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATATCTCGTATGCCGTCTTCTGCTTG
>TruSeqUniversalAdapter
AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
>seq1_revcomp
CATACGAGATATCCACTCGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT
>seq2_revcomp
ACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGAT
>seq3_revcomp
CATACGAGATATATCAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT
>TruSeqAdapterIndex22_revcomp
CAAGCAGAAGACGGCATACGAGATATCCACTCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
>TruSeqAdapterIndex25_revcomp
CAAGCAGAAGACGGCATACGAGATATATCAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATC
>TruSeqUniversalAdapter_revcomp
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
```

# Trimming

Trim these data like this:
```
java -jar /home/ben/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog AO245.txt AO245_S7_L007_R1_001.fastq.gz AO245_S7_L007_R2_001.fastq.gz AO245_R1_trim_paired.fastq.gz AO245_R1_trim_single.fastq.gz AO245_R2_trim_paired.fastq.gz AO245_R2_trim_single.fastq.gz ILLUMINACLIP:/home/ben/Trimmomatic-0.36/adapters/HiSeqX_overrep_new.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36

java -jar /home/ben/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog AO248.txt AO248_S8_L007_R1_001.fastq.gz AO248_S8_L007_R2_001.fastq.gz AO248_R1_trim_paired.fastq.gz AO248_R1_trim_single.fastq.gz AO248_R2_trim_paired.fastq.gz AO248_R2_trim_single.fastq.gz ILLUMINACLIP:/home/ben/Trimmomatic-0.36/adapters/HiSeqX_overrep_new.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
```

new trim

```
java -jar /home/ben/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog AO245new.txt AO245_S7_L007_R1_001.fastq.gz AO245_S7_L007_R2_001.fastq.gz AO245_R1_newtrim_paired.fastq.gz AO245_R1_newtrim_single.fastq.gz AO245_R2_newtrim_paired.fastq.gz AO245_R2_newtrim_single.fastq.gz ILLUMINACLIP:/home/ben/Trimmomatic-0.36/adapters/HiSeqX_overrep_new.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36

java -jar /home/ben/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 -trimlog AO248new.txt AO248_S8_L007_R1_001.fastq.gz AO248_S8_L007_R2_001.fastq.gz AO248_R1_newtrim_paired.fastq.gz AO248_R1_newtrim_single.fastq.gz AO248_R2_newtrim_paired.fastq.gz AO248_R2_newtrim_single.fastq.gz ILLUMINACLIP:/home/ben/Trimmomatic-0.36/adapters/HiSeqX_overrep_new.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:36
```


# FastQC Again
This hopefully now will find no overrepresented adaptor seqs, which is accurate

# Quake
I am working in this directory `/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS`

```
COUNT
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/AO245_R*_trim_*.fastq.gz | jellyfish count /dev/fd/0 -m 19 -s 1000000000 -t 16 -C -o AO245_jelly_count_all_19mers
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO248_R*_trim_*.fastq.gz | jellyfish count /dev/fd/0 -m 19 -s 1000000000 -t 16 -C -o AO248_jelly_count_all_19mers

MERGE (merging actually was not necessary for AO245 or AO248)
jellyfish merge -o AO245_jelly_count_all_19mers.jf AO245_jelly_count_all_19mers\_*
jellyfish merge -o AO248_jelly_count_all_19mers.jf AO248_jelly_count_all_19mers\_*

DUMP
because merging was not necessary, commandline changed from this:
jellyfish dump -c -t AO245_jelly_count_all_19mers.jf -o AO245_jelly_dump_all_19mers
to this:
jellyfish dump -c -t AO245_jelly_count_all_19mers -o AO245_jelly_dump_all_19mers

because merging was not necessary, commandline changed from this:
jellyfish dump -c -t AO248_jelly_count_all_19mers.jf -o AO248_jelly_newdump_all_19mers
to this:
jellyfish dump -c -t AO248_jelly_count_all_19mers -o AO248_jelly_dump_all_19mers

INTERPRET
/usr/local/quake/bin/cov_model.py --int AO245_jelly_newdump_all_19mers
/usr/local/quake/bin/cov_model.py --int AO248_jelly_newdump_all_19mers

CORRECT
/usr/local/quake/bin/correct -f tymp_newuncorrected_data -z -k 19 -c XXX -m AO245_jelly_newdump_all_19mers -p 4
/usr/local/quake/bin/correct -f oct_newuncorrected_data -z -k 19 -c XXX -m AO248_jelly_newdump_all_19mers -p 4

```
# Update

```
COUNT
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/AO245_R*_newtrim_*.fastq.gz | jellyfish count /dev/fd/0 -m 19 -s 1000000000 -t 16 -C -o AO245_jelly_newcount_all_19mers
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS?AO248_R*_newtrim_*.fastq.gz | jellyfish count /dev/fd/0 -m 19 -s 1000000000 -t 16 -C -o AO248_jelly_newcount_all_19mers

MERGE (actually this was not necessary)
jellyfish merge -o AO245_jelly_newcount_all_19mers.jf AO245_jelly_newcount_all_19mers\_*
jellyfish merge -o AO248_jelly_newcount_all_19mers.jf AO248_jelly_newcount_all_19mers\_*

DUMP (because merging was not necessary these changed from this:
jellyfish dump -c -t AO245_jelly_newcount_all_19mers.jf -o AO245_jelly_newdump_all_19mers
jellyfish dump -c -t AO248_jelly_newcount_all_19mers.jf -o AO248_jelly_newdump_all_19mers
to this:
jellyfish dump -c -t AO245_jelly_newcount_all_19mers -o AO245_jelly_newdump_all_19mers
jellyfish dump -c -t AO248_jelly_newcount_all_19mers -o AO248_jelly_newdump_all_19mers



INTERPRET
/usr/local/quake/bin/cov_model.py --int AO245_jelly_newdump_all_19mers
/usr/local/quake/bin/cov_model.py --int AO248_jelly_newdump_all_19mers

CORRECT
/usr/local/quake/bin/correct -f tymp_newuncorrected_data -z -k 19 -c 1 -m AO245_jelly_newdump_all_19mers -p 4
/usr/local/quake/bin/correct -f oct_newuncorrected_data -z -k 19 -c 1 -m AO248_jelly_newdump_all_19mers -p 4

where the -c flag is the cutoff (equal to one in both cases based on the interpret command)

```

# Assembly with Abyss

I made a directory on iqaluk for each genome:
```
/work/ben/2017_Tymp_Oct_HiSeqX/AO248_quaked_data
```
I renamed the paired files as recommended in the abyss manual. Here is the commandline for oct:
```
abyss-pe np=8 name=AO248 lib='pea' k=64 pea='AO248_R1_trim_paired.cor_1.fq.gz AO248_R2_trim_paired.cor_2.fq.gz' se='AO248_R1_trim_paired.cor_single.fastq.gz AO248_R2_trim_paired.cor_single.fastq.gz AO248_R1_trim_single.cor.fastq.gz AO248_R2_trim_single.cor.fastq.gz'
```
and for tympa:
```
/work/ben/2016_XL_ST_liver_RNAseq/Sample_BenEvansBJE4168cDNA_Library
 ```
```
abyss-pe np=8 name=AO248 lib='pea' k=64 pea='AO245_R1_trim_paired.cor_1.fq.gz AO245_R2_trim_paired.cor_2.fq.gz' se='AO245_R1_trim_paired.cor_single.fastq.gz AO245_R2_trim_paired.cor_single.fastq.gz AO245_R1_trim_single.cor.fastq.gz AO245_R2_trim_single.cor.fastq.gz'
```
# Update (on iqaluk)

first unload this module:
`module unload intel mkl openmpi`

then load other modules
```
module load gcc/4.9.2
module load openmpi/gcc492-std/1.8.7
module load boost/gcc492-openmpi187std/1.59.0
export PATH=/work/ben/abyss/bin:$PATH
```

I made a directory on iqaluk for each genome:
```
/work/ben/2017_Tymp_Oct_HiSeqX/AO248_newquaked_data
```
I renamed the paired files as recommended in the abyss manual. 
```
mv AO248_R1_newtrim_paired.cor.fastq.gz AO248_R1_newtrim_paired.cor_1.fq.gz
mv AO248_R2_newtrim_paired.cor.fastq.gz AO248_R2_newtrim_paired.cor_2.fq.gz
```

# Reset path for tmp directory
Abyss crashed due to not enough disk space on /tmp. So reset the TMPDIR variable like this before executing the command. Hopefully this works!

```
export TMPDIR=/scratch/ben/AO245_new_quaked_data
```

Here is the commandline for oct:
```
abyss-pe np=24 name=AO248 lib='pea' k=64 pea='AO248_R1_newtrim_paired.cor_1.fq.gz AO248_R2_newtrim_paired.cor_2.fq.gz' se='AO248_R1_newtrim_paired.cor_single.fastq.gz AO248_R2_newtrim_paired.cor_single.fastq.gz AO248_R1_newtrim_single.cor.fastq.gz AO248_R2_newtrim_single.cor.fastq.gz'
```
and for tympa:
```
/work/ben/2016_XL_ST_liver_RNAseq/Sample_BenEvansBJE4168cDNA_Library
 ```
```
abyss-pe np=24 name=AO248 lib='pea' k=64 pea='AO245_R1_newtrim_paired.cor_1.fq.gz AO245_R2_newtrim_paired.cor_2.fq.gz' se='AO245_R1_newtrim_paired.cor_single.fastq.gz AO245_R2_newtrim_paired.cor_single.fastq.gz AO245_R1_newtrim_single.cor.fastq.gz AO245_R2_newtrim_single.cor.fastq.gz'
```

# Kmer on postquake data
Working in this directory on info.
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS
```
Count with jellyfish
```
19mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/*cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 19 -s 1000000000 -t 16 -C -o AO245_jelly_count_afterquake_19mers -Q 5 
25mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/*cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 25 -s 1000000000 -t 16 -C -o AO245_jelly_count_afterquake_25mers -Q 5 
35mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/*cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 35 -s 1000000000 -t 16 -C -o AO245_jelly_count_afterquake_35mers -Q 5 

19mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/*cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 19 -s 1000000000 -t 16 -C -o AO248_jelly_count_afterquake_19mers -Q 5 
25mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/*cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 25 -s 1000000000 -t 16 -C -o AO248_jelly_count_afterquake_25mers -Q 5
35mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/*cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 35 -s 1000000000 -t 16 -C -o AO248_jelly_count_afterquake_35mers -Q 5

```

Update count with jellyfish
```
19mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/*newtrim*cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 19 -s 1000000000 -t 16 -C -o AO245_newtrim_jelly_count_afterquake_19mers -Q 5 
25mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/*newtrim**cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 25 -s 1000000000 -t 16 -C -o AO245_newtrim_jelly_count_afterquake_25mers -Q 5 
35mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/*newtrim**cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 35 -s 1000000000 -t 16 -C -o AO245_newtrim_jelly_count_afterquake_35mers -Q 5 

19mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/*newtrim**cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 19 -s 1000000000 -t 16 -C -o AO248_newtrim_jelly_count_afterquake_19mers -Q 5 
25mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/*newtrim**cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 25 -s 1000000000 -t 16 -C -o AO248_newtrim_jelly_count_afterquake_25mers -Q 5
35mer
zcat /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/*newtrim**cor*gz | /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish count /dev/fd/0 -m 35 -s 1000000000 -t 16 -C -o AO248_newtrim_jelly_count_afterquake_35mers -Q 5

```


# Update histogram with jellyfish
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish histo -o AO248_newtrim_jelly_count_afterquake_19mers.histo AO248_newtrim_jelly_count_afterquake_19mers
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish histo -o AO248_newtrim_jelly_count_afterquake_25mers.histo AO248_newtrim_jelly_count_afterquake_25mers
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/jellyfish-2.2.6/bin/jellyfish histo -o AO248_newtrim_jelly_count_afterquake_35mers.histo AO248_newtrim_jelly_count_afterquake_35mers

```

# Reparc

It seems like this program is only installed on the evanslab account:

```
/home/evanslab/Hymeno_fastqc/RepArk_analysis/RepARK-master/RepARK.pl
```

I copied this to my account and changed the path of jellyfish in the perl script. Then I issues this command:
```
./RepARK.pl -l AO248_R1_newtrim_paired.cor.fastq.gz -l AO248_R2_newtrim_paired.cor.fastq.gz -l AO248_R2_newtrim_paired.cor_single.fastq.gz -l AO248_R1_newtrim_paired.cor_single.fastq.gz -l AO248_R1_newtrim_single.cor.fastq.gz -l AO248_R2_newtrim_single.cor.fastq.gz -o AO248_newtrim_kmer_31
```

This finished smoothly and results are in this directory:
`/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO248_newtrim_kmer_31`
and expecially here:
`/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO248_newtrim_kmer_31/velvet_repeat_lib`


(Old stuff) The kmer approach using RepArk.pl used a default kmer size of 31 bases. I am going to try a larger kmer size as follows:

```
./RepARK.pl -l AO245_R1_trim_paired.fastq.gz -l AO245_R2_trim_paired.fastq.gz -k 70 -o AO245_kmer_70
./RepARK.pl -l AO248_R1_trim_paired.fastq.gz -l AO248_R2_trim_paired.fastq.gz -k 70 -o AO248_kmer_70
```

The kmer densities can be plotted using this R script:


```
library(ggplot2)

# make a pdf
pdf("kmer_density plot.pdf",w=6, h=2, version="1.4", bg="transparent")

# load the data
tymp<-read.table("jf_RepARK.tymp_histo")
tymp$species <- 'tymp'
oct<-read.table("jf_RepARK.octomys_histo")
oct$species <- 'oct'

# combine the data
data <- rbind(tymp, oct)

# make a density plot (data$V2 has the counts of each kmer)
ggplot(data, aes(V2, fill = species)) +
  # make it transparent and add a limit to the X and Y axes for clarity
  geom_density(alpha = 0.2) + xlim(0,5000) + ylim(0,0.002) +
  # modify the labels
  xlab("Occurance") + ylab("Density") +
  # modify the title
  ggtitle(expression(paste("Densities of 31-mers for ",italic("T. barrerae")," and ",italic("O. mimax")," transcriptomes"))) +
  # get rid of the background
  theme_classic() +
  # fix the legend
  scale_fill_manual(values=c("red", "blue"),
                    name="Species",
                    labels=c(expression(paste(italic("O. mimax"))), expression(paste(italic("T. barrerae")))))+
  # get the legend to be left justified
  theme(legend.text.align   =0) +
  # move the legend over
  theme(legend.position = c(.8, .5))+
  # make the title smaller
  theme(plot.title = element_text(size = 10)) +
  # make legend title smaller too
  theme(legend.title = element_text(size = 8)) +
  # make legend text smaller too
  theme(legend.text = element_text(size = 8)) +
  # make x axis smaller too
  theme(axis.title.x = element_text(size = 8)) +
  # make y axis smaller too
  theme(axis.title.y = element_text(size = 8)) +
  # make axis text smaller too
  theme(axis.text = element_text(size = 8)) +
  geom_segment(aes(x = 500, y = 0.0005, xend = 3000, yend = 0.0005))

dev.off()
# DONE!
```

# Make a blast db out of the high abundance kmer contig

/usr/local/blast/2.3.0/bin/makeblastdb -in /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO245_newtrim_kmer_31/velvet_repeat_lib/contigs.fa -dbtype nucl -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO245_newtrim_kmer_31/velvet_repeat_lib/contigs.fa_blastable


## blast mtdna against octomys high abundance kmer contigs
/usr/local/blast/2.6.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/Tymp_mtDNA_HM544132.1.fa -db /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO245_newtrim_kmer_31/velvet_repeat_lib/contigs.fa_blastable -outfmt 6 -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO245_newtrim_kmer_31/velvet_repeat_lib/Tymp_mtDNA_HM544132.1_to_AO248_highabundancekmercontigs -evalue 1e-20 -task megablast

Two hits for octomys high abundance kmer contigs
```
HM544132.1	NODE_180_length_16302_cov_3.096000	88.508	16289	1794	
62	1	16262	16238	1	0.0	19643
HM544132.1	NODE_16084_length_208_cov_3.980769	91.102	236	21	
0	16534	16769	236	1	2.21e-86	320
```
Total length is 16,566 bp because each of the lengths above (16302 and 208) is 28 bp too short.

## blast octomys high abundance kmer contigs against themselves
/usr/local/blast/2.6.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO245_newtrim_kmer_31/velvet_repeat_lib/contigs.fa -db /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO245_newtrim_kmer_31/velvet_repeat_lib/contigs.fa_blastable -outfmt 6 -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO245_newtrim_kmer_31/velvet_repeat_lib/AO248_highabundancekmercontigs_to_AO248_highabundancekmercontigs -evalue 1e-20 -task megablast

## check how many high abundance kmer contigs have more than one blast hit (`parses_highabund_blasted_to_highabund.pl`):

```
#!/usr/bin/perl
use warnings;
use strict;
use List::MoreUtils qw/ uniq /;

# This program reads in a blast output from a highabund assembly to itself.
# for each portion of each query, I check if that same portion matches more than one place in the genome
# assembly

# here is the format of the blast output

# Column headers:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

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

my $outputfile = "./highabund_doublehits.out";

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile  $!\n\n";
	exit;
}
print "Creating output file: $outputfile\n";

# open blast results

open (DATA, "./AO248_highabundancekmercontigs_to_AO248_highabundancekmercontigs") or die "Failed to open laevis Blast results";

my $genez=0;

my @temp;
my @multimatch;

while ( my $line = <DATA>) {
	@temp = split("\t",$line);
	if($temp[0] ne $temp[1]){
		# this is a match to another 
		push (@multimatch, $temp[0]);
	}
	else{
		$genez+=1;
	}
}
my $x_uniq = uniq @multimatch;

print "The number of contigs that matched at least one other is ",$x_uniq,"\n";
print "The number of genes is ",$genez,"\n";
print "The proportion is ",$x_uniq/$genez,"\n";
```


# Parsing the velvet assembly
I'm in this directory: `/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_kmer_51/velvet_repeat_lib`

To out put the coverage of the top ten contigs assembled by velvet after the Repark.pl analysis, I typed this:
```
grep -o -P '(?<=cov_).*(?=)' contigs.fa | sort -rn | head -n 10
```

# Genome size stuff with jellyfish

### genome size estimation (http://koke.asrc.kanazawa-u.ac.jp/HOWTO/kmer-genomesize.html)
### for tympa (in R)
```
data245 = read.table("jf_RepARK.tymp_histo", header=FALSE)
plot(data245[5:200,],type="h", main="tymp_31mer")
```
### sum the counts above the error threshold (here it is inaccurate because of low coverage; we are only skipping the first one)
sum(as.numeric(as.numeric(data245[2:10000,1])*as.numeric(data245[2:10000,2])))
### find out where the peak is (it is ~12):
data245[15:20,]
### genome size is sum of counts (66,664,660,599) divided by peak count class (12)
sum(as.numeric(as.numeric(data245[2:10000,1])*as.numeric(data245[2:10000,2])))/12
### genome size: 5,555,388,383 (not bad for the 31mer) for the 70mer, the estimated genome size is 8,077,942,850 using a peak value of 6
### and a kmer sum of 48,467,657,099

### size of single copy region
sum(as.numeric(as.numeric(data245[2:50,1])*as.numeric(data245[2:50,2])))/12
### 3166762086
### 3166762086/5555388383
### 0.5700343 of the genome is single copy (for 70mer, ~0.60 is single copy using 2:25 as the range for single copy)


### for octomys
```R
data248 = read.table("jf_RepARK.octomys_histo", header=FALSE)
plot(data248[5:200,],type="h", main="oct_31mer")
sum(as.numeric(as.numeric(data248[2:9995,1])*as.numeric(data248[2:9995,2])))/11
```
#Genome size: 3,283,950,778 (should be ~4gb I think; also not too bad) for 70mer it is 3,417,128,968 (quite similar)

sum(as.numeric(as.numeric(data248[2:45,1])*as.numeric(data248[2:45,2])))/11

(sum(as.numeric(as.numeric(data248[2:45,1])*as.numeric(data248[2:45,2])))/11)/(sum(as.numeric(as.numeric(data248[2:9995,1])*as.numeric(data248[2:9995,2])))/11)

#0.8665599 of the genome is single copy or for 70mer 0.9233744 of the genome is single copy

### now make a plot with poisson expectaton
```R
singleC <- sum(as.numeric(as.numeric(data248[2:45,1])*as.numeric(data248[2:45,2])))/11
pdf("octo_jellyfishkmer_hist_withpoisson.pdf",w=6, h=5, version="1.4", bg="transparent")
plot(1:200,dpois(1:200, 11)*singleC, type = "l", col=3, lty=2)
lines(data248[2:200,],type="l")
dev.off()
```
```R
pdf("tympa_jellyfishkmer_hist_withpoisson.pdf",w=6, h=5, version="1.4", bg="transparent")
plot(1:200,dpois(1:200, 12)*singleC, type = "l", col=3, lty=2)
lines(data245[2:200,],type="l")
dev.off()
```
