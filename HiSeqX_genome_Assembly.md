# Genome assembly for Tymp and Oct using Abyss

(On iqaluk) I am working in this directory:
```
/work/ben/2017_Tymp_Oct_HiSeqX
```
First unload this module:
```
module unload intel mkl openmpi
```
Then load other modules:
```
module load gcc/4.9.2
module load openmpi/gcc-4.9.2/std/1.8.7
module load boost/gcc492-openmpi187/1.59.0
export PATH=/work/ben/abyss/bin:$PATH
```
This worked:
```
module load gcc/4.9.2
module load openmpi/gcc492-std/1.8.7
module load boost/gcc492-openmpi187std/1.59.0
export PATH=/work/ben/abyss/bin:$PATH
```

Assemble genome of each species using only paired end reads. 
```
abyss-pe np=8 name=AO245 k=64 in='AO245_R1_trim_paired.fastq.gz AO245_R2_trim_paired.fastq.gz' 
abyss-pe np=8 name=AO248 k=64 in='AO248_R1_trim_paired.fastq.gz AO248_R2_trim_paired.fastq.gz' 
```

Assembly of Octomys genome is completed.  Here are the stats:
```
The minimum coverage of single-end contigs is 0.472727.
The minimum coverage of merged contigs is 2.83636.
Consider increasing the coverage threshold parameter, c, to 2.83636.
ln -sf AO248-8.fa AO248-scaffolds.fa
PathOverlap --overlap   -k64 --dot AO248-7.dot AO248-7.path >AO248-8.dot
ln -sf AO248-8.dot AO248-scaffolds.dot
abyss-fac  AO248-unitigs.fa AO248-contigs.fa AO248-scaffolds.fa |tee AO248-stats.tab
n       n:500   L50     min     N80     N50     N20     E-size  max     sum    name
6294392 962523  180653  500     1677    3812    7731    5100    68869   2.404e9AO248-unitigs.fa
5585105 819908  139635  500     2035    4872    10309   6740    73917   2.44e9 AO248-contigs.fa
5515621 783585  129154  500     2143    5223    11254   7336    74984   2.448e9AO248-scaffolds.fa
ln -sf AO248-stats.tab AO248-stats
tr '\t' , <AO248-stats.tab >AO248-stats.csv
abyss-tabtomd AO248-stats.tab >AO248-stats.md
```

directory:
`/work/ben/2017_Tymp_Oct_HiSeqX/AO248_quaked_data`
assembly is either `AO248-contigs.fa` or `AO248-unitigs.fa`


# Assembling repeat elements using Repark from wgs reads

From within this directory: `/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS`
```
../RepARK.pl -l AO245_R1_trim_paired.cor.fastq.gz -l AO245_R2_trim_paired.cor.fastq.gz -l AO245_R1_trim_paired.cor_single.fastq.gz -l AO245_R1_trim_single.cor.fastq.gz -l AO245_R2_trim_paired.cor_single.fastq.gz -l AO245_R2_trim_single.cor.fastq.gz -k 31 -o repArc_kmer_31
```
and this directory `/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS`
```
../RepARK.pl -l AO248_R1_trim_paired.cor.fastq.gz -l AO248_R2_trim_paired.cor.fastq.gz -l AO248_R1_trim_paired.cor_single.fastq.gz -l AO248_R2_trim_paired.cor_single.fastq.gz -l AO248_R1_trim_single.cor.fastq.gz -l AO248_R2_trim_single.cor.fastq.gz -k 31 -o repArc_AO248_kmer_31
```
I also did this with a k-mer size of 45. I used this one liner to output the highest coverage contig in the velvet_repeat_lib directory:
```
grep -o -P '(?<=cov_).*(?=)' contigs.fa | sort -rn | head -n 1
```
This gets the headers of the entries with the highest coverage:
```
grep -o -P '(?<=cov_).*(?=)' contigs.fa | sort -rn | head -n 10
```

print out a fasta file
```
awk -v seq="NODE_9_length_29_cov_3.000000" -v RS='>' '$1 == seq {print RS $0}' contigs.fa
```

# Make a blast database out of the assembly
```
/usr/local/blast/2.3.0/bin/makeblastdb -in /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248-scaffolds.fa -dbtype nucl -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248-scaffolds.fa_blastable
```
and to blast the repeats from Octomys against the genome assembly
```
/usr/local/blast/2.3.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa -db /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248-scaffolds.fa_blastable -outfmt 6 -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/repArc_AO248_kmer_31_to_AO248-scaffolds.fa_blastable -evalue 1e-20 -task megablast 
```
In the above command, could add `-max_target_seqs 1` if we want to restrict the output to some number of entries.

# Examining the incidence of repetitive elements in the Assembly

Here is the repeat element assembly:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa
```
Here is the results of the blast output
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/repArc_AO248_kmer_31_to_AO248-scaffolds.fa_blastable
```

So, how many times does the most repetitive element show up in the assembly?  First find the coverage of the 10 most repetive elements like this:
```
grep -o -P '(?<=cov_).*(?=)' /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/repArc_AO248_kmer_31/velvet_repeat_lib/contigs.fa | sort -rn | head -n 10
```
```
82.937500
82.903229
77.638885
74.660004
73.687500
73.628571
73.393936
72.500000
71.366669
70.621620
```

Now, get the header of the highest coverage one.
```
grep '82.937500' /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/repArc_kmer_31/velvet_repeat_lib/contigs.fa
```
```
>NODE_251811_length_32_cov_82.937500
```
print out a fasta file
```
awk -v seq="NODE_251811_length_32_cov_82.937500" -v RS='>' '$1 == seq {print RS $0}' /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/repArc_kmer_31/velvet_repeat_lib/contigs.fa
```
