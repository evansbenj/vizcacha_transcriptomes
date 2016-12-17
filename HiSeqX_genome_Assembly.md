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
Assemble genome of each species using only paired end reads. 
```
abyss-pe np=8 name=AO245 k=64 in='AO245_R1_trim_paired.fastq.gz AO245_R2_trim_paired.fastq.gz' 
abyss-pe np=8 name=AO248 k=64 in='AO248_R1_trim_paired.fastq.gz AO248_R2_trim_paired.fastq.gz' 
```
