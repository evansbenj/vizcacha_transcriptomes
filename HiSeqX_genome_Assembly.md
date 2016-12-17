# Genome assembly for Tymp and Oct using Abyss

(On iqaluk), first unload this module:
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
Command for each species
```
abyss-pe np=8 name=BJE3814 k=64 in='BJE3814_S3_L003_R1_001_trim_paired.cor.fastq.gz BJE3814_S3_L003_R2_001_trim_paired.cor.fastq.gz' 
```
