# Make blast db

```
/usr/local/blast/2.3.0/bin/makeblastdb -in AO248_newtrim_scaffolds.fa -dbtype nucl -out AO248_newtrim_scaffolds.fa_blastable
```
```
/usr/local/blast/2.3.0/bin/makeblastdb -in AO245_newtrim_scaffolds.fa -dbtype nucl -out AO245_newtrim_scaffolds.fa_blastable
```

# Blast unique transcripts against genome assemblies

```
/usr/local/blast/2.3.0/bin/blastn -query /home/ben/2014_Tympanoctomys_transcriptomes/Octomys/Octomys_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Octomys_all_transcriptomes_assembled_together_unique.fasta -db /4/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248_newtrim_scaffolds.fa_blastable -outfmt 6 -out /4/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/Oct_RNAseq_blasted_to_Oct_WGS_newtrim_genome_assembly.out
```
```
/usr/local/blast/2.3.0/bin/blastn -query /home/ben/2014_Tympanoctomys_transcriptomes/Tympano/Tympano_joint_trinity_assembly_with_concatenated_reads/trinity_out_dir/Tympa_all_transcriptomes_assembled_together_unique.fasta -db /4/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/AO245_newtrim_scaffolds.fa_blastable -outfmt 6 -out /4/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/Tymp_RNAseq_blasted_to_Tymp_WGS_newtrim_genome_assembly.out
```

# Make tab delimited files

```
~/tabix-0.2.6/tabix -p vcf XXXf.vcf.gz
zcat XXX.vcf.gz | /usr/local/vcftools/src/perl/vcf-to-tab > XXX.vcf.gz.tab
```

# Blast high abundance kmer contigs against genome assemblies

oct newtrim high abundance kmer contigs are here:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO248_newtrim_kmer_31/velvet_repeat_lib/contigs.fa
```

tymp newtrim high abundance kmer contigs are here:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/AO245_newtrim_kmer_31/velvet_repeat_lib/contigs.fa
```

oct newtrim genome assembly is here:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248_newtrim_scaffolds.fa
```
tymp newtrim genome assembly is here:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/AO245_newtrim_scaffolds.fa
```


```
/usr/local/blast/2.3.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO248_newtrim_kmer_31/velvet_repeat_lib/contigs.fa -db /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248_newtrim_scaffolds.fa_blastable -outfmt 6 -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO248_newtrim_kmer_31/velvet_repeat_lib/Oct_highabundancekmercontig_blasted_to_Oct_WGS_newtrim_genome_assembly.out
```
```
/usr/local/blast/2.3.0/bin/blastn -query XXX -db /4/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/AO245_newtrim_scaffolds.fa_blastable -outfmt 6 -out /4/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/Tymp_highabundancekmercontig_blasted_to_Tymp_WGS_newtrim_genome_assembly.out
```
