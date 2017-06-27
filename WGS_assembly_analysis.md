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
/usr/local/blast/2.3.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/AO248_newtrim_kmer_31/velvet_repeat_lib/contigs.fa -db /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/AO248_newtrim_scaffolds.fa_blastable -outfmt 6 -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/Oct_highabundancekmercontig_blasted_to_Oct_WGS_newtrim_genome_assembly.out
```
```
/usr/local/blast/2.3.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/AO245_newtrim_kmer_31/velvet_repeat_lib/contigs.fa -db /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/AO245_newtrim_scaffolds.fa_blastable -outfmt 6 -out /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/Tymp_highabundancekmercontig_blasted_to_Tymp_WGS_newtrim_genome_assembly.out
```

# Counting high abundance kmer contigs in genome assembly
This script (code on another page) counts how many times unique high abundance contigs are found in each scaffold (parses_highabundancecontigs_blasted_to_WGSgenome_assembly.pl)

it is located here:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/repArc_kmer_29/velvet_repeat_lib/parses_highabundancecontigs_blasted_to_WGSgenome_assembly.pl
```

Results for oct are here:

```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO248_WGS/abyss_genome_assembly/oct_newtrim_highcontigs_to_oct_newtrim_genome.out
```

Results for tymp are here:
```
/net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/tymp_newtrim_highcontigs_to_tymp_newtrim_genome.out
```
According to this link: https://edwards.sdsu.edu/research/perl-one-liner-to-extract-sequences-by-their-identifer-from-a-fasta-file/
Get fasta entries like this:
```
perl -ne 'if(/^>(\S+)/){$c=grep{/^$1$/}qw(42732422)}print if $c' /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/AO245_newtrim_scaffolds.fa
```
or with a file:
```
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' ids.file fasta.file
```
```
awk -v seq="42732422 10004 289064 42567259-,...,42597465-" -v RS='>' '$1 == seq {print RS $0}' /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/AO245_newtrim_scaffolds.fa
```
