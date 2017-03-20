# Testing whether mtDNA is in genome or transcriptome assemblies.

I am going to blast the entire Tymp mtDNA into the genome and transcriptome assemblies of Tymp and Oct.

## first the genome assemblies

/usr/local/blast/2.6.0/bin/blastn -query Tymp_mtDNA_HM544132.1.fa -db AO245-scaffolds.fa_blastable -outfmt 6 -out Tymp_mtDNA_HM544132.1_to_AO245_genome_assembly -evalue 1e-20 -task megablast 


/usr/local/blast/2.6.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/Tymp_mtDNA_HM544132.1.fa -db AO248-scaffolds.fa_blastable -outfmt 6 -out Tymp_mtDNA_HM544132.1_to_AO248_genome_assembly -evalue 1e-20 -task megablast 

## Now the transcriptime assemblies

/usr/local/blast/2.6.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/Tymp_mtDNA_HM544132.1.fa -db Tympa_all_transcriptomes_assembled_together.fasta_blastable -outfmt 6 -out Tymp_mtDNA_HM544132.1_to_AO245_transcriptome_assembly -evalue 1e-20 -task megablast 


/usr/local/blast/2.6.0/bin/blastn -query /net/infofile4-inside/volume1/scratch/ben/2016_Tympa_and_Octomys_WGS/AO245_WGS/abyss_genome_assembly/Tymp_mtDNA_HM544132.1.fa -db Octomys_all_transcriptomes_assembled_together.fasta_blastable -outfmt 6 -out Tymp_mtDNA_HM544132.1_to_AO248_transcriptome_assembly -evalue 1e-20 -task megablast



