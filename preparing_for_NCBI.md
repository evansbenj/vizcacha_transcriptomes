# Preparing data for NCBI

For genome assemblies, NSBI only will take assembled bits that are greater than 200bp.

I got this script off of the internet (removessmall.pl). I think the command line is something like ./removessmall.pl 200 fasta.in > fasta_filtered.out

```perl
## removesmalls.pl
#!/usr/bin/perl
use strict;
use warnings;

my $minlen = shift or die "Error: `minlen` parameter not provided\n";
{
    local $/=">";
    while(<>) {
        chomp;
        next unless /\w/;
        s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
        my $seqlen = length join "", @chunk;
        print ">$_" if($seqlen >= $minlen);
    }
    local $/="\n";
}
```

# finding the adapter in the assembly
/usr/local/blast/2.6.0/bin/blastn -query ../../AO245_WGS/abyss_genome_assembly/illumina_adaptor.fa -db AO248-scaffolds.fa_blastable -outfmt 6 -out illumina_adaptor_to_AO248WGS_assembly -evalue 1e-5 -max_target_seqs 1000000

/usr/local/blast/2.6.0/bin/blastn -query illumina_adaptor.fa -db AO245-scaffolds.fa_blastable -outfmt 6 -out illumina_adaptor_to_AO245WGS_assembly -evalue 1e-5 -max_target_seqs 1000000


# Masking 
A very small proportion of the assembly contained adaptors sequences, so these need to be masked in order to be able to submit them to NCBI:

bedtools maskfasta -fi AO248-scaffolds_greaterthan_200bp.fa -bed adapter.bed -fo AO248-scaffolds_greaterthan_200bp_masked.fa
