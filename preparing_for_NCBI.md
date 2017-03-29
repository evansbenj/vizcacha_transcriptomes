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

A very small proportion of the assembly contained adaptors sequences, so these need to be masked in order to be able to submit them to NCBI:

bedtools maskfasta -fi AO248-scaffolds_greaterthan_200bp.fa -bed adapter.bed -fo AO248-scaffolds_greaterthan_200bp_masked.fa
