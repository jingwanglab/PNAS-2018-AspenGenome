#! usr/bin/perl

use strict;
use File::Basename;

##############################################################################################################
#####Main aim: this script is used to transfer the counts file from the .vcf file called for all sites using UG of P. trichocarpa  to bed file 


my $original=basename($ARGV[0],"\.count");  
my $dirname=dirname($ARGV[0]);
my $bed=$original.".bed";
my $Output=join "/", ($dirname,$bed); 


open COUNT, "<", $ARGV[0] or die "cannot open the .count file: $!";
open OUT, ">", $Output or die "cannot produce the OUT .bed file: $!";

while(<COUNT>) {
	chomp;
	if ($_=~/^CHROM/){print OUT "#chrom\tstart\tend\tn_allele\tref_allele\n";next;}

	my @line = split(/\t/, $_);
	my $chrom=$line[0];
	my $start=$line[1]-1;
	my $end=$line[1];
	my $n_allele=$line[2];
	my $ref_allele=$line[4];

	print OUT "$chrom\t$start\t$end\t$n_allele\t$ref_allele\n";
}

