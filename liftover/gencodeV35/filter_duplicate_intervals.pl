#!/usr/bin/perl -w
#
use strict;

my $Usage = "filter_duplicate_intervals.pl <bed file>\n";

$#ARGV==0
    or die "$Usage";

my $bedfile = $ARGV[0];
#my $metadata = '../../genes_and_metadata/gencodeV35.meta.tsv';
#my $ENSTbed = 'gencodeV35.primary.chm13v2.ucsc.sort.bed';

open BED, $bedfile
    or die "Couldn\'t open $bedfile: $!\n";

my $rh_observed_intervals = {};
while (<BED>) {
    chomp;
    my @fields = split /\s/, $_;
    if (!$rh_observed_intervals->{"$fields[0]:$fields[1]-$fields[2]"}) {
        print "$_\n";
        $rh_observed_intervals->{"$fields[0]:$fields[1]-$fields[2]"} = 1;
    }
}
close BED;

