#!/usr/bin/perl

use strict;
use List::Util qw(min max);

my @sump;
my $maxp = 0;

print STDERR "\n";
print STDERR "Sum per-base inputs.\n";

foreach my $f (@ARGV) {
    print STDERR "  $f\n";

    open(F, "< $f") or die "failed to open $f: $!\n";
    while (<F>) {
        my ($pos, $val) = split '\s+', $_;

        $sump[$pos] += $val;
        $maxp = max($maxp, $pos);
    }
}

print STDERR "\n";
print STDERR "Report.\n";

for (my $pos=0; $pos<=$maxp; $pos++) {
    if ($sump[$pos] > 0) {
        print "$pos\t$sump[$pos]\n";
    }
}

