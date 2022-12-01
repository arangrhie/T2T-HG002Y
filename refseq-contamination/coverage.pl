#!/usr/bin/perl

use strict;
use List::Util qw(min max);

my $winSiz = 10000.0;

my @allSum;

my @winMin;
my @winSum;
my @winMax;
my @winZer;

my $maxp = 0;
my $maxw = 0;

#  Sum all the inputs.

print STDERR "\n";
print STDERR "Sum query-per-base inputs.\n";

foreach my $f (@ARGV) {
    print STDERR "  $f\n";

    open(F, "< $f") or die "failed to open $f: $!\n";
    while (<F>) {
        my ($pos, $val) = split '\s+', $_;

        $allSum[$pos] += $val;
        $maxp = max($maxp, $pos);
    }
}

#  Find min/ave/max for windows.

print STDERR "\n";
print STDERR "Find min/ave/max for $winSiz bp windows.\n";

for (my $pos=0; $pos<=$maxp; $pos++) {
    my $val = $allSum[$pos];

    my $wp  = int($pos / $winSiz + 0.5);   #  int() truncates

    $winMin[$wp]  = $val   if (!defined($winMin[$wp]) && ($val > 0));
    $winSum[$wp]  = 0      if (!defined($winSum[$wp]));
    $winMax[$wp]  = $val   if (!defined($winMax[$wp]));
    $winZer[$wp]  = 0      if (!defined($winZer[$wp]));

    $winMin[$wp]  = min($winMin[$wp], $val)   if ($val > 0);;
    $winSum[$wp] += $val;
    $winMax[$wp]  = max($winMax[$wp], $val);
    $winZer[$wp] += 1   if ($val == 0);

    $maxw = max($maxw, $wp);
}


#open(F, "> mer-per-position") or die "failed to open 'mer-per-position': $!\n";
print STDERR "\n";
print STDERR "Report.\n";

print "#window-center  average-queries-per-base   num-zero  min-per-base  max-per-base\n";
for (my $wp=0; $wp<=$maxw; $wp++) {
    my $l = $winSiz;

    #  The first window covers bases 0 .. winSiz/2
    if ($wp == 0) {
        $l = $winSiz / 2;
    }

    #  The last window covers (wp-1)*winSiz .. maxp
    if ($wp == $maxw) {
        $l = $maxp - ($wp * $winSiz - $winSiz/2)
    }

    if ($winSum[$wp] > 0) {
        printf "%d\t%.2f\t%d\t%d\t%d\n", $wp * $winSiz + $winSiz/2, $winSum[$wp] / $l, $winZer[$wp], $winMin[$wp], $winMax[$wp];
    }
}
#close(F);
