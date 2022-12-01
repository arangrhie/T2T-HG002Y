#!/usr/bin/env perl
use strict;

my %names;

while (<STDIN>) {
    chomp;
    $names{$_} = 0;
}

print STDERR "Loaded ", scalar(keys %names), " names.\n";

open(F, "< ../release213/zzzNames") or die;
while (<F>) {
    chomp;

    my $ln = $_;
    my $id;
    my $ds;

    if ($ln =~ m/^>(.*?)\s(.*)$/) {
        $id = $1;
        $ds = $2;
    } else {
        die "Nope: '$ln'\n";
    }

    if (exists($names{$id})) {
        $names{$id}++;
        print "$id -> $ds\n";
    #} else {
    #    print "$id\n";
    }
}
close(F);


foreach my $id (keys(%names)) {
    if ($names{$id} != 1) {
        print "$id found $names{$id} times.\n";
    }
}
