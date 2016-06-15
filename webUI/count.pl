#!/usr/bin/perl -w
use strict;
use Scalar::Util qw(looks_like_number);
#count the RGA type and number from RGA.candidates.fas

my %count = ();
my %gene = ();
my $total = 0;

open(IN,$ARGV[0]); # generally the input.RGA.candidates.fasta
while (<IN>) {
    chomp;
    next unless (/^>/);
    my ($gene,$type) = $_ =~ /(.*)\|(.*)/;
    #print "$gene $type\n";
    $count{$type}++;
    $total++;
}
close IN;
my @array  = keys %count;

my $printed_count = 0;
print join("\t","NBS","CNL","TNL","CN","TN","NL","TX","Others","RLP","RLK","TM-CC\n");
foreach my $type ("NBS","CNL","TNL","CN","TN","NL","TX","OTHER","RLP","RLK","TM-CC") {
    if (looks_like_number($count{$type})) {
        print "$count{$type}\t";
        $printed_count += $count{$type};
    }
    else {
        print "0\t";
    }                 
}
print "\n";
#in case there is other type that is not thoughtfully considered. thus an comparison is set up as ;
my $notication = ($total == $printed_count) ? "\n=>looks good<=\n"  : "\n!!!!!!!!!!!! The total and printed number is unequal !!!!!!!!!!!!\n";
print "$notication\nRGA total = $total;\tPrinted_count = $printed_count\n";

