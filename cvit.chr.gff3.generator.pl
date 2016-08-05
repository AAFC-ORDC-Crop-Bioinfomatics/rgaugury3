#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

# generatee the chromosome gff3 for cvit

#----------------------------------RGAugury pipeline---------------------------------------

my $USAGE = <<USAGE;
Scripts: to create a chromosome gff3 based on genomics gff3 input

 Programmed by Pingchuan Li @ AAFC - Frank You Lab

Usage :perl cvit.chr.gff3generator.pl <xxxxx.gff3> <cvit.gff3>

arguments: 

        -i           gff3 file
        -o           output for cvit package
        
USAGE


GetOptions(my $options = {},
                    "-i=s","-o=s"
                  );

my $input = $options->{i};
my $output = $options->{o};
my %chr = ();


open(IN,$input) or die "unable to open $input";
while (<IN>) {
    chomp;
    next if (/^#/ or /^\s*$/);
    my @array = split/\t/,$_;
    
    #
    next unless(/^chr\d+/i);
    my ($chr) = $array[0] =~ /chr(\d+)/i;
    
    if (looks_like_number($chr) and $chr>=1) {
        if ($chr{$chr}) {
            if ($chr{$chr} < max($array[3],$array[4])) {
                $chr{$chr} = max($array[3],$array[4]);
            }
        }
        else {
            $chr{$chr} = max($array[3],$array[4]);
        }
    }
}
close IN;

open(OUT,">$output");

my $maxlen = 0;
foreach my $chr (sort {$a <=> $b} keys %chr) {
    print OUT join("\t","Chr".$chr, "RGAugury", "chromosome", 1, $chr{$chr},".",".",".","ID=Chr".$chr);
    print OUT "\n";
    
    if ($chr{$chr}>=$maxlen) {
        $maxlen = $chr{$chr};
    }
}
close OUT;

print "max length of chr = $maxlen";

sub max {
    my ($s1,$s2) = @_;
    if ($s1>=$s2) {
        return $s1;
    }
    elsif ($s1 < $s2) {
        return $s2;
    }
    else {
        die "wrong code 2111";
    }
}