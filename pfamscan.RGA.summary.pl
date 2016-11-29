#!/usr/bin/perl -w
use strict;
use Getopt::Long;


my $USAGE = <<USAGE;
Scripts: pfam_scan RGA summary script

arguments:

        -i           inputfile, output by pfam_scan.pl script
        -pfx         prefix

USAGE

#---------- parameters which should be defined in advance -----------
GetOptions(my $options = {},
                    "-i=s","-pfx=s","-s=s"
           );

die "unable to open iputfile\n$USAGE" unless ($options->{i});

my $input = $options->{i};
my $prefix = ($options->{pfx})? $options->{pfx} : "default." ;

#my $summary = $prefix."summary.txt";
#my $summary = $prefix.$options->{s};

open(IN,"$input");
my %summary = ();
my %domain = ();
while (<IN>) {
    chomp;
    my @array = split/\s+/,$_;
    next if (/^#/);
    next if (/^\s/);
    next if (/^>/);
    next unless (/^\w/);

    my $start = $array[1];
    my $end   = $array[2];

    my ($id) = $array[0] =~ /(\S+)/;
    if ($array[6] =~ /LRR/i) {
        $summary{$id}->{LRR}++;
        push(@{$domain{LRR}->{$id}},join("|","domain_LRR","$start-$end"));
    }

    if ($array[6] =~ /NB-ARC/i) {
        $summary{$id}->{NBS}++;
        push(@{$domain{NBS}->{$id}},join("|","domain_NBS","$start-$end"));
    }

    if ($array[6] =~ /TIR/) {
        $summary{$id}->{TIR}++;
        push(@{$domain{TIR}->{$id}},join("|","domain_TIR","$start-$end"));
    }

    # -----------more motif dissect------------------
    #if ($array[6] =~ /PPR/i) {
    #    $summary{$id}->{PPR}++;
    #    push(@{$domain{PPR}->{$id}},join("|","domain_PPR","$start-$end"));
    #}
}
close IN;


my %amount = ();
foreach my $id (sort {$a cmp $b} keys %summary) {
    foreach my $type (keys %{$summary{$id}}) {
        $amount{$type}++;
    }
}

foreach my $domain (keys %domain) {
    my $DOMAIN = $prefix."$domain.res.pfam.txt";
    open(DOMAIN,">$DOMAIN");
    foreach my $id (keys %{$domain{$domain}}) {
        print DOMAIN "$id";
        foreach my $rec (@{$domain{$domain}->{$id}}) {
             print DOMAIN "\t$rec";
        }
        print DOMAIN "\n";
    }
    close DOMAIN;
}
