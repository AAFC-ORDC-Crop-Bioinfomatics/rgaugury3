#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

GetOptions(my $options = {},
              "-seq=s","-nbs=s","-lrr=s","-tir=s","-cc=s"
);


my $USAGE = <<USAGE;
-------------------------------------------------------------------------------------------------------------
Script: to merge the prediction file

Usage:   perl $0 <options> <files>
-------------------------------------------------------------------------------------------------------------
 coded by Pingchuan Li
Arguments:
        -seq       the fasta seq, to acquire seq length
        -nbs       nbs filename
        -lrr       lrr filename
        -tir       tir filename
        -cc        cc filename  

enjoy it!
USAGE

my %domain = ();

die $USAGE unless (defined $options->{nbs});

my $fasta   = $options->{seq};
my $nbsfile = $options->{nbs};
my $lrrfile = $options->{lrr};
my $tirfile = $options->{tir};
my $ccfile  = $options->{cc};

domain_sort($nbsfile, "nbs") if (-e $nbsfile);
domain_sort($lrrfile, "lrr") if (-e $lrrfile);
domain_sort($tirfile, "tir") if (-e $tirfile);
domain_sort($ccfile,  "cc" ) if (-e $ccfile );

my %seqlen = ();
open(FASTA,"$fasta");
local $/ = ">";
while (<FASTA>) {
    chomp;
    my ($title,$seq) = split/\n/,$_,2;
    next unless ($title and $seq);
    my ($id) = $title =~ /(\S+)/;
    $seq =~ s/\s+//g;
    my $len = length($seq);
    $seqlen{$id} = $len;
}
close FASTA;

local $/ = "\n";

print "id\tLen\tNBS\tLRR\tTIR\tCC\n";
foreach my $gene (sort {$a cmp $b} keys %domain) {
    my $nbs = ($domain{$gene}->{nbs} and $domain{$gene}->{nbs} =~ /domain/) ? $domain{$gene}->{nbs} : '.' ;
    my $lrr = ($domain{$gene}->{lrr} and $domain{$gene}->{lrr} =~ /domain/) ? $domain{$gene}->{lrr} : '.' ;
    my $tir = ($domain{$gene}->{tir} and $domain{$gene}->{tir} =~ /domain/) ? $domain{$gene}->{tir} : '.' ;
    my $cc  = ($domain{$gene}->{cc}  and $domain{$gene}->{cc}  =~ /domain/) ? $domain{$gene}->{cc}  : '.' ;

    #make sure all identified genes go through below process as NBS-encoding genes
    next unless ($nbs =~ /domain/i or $tir =~ /domain/i);
    
    print join("\t", $gene,$seqlen{$gene},$nbs,$lrr,$tir,$cc);
    print "\n";
    
}

# -----------subfunctions-----------
sub domain_sort {
    my ($filename,$type) = @_;
    open(IN,$filename);
    while (<IN>) {
        chomp;
        my ($id,@content) = split/\s+/,$_;
        $domain{$id}->{$type} = join(" ",@content);
    }
    close IN;
}