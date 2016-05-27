#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $USAGE = <<USAGE;
Scripts: pre-NBS encoding summary

arguments: 

        -i           inputfile, merged NBS relevant domain/motif file
        -pfx         prefix
        -o           output file
        
USAGE

#---------- parameters which should be defined in advance -----------
GetOptions(my $options = {},
                    "-i=s","-pfx=s","-o=s"
           );

die "unable to open iputfile\n$USAGE" unless ($options->{i});

my $i = 0;
my %total = ();
my %gene  = ();

my $input = $options->{i};
my $prefix = ($options->{pfx}) ? $options->{pfx}: "default.";
my $out = ($options->{o}) ? $options->{o}: $prefix."NBS.pre.candidates.lst" ;
my $othertype = $prefix."other.type.candidates.txt";

open(IN,"$input");
open(OTHER,">$othertype");
open(TYPE, ">$out");
while (<IN>) {
    chomp;
    $i++;
    next if ($i ==1); # ignore the first line which is usually the title of each column.
    
    my ($id,$len,$nbs,$lrr,$tir,$cc,@other)  = split/\t/,$_;
    my $v_nbs = domain_detect($nbs);
    my $v_lrr = domain_detect($lrr);
    my $v_cc  = domain_detect($cc);
    my $v_tir = domain_detect($tir);
    
    if ($v_nbs == 1                                          and $v_lrr + $v_cc + $v_tir          == 0) {
        $total{NBS}++;
        push(@{$gene{NBS}},$id);
    }
    
    elsif ($v_nbs == 1  and  $v_lrr == 1                     and $v_nbs + $v_lrr + $v_cc + $v_tir == 2) {
        $total{NL}++;
        push(@{$gene{NL}},$id);
    }
    
    elsif ($v_cc == 1  and  $v_nbs == 1  and $v_lrr == 1     and $v_nbs + $v_lrr + $v_cc + $v_tir == 3) {
        $total{CNL}++;
        push(@{$gene{CNL}},$id);
    }

    elsif ($v_tir == 1 and  $v_nbs == 1  and $v_lrr == 1     and $v_nbs + $v_lrr + $v_cc + $v_tir == 3) {
        $total{TNL}++;
        push(@{$gene{TNL}},$id);
    }
    
    elsif ($v_tir == 1 and  $v_nbs == 1                      and $v_nbs + $v_lrr + $v_cc + $v_tir == 2) {
        $total{TN}++;
        push(@{$gene{TN}},$id);
    }
    
    elsif ($v_tir == 1                                       and $v_nbs +          $v_cc + $v_tir == 1) {
        $total{TX}++;
        push(@{$gene{TX}},$id);
    }
    
    elsif ($v_cc == 1  and  $v_nbs == 1                      and $v_nbs + $v_lrr + $v_cc + $v_tir == 2) {
        $total{CN}++;
        push(@{$gene{CN}},$id);
    }

#   other type summary    
    else {
        $total{other}++;
        push(@{$gene{OTHER}},$id);
        
        # for debug purpose
        print OTHER "$_\t";
        print OTHER join("|",$v_nbs,$v_lrr,$v_cc,$v_tir,$v_nbs + $v_lrr + $v_cc + $v_tir,"\n");
    }
}

close IN;

foreach my $type (sort {$a cmp $b} keys %total) {
    #print "$type\t$total{$type}\n";
}

foreach my $type (keys %gene) {
    my @list = @{$gene{$type}};
    foreach my $id (@list) {
        print TYPE "$id\t$type\n";
    }
}

close TYPE;
close OTHER;
# ---------------sub-------------------
sub domain_detect {
    my $value = shift;
    my $res = 0;
    if ($value =~ /domain/i) {
        $res = 1;
    }
    return $res;
}

sub Ptime{
          my $time = localtime;
          my ($msg)= @_;
          print STDERR "$time: $msg\n";
}