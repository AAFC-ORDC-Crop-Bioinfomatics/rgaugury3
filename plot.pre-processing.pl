#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

# this script is to format input as standard input for gene motif and domain plot

# NBS.merged RLPorRLP.merged need to be standardize


my $version = 0.1;
# -
GetOptions(my $options = {},
              "-l1=s","-l2=s","-l3=s","-l4=s","-nd=s","-rd=s","-o=s"
);

my $USAGE = <<USAGE;
-------------------------------------------------------------------------------------------------------------
gene structure plot pre processing..
 
-------------------------------------------------------------------------------------------------------------
Version: $version, coded by Pingchuan Li @ Frank You Lab

Arguments:

        -l1     nbs valid list
        -l2     rlp valid list
        -l3     rlk valid list
        -l4     tmcc valid list
        
        -nd     nbs merged motif
        -rd     rlp or rlk merged motif 
        
        -o      output file name
        
    enjoy it!
USAGE

die $USAGE unless (defined $options->{nd});

my %hash = (
    'NBS' => 1,
    'LRR' => 2,
    'TIR' => 3,
    'CC'  => 4,
    'stk' => 5,
    'tm'  => 6,
    'sp'  => 7,
    'LysM'=> 8
);

my $output =($options->{o}) ? $options->{o} : "dm.compile.txt" ;

my $nbs_domain_file = $options->{nd};
my $rd_domain_file  = $options->{rd};

my $cvt_l1_ref = convertdata($options->{l1}, $options->{nd}, 2);
my $cvt_l2_ref = convertdata($options->{l2}, $options->{rd}, 1);
my $cvt_l3_ref = convertdata($options->{l3}, $options->{rd}, 1);
my $cvt_l4_ref = convertdata($options->{l4}, $options->{rd}, 1);

open(OUT,">$output");
foreach my $ref ($cvt_l1_ref, $cvt_l2_ref, $cvt_l3_ref, $cvt_l4_ref) {
    my %hash = %{$ref};
    foreach my $id (keys %hash) {
        my @values = @{$hash{$id}};
        print OUT "$id";
        foreach my $value (@values) {
            print OUT "\t$value";
        }
        print OUT "\n";
    }
}
close OUT;


# ------------------------
sub convertdata {
    my ($list,$merged_domain,$init_start) = @_;
    
    my %valid  = ();
    my %result = ();
    
    open(IN, $list) or die "unable to open $list\n";
    while (<IN>) {
        chomp;
        my ($id, $type) = split/\t/,$_;
        $valid{$id}    = 1;
    }
    close IN;
    
    my $i = 0;
    my %number_to_id = ();
    open(IN,$merged_domain) or die "unable to open file $merged_domain";
    while (<IN>) {
        chomp;
        $i++;
        
        my @array = split/\t/,$_;
        
        my @domain;
        if ($i == 1) {
            for my $value ($init_start..$#array) {
                $number_to_id{$value} = $array[$value];
                #print "$value\t $array[$value]\n";
            }
        }
        else {
            next unless ($valid{$array[0]});
            for my $value ($init_start..$#array) {
                my $type = $number_to_id{$value};
                #print "$type\n";
                next if ($array[$value] eq '.');
                my @rangearray = split/\s/,$array[$value];
                
                foreach my $member (@rangearray) {
                    my ($start,$end) = $member =~ /\|(.*?)-(.*)/;
                    ($start,$end) = ($end,$start) if ($start>$end);
                    my $img_value = $hash{$type};
                    push(@{$result{$array[0]}},"motif_$img_value|$start-$end");
                }
            }
        }
    }
    close IN;
    
    return (\%result);    
}







































