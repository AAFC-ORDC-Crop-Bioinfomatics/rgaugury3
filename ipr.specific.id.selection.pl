#!/usr/bin/perl -w
use strict;
use Scalar::Util qw(looks_like_number);
use FindBin;
use Getopt::Long;


my $USAGE = <<USAGE;
Scripts: for furhter analysis of iprscan

Version: , Written by Pingchuan Li, Frank You Lab

arguments: 

        -i           protein fasta file
        -o_n         
        -o_l         
        -o_t         

USAGE

#---------- parameters which should be defined in advance -----------
GetOptions(my $options = {},
                    "-i=s","-o_n=s","-o_l=s","-o_t=s"
                  );

my $input   = $options->{i};
my $out_nbs = ($options->{o_n}) ? $options->{o_n} : "NBS.res.ipr.txt";
my $out_lrr = ($options->{o_l}) ? $options->{o_l} : "LRR.res.ipr.txt";
my $out_tir = ($options->{o_t}) ? $options->{o_t} : "TIR.res.ipr.txt" ;

die $USAGE unless ($input);

my %result = ();
my @interested_NBS = keys %{nbarc_parser("$FindBin::Bin/configuration_nbs.txt")};

open(IN,$input) or die "unable to open $input\n";
while (<IN>) {
    chomp;
    my ($id,$md5,$len,$database,$hitid,$desc,$start,$end,$evalue,$true,$date,$ipr,$domain_name) = split/\t/,$_;
    push(@{$result{$id}->{$database}},join("|",$len,$database,$hitid,($desc)?$desc:'na',$start,$end,$evalue,$true,($ipr)?$ipr:'na',($domain_name)?$domain_name:'na'));
}
close IN;

my $i = 0;
#use smart to detect TIR  and LRR. use pfam to detect NBS
my %TIR = ();
my %LRR = ();
my %NBS = ();
my %COIL = ();

open(INTER,">intermediate.txt");
foreach my $id (sort {$a cmp $b} keys %result){
    
    my @array = ();
    my $len = 0;
    my $nbsflag = 0;  
    my $lrrflag = 0;
    my $tirflag = 0;
    my $coilflag = 0;
    foreach my $database (sort {$a cmp $b} keys %{$result{$id}}) {
        foreach my $rec (@{$result{$id}->{$database}}) {
            foreach my $config (@interested_NBS) {
                if ($rec =~/$config/i) {
                    $nbsflag = 1;
                    ($len) = $rec =~ /^(\d+?)\|/;
                    
                    my ($len,$database,$hitid,$desc,$start,$end,$evalue,$true,$ipr,$domain_name) = split/\|/,$rec;
                    push(@{$NBS{$id}},join("|","domain_NBS","$start-$end"));
                }
            }
        }
    }
    
    #--------------- only nbs domain containing protein will be further processed by the subfunction loop ----------------
    if ($nbsflag == 1) {
        print INTER "$id\t$len";
        foreach my $database (sort {$a cmp $b} keys %{$result{$id}}) {
            foreach my $rec (sort {$a cmp $b} @{$result{$id}->{$database}}) {
                print INTER "\t$rec";
            }
        }
        print INTER "\n";
    }
    
    
    # --------------------------TIR domain identification ----------------------------
    eval {
        if ($nbsflag == 1) {
            foreach my $rec (@{$result{$id}->{SMART}}) {
                my ($len,$database,$hitid,$desc,$start,$end,$evalue,$true,$ipr,$domain_name) = split/\|/,$rec;
                if ($domain_name =~ /TIR/i) {
                    $tirflag = 1;
                    push(@{$TIR{$id}},join("|","domain_TIR","$start-$end"));
                }
            }
            if ($tirflag == 0) {
                push(@{$TIR{$id}},".");
            }
        }
    };
    
    # ------------------------------Leucien rich repeat -----------------------------      
    eval {
        if ($nbsflag == 1) {
            foreach my $rec (@{$result{$id}->{SMART}}) {
                my ($len,$database,$hitid,$desc,$start,$end,$evalue,$true,$ipr,$domain_name) = split/\|/,$rec;
                if ($desc =~ /Leucine-rich/i or $desc =~ /leucine rich/i) {
                    $lrrflag = 1;
                    push(@{$LRR{$id}},join("|","domain_LRR","$start-$end"));
                }
            }
            if ($lrrflag == 0) {
                push(@{$LRR{$id}},".");
            }
        }
    };
}
close INTER;

open(LRR,">$out_lrr");
foreach my $id (sort {$a cmp $b} keys %LRR) {
    print LRR "$id\t";
    my $value = join(" ",@{$LRR{$id}});
    print LRR "$value\n";
}
close LRR;

open(TIR,">$out_tir");
foreach my $id (sort {$a cmp $b} keys %TIR) {
    print TIR "$id\t";
    my $value = join(" ",@{$TIR{$id}});
    print TIR "$value\n";
}
close TIR;

open(NBS,">$out_nbs");
foreach my $id (sort {$a cmp $b} keys %NBS) {
    print NBS "$id\t";
    my $value = join(" ",@{$NBS{$id}});
    print NBS "$value\n";
}
close NBS;



#  ---------------------------sub------------------------------
sub nbarc_parser {
    my $file = shift;
    my %config = ();
    
    open(IN,"$file") or die "unable to open $file\n";
    while (<IN>) {
        chomp;
        my ($config,@other) = split/\s+/,$_;
        #next unless ($config =~ /^pf/i);
        $config{uc($config)} = 1;
    }
    close IN;
    return \%config;    
}


