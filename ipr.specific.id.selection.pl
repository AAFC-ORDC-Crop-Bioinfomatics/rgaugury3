#!/usr/bin/perl -w
use strict;
use Scalar::Util qw(looks_like_number);

die "to extract the id of NBS, LRR, etc. \nperl script.pl iprscan.tsv" unless ($#ARGV >= 0);

my $out_nbs = ($ARGV[1]) ? $ARGV[1] : "NBS.res.ipr.txt";
my $out_lrr = ($ARGV[2]) ? $ARGV[2] : "LRR.res.ipr.txt";
my $out_tir = ($ARGV[3]) ? $ARGV[3] : "TIR.res.ipr.txt" ;
my $out_cc  = ($ARGV[4]) ? $ARGV[4] : "CC.res.ipr.txt";



open(IN,"$ARGV[0]");
my %result = ();
while (<IN>) {
    chomp;
    my ($id,$md5,$len,$database,$hitid,$desc,$start,$end,$evalue,$true,$date,$ipr,$domain_name) = split/\t/,$_;
    #$result{$id}->{$database}->{len} = $len;
    #$result{$id}->{$database}->{hitid} = $hitid;
    #$result{$id}->{$database}->{desc} = ($desc)?$desc:'na';
    #$result{$id}->{$database}->{start} = $start;
    #$result{$id}->{$database}->{end} = $end;
    #$result{$id}->{$database}->{evalue} = $evalue;
    #$result{$id}->{$database}->{true} = $true;
    #$result{$id}->{$database}->{date} = $date;
    #$result{$id}->{$database}->{ipr} = ($ipr)?$ipr:'na';
    #$result{$id}->{$database}->{domain_name} = ($domain_name)?$domain_name:'na';
    
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
            if ($rec =~/PF00931/i) {
                $nbsflag = 1;
                ($len) = $rec =~ /^(\d+?)\|/;
                
                my ($len,$database,$hitid,$desc,$start,$end,$evalue,$true,$ipr,$domain_name) = split/\|/,$rec;
                push(@{$NBS{$id}},join("|","domain_NBS","$start-$end"));
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
                if ($desc =~ /Leucine-rich/i) {
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

#open(COIL,">$out_cc");
#foreach my $id (sort {$a cmp $b} keys %COIL) {
#    print COIL "$id\t";
#    my $value = join(" ",@{$COIL{$id}});
#    print COIL "$value\n";
#}
#close COIL;
