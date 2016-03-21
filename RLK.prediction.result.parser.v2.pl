#!/usr/bin/perl -w
use strict;
use Scalar::Util qw(looks_like_number);




die "to parse the output of RLKorP.assay.Alll.genes.script" unless ($#ARGV >= 0);

my %rlk  = ();
my %rlp  = ();
my %tmcc = ();
my %lrr_lysm = ();  # count those proportion of lrr and lysm in rlp and rlp;
my %receptor_type = ();
my %stk = ();

# acquire the RLK kinase definition from configuration file


my $merged_domain =  $ARGV[0];              # file like $prefix."RLKorRLP.merged.domains.txt";
my $exception     = ($ARGV[1]) ? $ARGV[1]: "input.NBS.candidates.lst";

my $rlk_output    = ($ARGV[2]) ? $ARGV[2]: "input.RLK.candidates.lst";
my $rlp_output    = ($ARGV[3]) ? $ARGV[3]: "input.RLP.candidates.lst";
my $tmcc_output   = ($ARGV[4]) ? $ARGV[4]: "input.TMCC.candidates.lst";


#  processing some of exception -------------
my %exception = ();
my @overlap = ();
open(IN,$exception) or die "unable to open $exception \n";
while (<IN>) {
    chomp;
    my ($geneid,@other) = split/\t/,$_;
    $exception{$geneid} = 1;
}
close IN;

#  processing merged domain list -----------
open(IN,"$merged_domain") or die "unable to open ARGV[0]\n"; #merged domain
while (<IN>) {
    chomp;
    next if ($_ =~ /LysM\tLRR\tCC/);

    my ($id, $stk, $tm, $sp, $lysm, $lrr, $cc) = split/\t/,$_;

    #in extreme situation, the domain results between NBS encoding and RLP/RLK would have overlapping, thus they can be filtered out by below statement.
    # thus each of RLK and RLP classificaiton needs a futher examination for NBS redundancy.

    if ($stk ne '.') {
        $stk{$id} = 1;
    }

    if ($stk ne '.' and $tm ne '.' and $lrr ne '.' and $lysm ne '.') {
        my $flag = redundant_check($id,1);
        next if ($flag == 1);
        $rlk{$id} = 1;
        
        $lrr_lysm{'rlk'}->{'lrr_lysm'}++;
        $receptor_type{$id} = 'lrr_lysm';
    }
    elsif($stk ne '.' and $tm ne '.' and $lysm ne '.') {
        my $flag = redundant_check($id,1);
        next if ($flag == 1);
        $rlk{$id} = 1;
        
        $lrr_lysm{'rlk'}->{'lysm'}++;
        $receptor_type{$id} = 'lysm';
    }
    elsif($stk ne '.' and $tm ne '.' and $lrr ne '.') {
        my $flag = redundant_check($id,1);
        next if ($flag == 1);
        $rlk{$id} = 1;
        
        $lrr_lysm{'rlk'}->{'lrr'}++;
        $receptor_type{$id} = 'lrr';
    }
    elsif ($stk ne '.' and $tm ne '.') {
        my $flag = redundant_check($id,1);
        next if ($flag == 1);
        $rlk{$id} = 1;
    }
    
    
    elsif($stk eq '.' and $tm ne '.' and  $lrr ne '.' and  $lysm ne '.') {
        my $flag = redundant_check($id,2);
        next if ($flag == 1);
        $rlp{$id} = 1;
        
        $lrr_lysm{'rlp'}->{'lrr_lysm'}++;
        $receptor_type{$id} = 'lrr_lysm';
    }
    elsif($stk eq '.' and $tm ne '.' and  $lysm ne '.') {
        my $flag = redundant_check($id,2);
        next if ($flag == 1);
        $rlp{$id} = 1;
        $lrr_lysm{'rlp'}->{'lysm'}++;
        $receptor_type{$id} = 'lysm';
    }
    elsif($stk eq '.' and $tm ne '.' and  $lrr ne '.') {
        my $flag = redundant_check($id,2);
        next if ($flag == 1);
        $rlp{$id} = 1;
        $lrr_lysm{'rlp'}->{'lrr'}++;
        $receptor_type{$id} = 'lrr';
    }
    

    elsif($stk eq '.' and $sp ne '.' and  $lrr ne '.' and  $lysm ne '.') { # either tm or sp is not equal to '.'
        my $flag = redundant_check($id,2);
        next if ($flag == 1);
        $rlp{$id} = 1;
        
        print STDERR "okokok1\n";
        
        $lrr_lysm{'rlp'}->{'lrr_lysm'}++;
        $receptor_type{$id} = 'lrr_lysm';
    }
    elsif($stk eq '.' and $sp ne '.' and  $lysm ne '.') {
        my $flag = redundant_check($id,2);
        next if ($flag == 1);
        $rlp{$id} = 1;
        
        
        print STDERR "okokok2\n";
        
        $lrr_lysm{'rlp'}->{'lysm'}++;
        $receptor_type{$id} = 'lysm';
    }
    elsif($stk eq '.' and $sp ne '.' and  $lrr ne '.') {
        my $flag = redundant_check($id,2);
        next if ($flag == 1);
        $rlp{$id} = 1;
        
        print STDERR "okokok3\n";
        
        $lrr_lysm{'rlp'}->{'lrr'}++;
        $receptor_type{$id} = 'lrr';
    }
    
    
    elsif($stk eq '.' and $tm ne '.' and $cc ne '.') {
        my $flag = redundant_check($id,3);
        next if ($flag == 1);
        $tmcc{$id} = 1;
    }
    else {
        next;
    }
}
close IN;

open(OUT1,">$rlk_output");
foreach my $id (sort {$a cmp $b} keys %rlk) {
    my $receptor_type = (exists $receptor_type{$id}) ? $receptor_type{$id} : "other_receptor";
    print OUT1 "$id\tRLK\t$receptor_type\n";
}
close OUT1;

open(OUT2,">$rlp_output");
foreach my $id (sort {$a cmp $b} keys %rlp) {
    my $receptor_type = $receptor_type{$id};
    print OUT2 "$id\tRLP\t$receptor_type\n";
}
close OUT2;

open(OUT3,">$tmcc_output");
my @tmcc_ids = keys %tmcc;
if ($#tmcc_ids >= 0) {
    foreach my $id (sort {$a cmp $b} keys %tmcc) {
        print OUT3 "$id\tTM-CC\n";
    }
}
close OUT3;

Ptime("WARNING!!!: @overlap,  might be NBS encoding, RLK, RLP or either TMCC") if ($#overlap>=0);

my $rlp_lrr  = (looks_like_number($lrr_lysm{'rlp'}->{'lrr'}))  ? $lrr_lysm{'rlp'}->{'lrr'}  : 0 ;
my $rlp_lysm = (looks_like_number($lrr_lysm{'rlp'}->{'lysm'})) ? $lrr_lysm{'rlp'}->{'lysm'} : 0 ;

my $rlk_lrr  = (looks_like_number($lrr_lysm{'rlk'}->{'lrr'}))  ? $lrr_lysm{'rlk'}->{'lrr'}  : 0 ;
my $rlk_lysm = (looks_like_number($lrr_lysm{'rlk'}->{'lysm'})) ? $lrr_lysm{'rlk'}->{'lysm'} : 0 ;


#print STDERR "RLP lrr = $rlp_lrr; RLP lysm = $rlp_lysm; RLK lrr = $rlk_lrr; RLK lysm = $rlk_lysm\n";
#my @stk = keys %stk;
#print STDERR "stk # = ",$#stk + 1,"\n";

# --------------------------------------sub ----------------------------------------------
sub Ptime{
          my $time = localtime;
          my ($msg)= @_;
          print STDERR "$time: $msg\n";
}

sub redundant_check{
    my ($id,$tag) = @_;
   
    if (exists $exception{$id}) {
        push(@overlap,$id);
        
        #Ptime("$id,$tag");  # for debug
        return 1;
    }
    else {
        return 0;
    }
}

sub file_parser {
    my $file = shift;
    my %pfam = ();
    
    open(IN,"$file") or die "unable to open $file\n";
    while (<IN>) {
        chomp;
        my ($pf,@other) = split/\s+/,$_;
        next unless ($pf =~ /^pf/i);
        $pfam{$pf} = 1;
    }
    close IN;
    my @pfam = keys %pfam;

    return @pfam;    
}
