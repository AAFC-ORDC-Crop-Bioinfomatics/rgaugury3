#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;
use File::Path qw(make_path remove_tree);
use FindBin;
use Log::Log4perl::Level;
use Log::Log4perl qw(:easy);
use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0);

# included in RGAugury package
use Tool qw(blastp_parallel pfamscan_parallel);

#----------------------------------RGAugury pipeline---------------------------------------

my $USAGE = <<USAGE;
Scripts: Resistance Gene Analogs (RGAs) prediction pipeline

 Programmed by Pingchuan Li @ AAFC - Frank You Lab

Usage :perl RGAugury.pl <options>

arguments: 

        -p           protein fasta file
        -n           corresponding cDNA/CDS nucleotide for -p   (optional)
        -g           genome file in fasta format   (optional)
        -d           database, by default value chooses only pfam and gene3d.
        -gff         gff3 file   (optional)
        -c           cpu or threads number, default = 2
        -pfx         prefix for filename, useful for multiple speices input in same folder   (optional)
        -e           E-value for alignment analysis
        
USAGE



#---------- parameters which should be defined in advance -----------
GetOptions(my $options = {},
                    "-p=s","-n=s","-g=s","-gff=s","-c=i","-pfx=s","-d=s","-e=i"
                  );

my $start_run = time();
die $USAGE unless defined ($options->{p});

my $aa_infile   = $options->{p};  
my $nt_infile   = $options->{n};  
my $g_infile    = $options->{g};
my $gff         = $options->{gff};
my $cpu         = ($options->{c})? $options->{c} : 2 ;
my ($prefix)    = ($options->{pfx}) ? "$options->{pfx}" : $aa_infile =~ /([a-zA-Z0-9]+)/; $prefix .= ".";
my $iprDB       = ($options->{d}) ? $options->{d} : "Pfam,gene3d";


# --------------------- cutoff, key interproScan ID and blastp evalue for pfamScan ---------------------------
my $e_seq    = 0.1;
my $e_dom    = 0.1;
my @interested_NBS = keys %{nbarc_parser("$FindBin::Bin/configuration_nbs.txt")};
my $blast_evalue   = (defined $options->{e} and $options->{e} <= 10) ? $options->{e} : "1e-5" ;

# make sure below folder contain pfam and preselected RGA database
my $pfam_index_folder = (-e $ENV{"HOME"}."/database/pfam") ? $ENV{"HOME"}."/database/pfam": die "unable to locate pfam DB";
my $RGD_index_file    = (-e $FindBin::Bin."/RGADB/plant.RGA.dataset.unique.fasta") ? $FindBin::Bin."/RGADB/plant.RGA.dataset.unique.fasta" : die "unalbe to locate RGADB file";

# --------set the directory of coils, be sure ncoils is under the path of RGAugury main directory ----------
$ENV{COILSDIR} = $FindBin::Bin."/coils";

# -------------------  main body -----------------------------
my %NBS_pfam_lst               = ();
my %RGA_blast_lst              = ();
my %protein_fasta              = ();
my %nt_fasta                   = ();
my %genome_fasta               = ();
my %overlap_RGAblast_pfam_lst  = ();
my %NBS_candidates_lst         = ();
my %coils                      = ();
my @deletion                   = ();

#otput file name
my $aa_formated_infile            = $prefix."formated.protein.input.fas";
my $pfam_out                      = $prefix."pfam.local.search.out";
#my $NBS_pfam_out                  = $prefix."NBS.pfam.out";
my $NBS_pre_candidates_lst        = $prefix."NBS.pre.candidates.lst";
my $NBS_merged_domain             = $prefix."NBS.merged.domains.txt";

my $RGA_blast_out                 = $prefix."RGA.blastp.$blast_evalue.out";
my $RGA_blast_fasta               = $prefix."preRGA.candidates.by.Blast.fasta";
my $RGA_blast_lst                 = $prefix."preRGA.candidates.by.Blast.lst";
my $RGA_candidates_fasta          = $prefix."RGA.candidates.fasta";
my $RGA_candidates_lst            = $prefix."RGA.candidates.lst";
my $RGA_candidates_fasta_nt       = $prefix."RGA.candidates.cdna.fasta";
my $candidate_RGA_pfam_out        = $prefix."candidates_RGA_pfam_out";
my $iprscan_out                   = $prefix."iprscan_out.tsv";
my $iprscan_out_2nd               = $prefix."iprscan_out_further.tsv";

my $nbs_prediction                = $prefix."NBS.res.pfam.txt";   #NBS.res.pfam.txt -lrr LRR.res.pfam.txt -tir TIR.res.pfam.txt 
my $lrr_prediction                = $prefix."LRR.res.pfam.txt";
my $tir_prediction                = $prefix."TIR.res.pfam.txt";

my $cc_prediction                 = $prefix."coils.res.txt";
#my $candidate_RGA_lst             = $prefix."candidates_RGA_lst";
my $RLKorRLP_prediction_output    = $prefix."RLKorRLP.domain.prediction.txt";
my $RLKorRLP_merged_domain        = $prefix."RLKorRLP.merged.domains.txt";

my $NBS_candidates_lst            = $prefix."NBS.candidates.lst";
my $NBS_candidates_fas            = $prefix."NBS.candidates.fas";
my $RLK_candidates_lst            = $prefix."RLK.candidates.lst";
my $RLK_candidates_fas            = $prefix."RLK.candidates.fas";  
my $RLP_candidates_lst            = $prefix."RLP.candidates.lst";
my $RLP_candidates_fas            = $prefix."RLP.candidates.fas";
my $TMCC_candidates_lst           = $prefix."TMCC.candidates.lst";
my $TMCC_candidates_fas           = $prefix."TMCC.candidates.fas";

# log file
my $interproscan_log              = $prefix."interproscan.log";
my $status_log                    = $prefix."status.log";

#tmp files
my $tmp_nbsonly_fas               = $prefix."tmp.NBSonly.fas";
my $tmp_nbsonly_lst               = $prefix."tmp.NBSonly.lst";
my $tmp_nbs_ipr                   = $prefix."tmp_nbs_ipr.dissect.txt";
my $tmp_lrr_ipr                   = $prefix."tmp_lrr_ipr.dissect.txt";
my $tmp_tir_ipr                   = $prefix."tmp_tir_ipr.dissect.txt";
my $tmp_cc_ipr                    = $prefix."tmp_cc_ipr.dissect.txt ";  
my $error_report                  = $prefix."Error.logfile.txt";

# figure plot
my $valid_dm                      = $prefix."dm.compile.txt";
my $summary_table                 = $prefix."RGA.summaries.txt";

# --initializing the log4perl modules --------------------------------------------
Log::Log4perl->init("$FindBin::Bin/log4perl.conf");

# ----------------preprocessing protein/DNA fasta sequence-----------------------
DEBUG("step 1 -> Pipeline to predict the plant RGA...");
DEBUG("step 2 -> Make sure all other programs are ready...");
open(ERRREPORT, ">$error_report");

local $/ = ">";

DEBUG("step 3 -> formatting the input file as standard input file...");
%protein_fasta = %{format_fasta($aa_infile, $aa_formated_infile)};

if ($nt_infile and -s $nt_infile) {
    open(IN,$nt_infile) or die "unable to open $nt_infile\n";
    while (<IN>) {
        chomp;
        my ($title,$seq) = split/\n/,$_,2;
        next unless ($title and $seq);
        my ($id) = $title =~ /([a-zA-Z0-9\.\-\_]+)/;
        $nt_fasta{$id} = $seq;
    }
    close IN;
}

if ($g_infile and -s $g_infile) {
    open(IN, $g_infile) or die "unable to open $g_infile\n";
    while (<IN>) {
        chomp;
        my ($title,$seq) = split/\n/,$_,2;
        next unless ($title and $seq);
        my ($id) = $title =~ /([a-zA-Z0-9\.\-\_]+)/;
        $seq   =~ s/\s+//g;
        $genome_fasta{$id}= $seq;
    }
    close IN;
}
local $/ = "\n";

# -------------------------blastp to RGA database------------------
# -----------this step will get a potential candidates RGA---------
# ---all furhter analysis will be based on the preRGA fasta seq----

if ($RGA_blast_out and -s $RGA_blast_out) {
    DEBUG("step 4 -> $RGA_blast_out detected in current folder, pipeline will jumps to next step - code 003");
}
else {
    DEBUG("step 4 -> Blast with RGA DB...");
    blastp_parallel($aa_formated_infile, $RGD_index_file, $blast_evalue, $RGA_blast_out, $cpu);
}

DEBUG("step 5 -> Blast is done, now parsing the output...");

open(IN,$RGA_blast_out);
while (<IN>) {
    chomp;
    my ($geneid,@array) = split/\t/,$_;
    $RGA_blast_lst{$geneid} = 1;  #unique geneid for fasta export purpose
}
close IN;

output_protein_fasta_ram_manner(\%RGA_blast_lst, \%protein_fasta, $RGA_blast_fasta, __LINE__);     #output pre-candidates of RGA, which include NBS-encoding, RLP, RLK and other potential disease resistance genes analogs.-----------
output_lst_ram_manner(\%RGA_blast_lst, $RGA_blast_lst);

#-----------------using pfam to scan the pfam-------------------
#--------remove the $pfam_out if not correctly executed --------


if ($pfam_out and -s $pfam_out) {
    DEBUG("step 6 -> pfam_out detected in current folder, pipeline will jumps to next step - code 001");
}
else {
    DEBUG("step 6 -> running pfam scan procedure...");
    pfamscan_parallel($pfam_out, $e_seq, $e_dom, $cpu, $RGA_blast_fasta, $pfam_index_folder);
}

open(IN, $pfam_out) or die "cant open $pfam_out file\n";  #use @interested_NBS to screen pfam out and acquire all the NBS-encoding genes list
while (<IN>) {
    chomp;
    next if ($_ =~ /^#/ or $_ =~ /^\s/);
    
    my ($geneid, @array) = split/\s+/, $_;
    foreach my $PF (@interested_NBS) {
        if ($_ =~ /$PF/) {
            $NBS_pfam_lst{$geneid} = 1;
        }
    }
}
close IN;

# --------------------coiled coil prediction -------------------

if ($cc_prediction and -s $cc_prediction) {
    DEBUG("step 7 -> $cc_prediction detected in current folder, pipeline will jumps to next step - code 002");
}
else {
    DEBUG("step 7 -> to predict coiled coils...");    
    system("perl -S coils.identification.pl $RGA_blast_fasta $cpu $cc_prediction");
}

open(IN, $cc_prediction);  #$cc_prediction is output of coils.identification.pl
while (<IN>) {
   chomp;
   my ($id,@res) = split/\t/,$_;
   $coils{$id} = join(" ",@res);
}
close IN;

# -------------------------iprscan---------------------------------
#----using above outputed fasta file to do iprscan for 1st time----

if ($iprscan_out and -s $iprscan_out) {
    DEBUG("step 8 -> $iprscan_out detected in current folder, pipeline will jumps to next step - code 004");
}
else {
    DEBUG("step 8 -> initializing interproscan...");    
    system("perl -S iprscan.pl -i $RGA_blast_fasta -appl $iprDB -f tsv -o $iprscan_out -log $interproscan_log");
}


# --------------------------RLK and RLP prediction--------------

if ($RLKorRLP_prediction_output and -s $RLKorRLP_prediction_output) {#$RLKRLP_out_raw
    DEBUG("step 9 -> $RLKorRLP_prediction_output detected in current folder, pipeline will jumps to next step - code 005");
}
else {
    DEBUG("step 9 -> analysising RLK and RLP");
    system("perl -S RLK.prediction.pl  -i  $RGA_blast_fasta  -pfx  $prefix  -pfam  $pfam_out  -iprs $iprscan_out  -cpu $cpu  -lst $RGA_blast_lst   -o $RLKorRLP_prediction_output");
}

#----------- merge coilsed coil to above output<$RLKorRLP_prediciton_outputs ---------------
open(IN,    $RLKorRLP_prediction_output);
open(MERGE, ">$RLKorRLP_merged_domain");
print MERGE join("\t", "id", "stk", "tm", "sp", "LysM", "LRR", "CC\n");


#---------------------- this will add one more column in terms of cc to the $RLKorRLP_prediction_output ---------------------
while (<IN>) {
   chomp;
   next if($_ =~ /LysM\tLRR/); #ignore first headline prior to processing
   my ($id, @content) = split/\t/,$_;
   my $cc = ($coils{$id}) ? $coils{$id} : '.' ;
   print MERGE join("\t", $id, @content, "$cc\n");
}
close IN;
close MERGE;


# ----------dissect NBS.pfam.out------------generate NBS.res.pfam.out etc.------------
#output  NBS.res.pfam.txt LRR.res.pfam.txt, TIR.res.pfam.txt and PPR.res.pfam.txt totall 4 files
system("perl -S pfamscan.RGA.summary.pl        -i   $pfam_out      -pfx $prefix");              

# -extra leucine rich repeat analysis --------------
extra_LRR_analysis("$lrr_prediction");

#merge essential motif/domain to one file
system("perl -S nbs.domain.result.merge.pl     -nbs $nbs_prediction    -lrr $lrr_prediction     -tir $tir_prediction -cc $cc_prediction -seq $aa_formated_infile >$NBS_merged_domain");  

# summary
system("perl -S NBS-encoding.amount.summary.pl -i   $NBS_merged_domain -o   $NBS_pre_candidates_lst -pfx $prefix");

# -------------------------------------- ATTENTION --------------------------------------
# $NBS_merged_domain has some of the false positive NBS protein, because pfam_scan's evalue is 1e-5,
# thus NBS confering only protein sequence need further analysis by doub check of interproscan


# -------------- nbs further analysis --------------
# this session will furhter analysis those nbs-domain only proteins
# -------------------------iprscan further and NBS refine------------------------------
#-analysis 'NBS' type only protein with extra database

open(IN,$NBS_pre_candidates_lst);
open(TMP_NBS_LST,">$tmp_nbsonly_lst");
while (<IN>) {
    chomp;
    my ($id, $type) = split/\t/,$_;
    if ($type eq 'NBS') {
        # for those NBS type, it will be undertaken for further analysis.
        print TMP_NBS_LST "$id\n";#$protein_fasta{$id}\n
    }
    else {
        # for those not NBS tupe, they will be hashed as final NBS-encoding candidates
        $NBS_candidates_lst{$id} = $type;
    }
}
close IN;
close TMP_NBS_LST;

push(@deletion,$tmp_nbsonly_lst) if (-e $tmp_nbsonly_lst);
push(@deletion,$NBS_pre_candidates_lst);


if ($iprscan_out_2nd and -s $iprscan_out_2nd) {
    DEBUG("step 10 -> $iprscan_out_2nd detected in current folder, pipeline will jumps to next step - code 006");
}
else {
    # extract iprscan data from previous analyzed iprscan_out as iprscan_out_2nd
    DEBUG("step 10 -> extracting interproscan for 2nd round of additional NBS encoding genes...");
    iprscan_out_extraction($iprscan_out, $tmp_nbsonly_lst, $iprscan_out_2nd);
}

system("perl -S ipr.specific.id.selection.pl -i $iprscan_out_2nd -o_n $tmp_nbs_ipr -o_l $tmp_lrr_ipr -o_t $tmp_tir_ipr") if ($iprscan_out_2nd and -s $iprscan_out_2nd);#keep the order of output

push(@deletion,$tmp_nbs_ipr) if (-e $tmp_nbs_ipr);
push(@deletion,$tmp_lrr_ipr) if (-e $tmp_lrr_ipr);
push(@deletion,$tmp_tir_ipr) if (-e $tmp_tir_ipr);
push(@deletion,$tmp_cc_ipr ) if (-e $tmp_cc_ipr );

# if this candidates are further confirmed to contain NBS domain by interProScan. then them can be defined as NBS containing only genes,
# beause in previous pfam scan analysis, they have been proved to contain zero tir, cc or lrr motif or domain, then only nbs needs furhter confirmation.
if ($tmp_nbs_ipr and -s $tmp_nbs_ipr) {
    open(IN,$tmp_nbs_ipr) or warn "wrong : $!\n";
    while (<IN>) {
        chomp;
        my ($id,$type) = split/\t/,$_;
        $NBS_candidates_lst{$id} = "NBS";
    }
    close IN;
}

# ------------output and sort final NBS encoding RGAs. ----------------
open(OUT,">$NBS_candidates_lst");
foreach my $key (sort {$NBS_candidates_lst{$a} cmp $NBS_candidates_lst{$b}} keys %NBS_candidates_lst) {
   print OUT "$key\t$NBS_candidates_lst{$key}\n";
}
close OUT;

#----output lst of RLK, RLP and TMCC---
DEBUG("step 11 -> RLK and RLP prediction is done...");
system("perl -S RLK.prediction.result.parser.v2.pl $RLKorRLP_merged_domain $NBS_candidates_lst $RLK_candidates_lst $RLP_candidates_lst $TMCC_candidates_lst");


# ----------output RGA candidates aa sequence --------------
my %id = ();
open(OUT,">$RGA_candidates_fasta");
open(SUMMARY, ">$summary_table");
open(LST,">$RGA_candidates_lst");

foreach my $lst ($NBS_candidates_lst, $RLK_candidates_lst, $RLP_candidates_lst, $TMCC_candidates_lst) {
    open(IN,"$lst");
    while (<IN>) {
        chomp;
        next if (/^\s/);
        my @array = split/\s/,$_;

        my $id   = $array[0];
        my $type = $array[1];
        my $seq  = $protein_fasta{$id};
        
        
        my $new_id = join("|",$id,$type);
        print OUT ">$new_id\n$seq\n";
        print LST "$id\t$type\n";
        
        my $seqlen = length($seq);
        if ($lst eq $NBS_candidates_lst) {
            $type = join(">","NBS",$type);
        }
        
        # -----------------------------create the RGA candidates basic info. SUMMARY table ---------------------        
        if ($gff) {
            print SUMMARY join("\t",$id, $seqlen, $type, "img/$id.png");
        }
        else {
            print SUMMARY join("\t",$id, $seqlen, $type, ".");
        }
        print SUMMARY "\n";
        
        # ------record outputed RGA --------
        $id{$id} = 1;
    }
    close IN;
}
close OUT;
close LST;
close SUMMARY;

# export fasta sequence via candidate list files ------------------------------------------------------------------------------
output_protein_fasta_lst_manner($NBS_candidates_lst ,$NBS_candidates_fas)  if ($NBS_candidates_lst  and -s $NBS_candidates_lst);
output_protein_fasta_lst_manner($RLK_candidates_lst ,$RLK_candidates_fas)  if ($RLK_candidates_lst  and -s $RLK_candidates_lst);
output_protein_fasta_lst_manner($RLP_candidates_lst ,$RLP_candidates_fas)  if ($RLP_candidates_lst  and -s $RLP_candidates_lst);
output_protein_fasta_lst_manner($TMCC_candidates_lst,$TMCC_candidates_fas) if ($TMCC_candidates_lst and -s $TMCC_candidates_lst);

# ------------ output RGAs nucleotide seq if sepcificed in command line------------
if ($nt_infile) {
    open(OUT,">$RGA_candidates_fasta_nt") or warn "unable to open $nt_infile : $!\n";
    foreach my $id (sort {$a cmp $b} keys %id) {
        if (exists $nt_fasta{$id}) {
            print OUT ">$id\n$nt_fasta{$id}";
        }
        else {
            print ERRREPORT "failed to ouptut nucleotide sequence: $id \n";
        }
    }
    close OUT;
}

#------------------prepare CViT data and domain architecture plot figures-=----------------------
if ($gff and -s $gff) {
    DEBUG("step 12 -> creating domain and motif architecture figure...");
    system("perl -S plot.pre-processing.pl -l1 $NBS_candidates_lst -l2 $RLP_candidates_lst -l3 $RLK_candidates_lst -l4 $TMCC_candidates_lst -nd $NBS_merged_domain -rd $RLKorRLP_merged_domain -o $valid_dm");
    system("perl -S plot.gene.domain.motif.pl -gff3   $gff   -i $valid_dm  -s  n" );
    
    DEBUG("step 13 -> creating input file for CViT plotting script package...");
    system("perl -S cvit.input.generator.pl -l $NBS_candidates_lst  -p $aa_infile -f $gff -t gene -c blue   -t2 NBS  -pfx $prefix");
    system("perl -S cvit.input.generator.pl -l $RLK_candidates_lst  -p $aa_infile -f $gff -t gene -c green  -t2 RLK  -pfx $prefix");
    system("perl -S cvit.input.generator.pl -l $RLP_candidates_lst  -p $aa_infile -f $gff -t gene -c orange -t2 RLP  -pfx $prefix");
    system("perl -S cvit.input.generator.pl -l $TMCC_candidates_lst -p $aa_infile -f $gff -t gene -c black  -t2 TMCC -pfx $prefix");
}
else {
    DEBUG("step 12 -> skipped, due to lack of gff3 file...");
    DEBUG("step 13 -> skipped, due to lack of gff3 file...");
}
#------------------ clean files ----------------------------
push(@deletion, $error_report) if (-z $error_report);  #remove Error.log if its' empty logged.
close ERRREPORT;

foreach my $file (@deletion) {
    unlink "$file" or warn "couldnt delete $file: $!\n";
}

my $end_run = time();
DEBUG("step 14 -> You have successfully finished RGA identification, cheering!");
hhmmss_consumed($end_run - $start_run);

#-----------------------------------------------------------
#--------------------------sub functions--------------------
#-----------------------------------------------------------
sub Ptime{
          my $time = localtime;
          my ($msg)= @_;
          print STDERR "$time: $msg\n";
}

sub output_protein_fasta_lst_manner {
    my ($lst,$out) = @_;
    open(IN,$lst) or die "unable to open $lst\n";
    open(OUT,">$out");
    while (<IN>) {
        chomp;
        my ($id,@others) = split/\t/,$_;
        if (exists $protein_fasta{$id}) {
            print OUT ">$id\n$protein_fasta{$id}\n";
        }
        else {
            #DEBUG("not found $id while outputting fasta");
        }
    }
    close IN;
    close OUT;
}


sub output_protein_fasta_ram_manner {
    my ($lst_ref, $source_ref, $output, $line) = @_;
    my @error ;
    my $errorFlag = 0;
    
    my %lst    = %{$lst_ref};
    my %source = %{$source_ref};
    
    open(OUT,">$output");
    foreach my $id (sort {$a cmp $b} keys %lst) {
        if (exists $source{$id}) {
            print OUT ">$id\n$source{$id}\n";
        }
        else {
            $errorFlag++;
            push(@error, $id);
        }
    }
    close OUT;
    
    if ($errorFlag>0) {
        #DEBUG("unalbe to output @error in $line");
    }
}

sub output_lst_ram_manner {
    my ($hash_ref, $output) = @_;
    
    open(OUTPUT,">$output");
    my %hash  = %{$hash_ref};
    foreach my $key (sort {$a cmp $b} keys %hash) {
        print OUTPUT "$key\n";
    }
    close OUTPUT;
}

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

sub format_fasta {
    my ($input, $output) = @_;
    my %fasta = ();
    
    if ($output and -s $output) {
        open(IN,$output) or die "unable to open $output fasta";
        while (<IN>) {
            chomp;
            my ($title, $seq) = split/\n/,$_,2;
            next unless ($title and $seq);
            
            $seq =~ s/\s+//g;
            $fasta{$title} = $seq;
        }
        close IN;
        #DEBUG("$output existed in current folder, reading it to RAM and going ahead with next step...");
    }
    else {
        open(IN, $input) or die "unalbe to open $input\n";
        while (<IN>) {
            chomp;
            my ($title, $seq) = split/\n/,$_,2;
            next unless ($title and $seq);
            my ($id) = $title =~ /([a-zA-Z0-9\.\-\_]+)/;
            
            $seq   =~ s/\s+//g;
            $seq   =~ s/\*//g;
            
            $fasta{$id}= $seq;
        }
        close IN;
        
        open(OUT,">$output");
        foreach my $id (sort {$a cmp $b} keys %fasta) {
            print OUT ">$id\n$fasta{$id}\n";
        }
        close OUT;
    }
    
    return \%fasta;
}

sub hhmmss_consumed {
  my $hourz = int($_[0]/3600);
  my $leftover = $_[0] % 3600;
  my $minz = int($leftover/60);
  my $secz = int($leftover % 60);
  my $consumed = sprintf ("%02d:%02d:%02d", $hourz,$minz,$secz);
  DEBUG("step 15 -> input <$aa_infile> time taken -  $consumed");
}

sub extra_LRR_analysis{
    my $file = shift;
    my %result = ();
    
    # reanalysis on iprscan out for LRR motif
    open(IN, $iprscan_out) or die "unable to open $iprscan_out";
    while (<IN>) {
        chomp;
        my ($id,$md5,$len,$database,$hitid,$desc,$start,$end,$evalue2,$true,$date,$ipr,$domain_name) = split/\t/,$_;
        $desc        = ($desc)? $desc : 'na' ;
        $ipr         = ($ipr) ? $ipr : 'na' ;
        $domain_name = ($domain_name) ? $domain_name : 'na';
        
        if ($desc =~ /leucine.*rich/i or $domain_name =~ /leucine.*rich/i) {
            push(@{$result{$id}}, join("|","IPR_domain_LRR","$start-$end"));
        }
    }
    close IN;

    #  -----------merge with pfam_scan predicted e.g. $lrr_prediction generated by pfamscan.RGA.summary.pl------
    if (-s $file) {
        open(IN,  $file) or die "unable to open  $file";
        while (<IN>) {
            chomp;
            my ($id, @lrr) = split/\t/,$_;
            push(@{$result{$id}}, @lrr);
        }
        close IN;
    }
    
    # rewrite $file contents by >"
    system("rm $file") if (-e $file);
    open(OUT, ">$file") or die "unable to write to $file";
    foreach my $id (sort {$a cmp $b} keys %result) {
        my @content = @{$result{$id}};
        print OUT join("\t",$id, @content);
        print OUT "\n";
    }
    close OUT;
}

sub iprscan_out_extraction{
    my ($iprscan_out, $lst, $output) = @_;
    
    my %lst = ();
    open(IN,$lst) or die "unable to open $lst\n";
    while (<IN>) {
        chomp;
        my @lst = split/\t/,$_;
        $lst{$lst[0]} = 1;
    }
    close IN;
    
    open(IN, $iprscan_out) or die "unable to open $iprscan_out\n";
    open(OUT,">$output");
    while (<IN>) {
        chomp;
        my @array = split/\t/,$_;
        if (exists $lst{$array[0]}) {
            print OUT "$_\n";
        }
    }
    close OUT;
    close IN;
}

sub getLogFilename{
    my $filename = $status_log;
    return $filename;
} 
