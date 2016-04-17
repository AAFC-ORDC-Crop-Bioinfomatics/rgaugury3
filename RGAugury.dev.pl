#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;
use File::Path qw(make_path remove_tree);
use FindBin;

#----------------------------------RGAugury pipeline---------------------------------------

my $USAGE = <<USAGE;
Scripts: Resistance Gene Analogs (RGAs) prediction pipeline

 Programmed by Pingchuan Li @ AAFC - Frank You Lab

Usage :perl RGAugury.pl <options>

arguments: 

        -p           protein fasta file
        -n           corresponding cDNA/CDS nucleotide for -p   (optional)
        -g           genome file in fasta format   (optional)
        -gff         gff3 file   (optional)
        -c           cpu or threads number, default = 2
        -pfx         prefix for filename, useful for multiple speices input in same folder   (optional)
        
USAGE

#---------- parameters which should be defined in advance -----------
GetOptions(my $options = {},
                    "-p=s","-n=s","-g=s","-gff=s","-c=i","-pfx=s"
                  );
my $start_run = time();
die $USAGE unless defined ($options->{p});

my $aa_infile   = $options->{p};  
my $nt_infile   = $options->{n};  
my $g_infile    = $options->{g};
my $gff         = $options->{gff};
my $cpu         = ($options->{c})? $options->{c} : 2 ;
my ($prefix)    = ($options->{pfx}) ? "$options->{pfx}" : $aa_infile =~ /([a-zA-Z0-9]+)/; $prefix .= ".";

# --------------------- cutoff, key interproScan ID and blastp evalue for pfamScan ---------------------------
my $e_seq    = 0.1;
my $e_dom    = 0.1;
my @interested_NBS = keys %{nbarc_parser("$FindBin::Bin/configuration_nbs.txt")};
my $blast_evalue   = "1e-5";

#make sure below folder contain pfam and preselected RGA database
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
my %coils = ();
my @deletion                   = ();

#otput file name
my $aa_formated_infile            = $prefix."formated.protein.input.fas";
my $pfam_out                      = $prefix."pfam.local.search.out";
my $NBS_pfam_out                  = $prefix."NBS.pfam.out";
my $NBS_pre_candidates_lst        = $prefix."NBS.pre.candidates.lst";
my $NBS_candidates_lst            = $prefix."NBS.candidates.lst";
my $NBS_candidates_fas            = $prefix."NBS.candidates.fas";
my $NBS_merged_domain             = $prefix."NBS.merged.domains.txt";

my $RGA_blast_out                 = $prefix."RGA.blastp.$blast_evalue.out";
my $RGA_blast_fasta               = $prefix."preRGA.candidates.by.Blast.fasta";
my $RGA_blast_lst                 = $prefix."preRGA.candidates.by.Blast.lst";
my $RGA_candidates_fasta          = $prefix."RGA.candidates.fasta";
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

my $RLK_candidates_lst            = $prefix."RLK.candidates.lst";
my $RLK_candidates_fas            = $prefix."RLK.candidates.fas";  
my $RLP_candidates_lst            = $prefix."RLP.candidates.lst";
my $RLP_candidates_fas            = $prefix."RLP.candidates.fas";
my $TMCC_candidates_lst           = $prefix."TMCC.candidates.lst";
my $TMCC_candidates_fas           = $prefix."TMCC.candidates.fas";

#tmp files
my $tmp_nbsonly_fas               = $prefix."tmp.NBSonly.fas";
my $tmp_nbs_ipr                   = $prefix."tmp_nbs_ipr.dissect.txt";
my $tmp_lrr_ipr                   = $prefix."tmp_lrr_ipr.dissect.txt";
my $tmp_tir_ipr                   = $prefix."tmp_tir_ipr.dissect.txt";
my $tmp_cc_ipr                    = $prefix."tmp_cc_ipr.dissect.txt ";  
my $error_report                  = $prefix."Error.logfile.txt";

# ----------------preprocessing protein/DNA fasta sequence-----------------------
Ptime("Pipeline to predict the plant RGA...");
Ptime("Make sure all other programs are ready...");
open(ERROR, ">$error_report");

local $/ = ">";

Ptime("formatting the input file as standard input file...");
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
Ptime("Blast with RGA DB...");
if ($RGA_blast_out and -s $RGA_blast_out) {
    Ptime("$RGA_blast_out detected in current folder, pipeline will jumps to next step - code 003");
}
else {
    blastp_parallel($aa_formated_infile, $RGD_index_file, $blast_evalue, $RGA_blast_out, $cpu);
}

Ptime("Blast is done, now parsing the output...");

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
    Ptime("$pfam_out detected in current folder, pipeline will jumps to next step - code 001");
}
else {
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

#-----------selectively output all pfam scan data for NBS(true) coding gene only
open(IN,  "$pfam_out");
open(OUT, ">$NBS_pfam_out");
while (<IN>) {
    chomp;
    next unless ($_ =~ /\w/ or $_ =~ /\d/);
    my ($geneid, @array) = split/\s+/,$_;
    if ($NBS_pfam_lst{$geneid}) {
        print OUT "$_\n";                    
    }
}
close IN;
close OUT;

# --------------------coiled coil prediction -------------------
Ptime("to predict coiled coils...");
if ($cc_prediction and -s $cc_prediction) {
    Ptime("$cc_prediction detected in current folder, pipeline will jumps to next step - code 002");
}
else {
    Ptime("initializing coiled-coil prediction...");
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
    Ptime("$iprscan_out detected in current folder, pipeline will jumps to next step - code 004");
}
else {
    Ptime("initializing interproscan...");
    system("interproscan.sh -i $RGA_blast_fasta -appl Pfam,panther,smart,gene3d,superfamily -f tsv -iprlookup -o $iprscan_out 1>/dev/null");
}

Ptime("Interproscan is done...");

# --------------------------RLK and RLP prediction--------------
if ($RLKorRLP_prediction_output and -s $RLKorRLP_prediction_output) {#$RLKRLP_out_raw
    Ptime("$RLKorRLP_prediction_output detected in current folder, pipeline will jumps to next step - code 005");
}
else {
     system("perl -S RLK.prediction.pl -i $RGA_blast_fasta -pfx $prefix -pfam $pfam_out               -iprs $iprscan_out -cpu $cpu -lst $RGA_blast_lst     -o $RLKorRLP_prediction_output");
}

#-----------merge coilsed coil to above output<$RLKorRLP_prediciton_outputs---------------
open(IN,      $RLKorRLP_prediction_output);
open(MERGE,     ">$RLKorRLP_merged_domain");
print MERGE join("\t", "id", "stk", "tm", "sp", "LysM", "LRR", "CC\n");
#this will add one more column in terms of cc to the $RLKorRLP_prediction_output
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
system("perl -S pfamscan.RGA.summary.pl        -i   $NBS_pfam_out      -pfx $prefix");              

# -extra leucine rich repeat analysis --------------
extra_LRR_analysis("$lrr_prediction");

#merge essential motif/domain to one file
system("perl -S nbs.domain.result.merge.pl     -nbs $nbs_prediction    -lrr $lrr_prediction         -tir $tir_prediction -cc $cc_prediction -seq $aa_formated_infile >$NBS_merged_domain");  

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
open(TMP_NBS,">$tmp_nbsonly_fas");
push(@deletion,$tmp_nbsonly_fas) if (-e $tmp_nbsonly_fas);
#push(@deletion,$NBS_pre_candidates_lst);
#push(@deletion,"summary.txt");

while (<IN>) {
    chomp;
    my ($id,$type) = split/\t/,$_;
    if ($type eq 'NBS') {
        # for those NBS type, it will be undertaken further analysis.
        print TMP_NBS ">$id\n$protein_fasta{$id}\n";
    }
    else {
        # for those not NBS tupe, they will be hashed as final NBS-encoding candidates
        $NBS_candidates_lst{$id} = $type;
    }
}
close IN;
close TMP_NBS; 

if ($iprscan_out_2nd and -s $iprscan_out_2nd) {
    Ptime("$iprscan_out_2nd detected in current folder, pipeline will jumps to next step - code 006");
}
else {
    Ptime("initializing interproscan 2nd round...");
    system("interproscan.sh -i $tmp_nbsonly_fas -appl pfam,superfamily,coils -f tsv -iprlookup -o $iprscan_out_2nd 1>/dev/null");
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

Ptime("Interproscan further analysis is done...");


#----output lst of RLK, RLP and TMCC---
Ptime("RLK and RLP prediction is done...");
system("perl -S RLK.prediction.result.parser.v2.pl $RLKorRLP_merged_domain $NBS_candidates_lst $RLK_candidates_lst $RLP_candidates_lst $TMCC_candidates_lst");


# ----------output RGA candidates aa sequence --------------
my %id = ();
open(OUT,">$RGA_candidates_fasta");
foreach my $lst ($NBS_candidates_lst, $RLK_candidates_lst, $RLP_candidates_lst, $TMCC_candidates_lst) {
    open(IN,"$lst");
    while (<IN>) {
        chomp;
        next if (/^\s/);
        my @array = split/\s/,$_;

        my $id   = $array[0];
        my $type = $array[1]; 
        
        my $new_id = join("|",$id,$type);
        print OUT ">$new_id\n$protein_fasta{$id}\n";
        
        # ------record outputed RGA --------
        $id{$id} = 1;
    }
    close IN;
}
close OUT;

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
            print ERROR "failed to ouptut nucleotide sequence: $id \n";
        }
    }
    close OUT;
}


#------------------prepare CVIT data-=----------------------
if ($gff and -s $gff) {
    system("perl -S cvit.input.generator.pl -l $NBS_candidates_lst  -p $aa_infile -f $gff -t gene -c blue   -t2 NBS  ");
    system("perl -S cvit.input.generator.pl -l $RLK_candidates_lst  -p $aa_infile -f $gff -t gene -c green  -t2 RLK  ");
    system("perl -S cvit.input.generator.pl -l $RLP_candidates_lst  -p $aa_infile -f $gff -t gene -c orange -t2 RLP  ");
    system("perl -S cvit.input.generator.pl -l $TMCC_candidates_lst -p $aa_infile -f $gff -t gene -c black  -t2 TMCC ");
}

#------------------ clean files ----------------------------
push(@deletion, $error_report) if (-z $error_report);  #remove Error.log if its' empty logged.
close ERROR;

foreach my $file (@deletion) {
    unlink "$file" or warn "couldnt delete $file: $!\n";
}

my $end_run = time();
hhmmss_consumed($end_run - $start_run);


#-----------------------------------------------------------
#--------------------------sub functions--------------------
#-----------------------------------------------------------
sub Ptime{
          my $time = localtime;
          my ($msg)= @_;
          print STDERR "$time: $msg\n";
}

sub generate_rand_filename {

     my $length_of_randomstring=shift;

     my @chars=('a'..'z','A'..'Z','0'..'9');
     my $random_string;
     foreach (1..$length_of_randomstring)
     {
          # rand @chars will generate a random
          # number between 0 and scalar @chars
          $random_string.=$chars[rand @chars];
     }
     return $random_string;
}

sub fasta_file_split {
    my ($file, $thread) = @_;
    my @splited_files = ();
    my $i = 1;
    my $fileNumber = 1;
    local $/ = ">";
    
    my $seq_number       = `grep ">" $file |wc -l`; $seq_number =~ s/\s+//g;
    my $each_file_number = abs($seq_number/$thread) + 1;
    
    push(@splited_files,".splited_file_$fileNumber.txt");
    
    open(IN,$file);
    open (OUT,">.splited_file_$fileNumber.txt");
    while (<IN>) {
        chomp;
        my ($title,$seq) = split/\n/,$_,2;
        next unless ($title && $seq);
        
        $seq   =~ s/\s+//g;
        
        if ($i <= ($fileNumber*$each_file_number)) {
            print OUT ">$title\n$seq\n";
            $i++;
        }
        else {
            close OUT;
            $fileNumber++;
            push(@splited_files,".splited_file_$fileNumber.txt");
            
            open(OUT,">.splited_file_$fileNumber.txt");
            print  OUT ">$title\n$seq\n";
            $i++;
        }
    }
    close IN;
    close OUT;
    local $/ = "\n";
    return(@splited_files);
}

sub files_remove{
    my @files = @_;
    foreach my $file (@files) {
        unlink "$file" or warn "could not remove $file: $!\n";
    }
}

sub threads_finish_examination{
    my ($userID, @keywords) = @_;  #keywords is a key word that can distinguish it from other unrunning return from ps -ef
    while (1) {
        #print "ps -u $userID -f\n";
        my @commands = `ps -f -u $userID`;
        my $flag = 0;
        
        foreach my $line (@commands) {
               foreach my $keyword (@keywords) {
                   if ($line =~ /$keyword/i) {
                       $flag++;
                   }
               }
        }
        
        #print STDERR "flag = $flag\n";
        if ($flag>0) {
            sleep 10;
        }
        else {
            #all keywords can ben detected from @commands.
            last;
        }
    }
}


sub splitted_results_merge {
    my ($output,@splitted_out) = @_;
    
    open(OUT,">$output");
    foreach my $file (@splitted_out) {
        if (-s $file) {
            open(IN,$file);
            while (<IN>) {
                chomp;
                print OUT "$_\n";
            }
            close IN;
        }
        else {
            print STDERR "$file size is zero\n";
        }
    }
    close OUT;
}

sub pfamscan_parallel {
    my ($pfam_out, $e_seq, $e_dom, $cpu, $RGA_blast_fasta, $pfam_index_folder) = @_;
    my @splitted_out = ();
    my $userID = `echo \$USER`;
    my @fingerprints = ();
    
    my @split_files = fasta_file_split($RGA_blast_fasta, $cpu);
    foreach my $file (@split_files) {
        if (-e $file) {
            my $fingerprint = generate_rand_filename(10); #to replace PID to monitor the status of threads
            push(@fingerprints,$fingerprint);
            system("perl -S pfam_scan.pl -outfile .splited.$fingerprint.pfam$file.out -e_seq $e_seq -e_dom $e_dom --cpu 1 -as -fasta $file -dir $pfam_index_folder &");
            push(@splitted_out,".splited.$fingerprint.pfam$file.out");
        }else {
            print STDERR "unable to find splitted infile $file\n";
        }
    }
    
    #below section is to check if all pfam_scan.pl threads are done.
    threads_finish_examination($userID,@fingerprints);

    #splitted result merge to an intact output $tmp_file
    splitted_results_merge($pfam_out,@splitted_out);

    files_remove(@splitted_out);
    files_remove(@split_files);
}

sub blastp_parallel {
    my ($aa_formated_infile,$RGD_index_file,$blast_evalue,$RGA_blast_out,$cpu) = @_;
    my @splitted_out = ();
    my $userID = `echo \$USER`;
    my @fingerprints = ();
    
    my @split_files = fasta_file_split($aa_formated_infile, $cpu);
    foreach my $file (@split_files) {
        if (-e $file) {
            my $fingerprint = generate_rand_filename(10); #to replace PID to monitor the status of threads
            push(@fingerprints,$fingerprint);
            system("blastp -task blastp -query $file -db $RGD_index_file -evalue $blast_evalue -out .splited.$fingerprint.blast$file.out -num_threads 1 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle\" &");
            push(@splitted_out,".splited.$fingerprint.blast$file.out");
        }else {
            print STDERR "unable to find splitted infile $file\n";
        }
    }
    
    #below section is to check if all pfam_scan.pl threads are done.
    threads_finish_examination($userID, @fingerprints);

    #splitted result merge to an intact output $tmp_file
    splitted_results_merge($RGA_blast_out,@splitted_out);

    files_remove(@splitted_out);
    files_remove(@split_files);
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
            Ptime("not found $id while outputting fasta");
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
        Ptime("unalbe to output @error in $line");
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
        Ptime("$output existed in current folder, reading it to RAM and going ahead with next step...");
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
  Ptime("input <$aa_infile> time taken -  $consumed");
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
