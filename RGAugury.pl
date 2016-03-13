#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use File::Basename;
use File::Path qw(make_path remove_tree);

#RGAugury pipeline

Ptime("Pipeline to predict the plant RGA...");
Ptime("Make sure all other programs are ready...");


my $version = 0.2;

my $USAGE = <<USAGE;
Scripts: RGA prediction pipeline

Version: $version, Written by Pingchuan Li, Frank You Lab

Usage :perl RGAs.prediction.pl <option>

arguments: 

        -p           protein fasta file
        -n           corresponding cDNA/CDS nucleotide for -p (optional)
        -g           genome file in fasta format (optional)
        -gff         gff3 file (optional)
        -c           cpu or threads number, default = 2
        -pfx         prefix for filename, useful for multiple speices input in same folder (optional)
        
USAGE

#---------- parameters which should be defined in advance -----------
GetOptions(my $options = {},
                    "-p=s","-n=s","-g=s","-gff=s","-c=i","-pfx=s"
                  );

die $USAGE unless defined ($options->{p});

my $aa_infile   = $options->{p};  
my $nt_infile   = $options->{n};  
my $g_infile    = $options->{g};
my $gff         = $options->{gff};
my $cpu         = $options->{c};
my ($prefix)    = ($options->{pfx}) ? "$options->{pfx}" : $aa_infile =~ /([a-zA-Z0-9]+)/; $prefix .= ".";

# --------------------- cutoff, key interproScan ID and blastp evalue for pfamScan ---------------------------
my $e_seq    = 0.1;
my $e_dom    = 0.1;
my @interested_NBS = ("PF00931"); # incase of multiple domains, array is used here. make sure all domain code starting with 'PF'
my $blast_evalue = "1e-5";

#make sure below folder contain pfam and preselected RGA database
my $pfam_index_folder = (-e $ENV{"HOME"}."/database/pfam.v27") ? $ENV{"HOME"}."/database/pfam.v27": die "unable to locate pfam DB";
my $RGD_index_file    = (-e $ENV{"HOME"}."/database/RGADB/plant.RGA.dataset.unique.fasta") ? $ENV{"HOME"}."/database/RGADB/plant.RGA.dataset.unique.fasta" : die "unalbe to locate RGADB file";

# -------------------  main body -----------------------------
my %NBS_pfam_lst               = ();
my %RGA_blast_lst              = ();
my %protein_fasta              = ();
my %nt_fasta                   = ();
my %genome_fasta               = ();
my %overlap_RGAblast_pfam_lst  = ();
my %NBS_candidates_lst         = ();
my @deletion                   = ();
my %coils = ();



#otput file name
my $aa_formated_infile            = $prefix."formated.protein.input.fas";
my $pfam_out                      = $prefix."pfam.local.search.out";
my $NBS_pfam_out                  = $prefix."NBS.pfam.out";
my $NBS_pre_candidates_lst        = $prefix."NBS.pre.candidates.lst";
my $NBS_candidates_lst            = $prefix."NBS.candidates.lst";
my $NBS_candidates_fas            = $prefix."NBS.candidates.fas";
my $NBS_merged_domain             = $prefix."NBS.merged.domains.txt";

my $RGA_blast                     = $prefix."RGA.blastp.$blast_evalue.out";
my $RGA_blast_fasta               = $prefix."tmp.RGA.blast.lst.fasta";
my $RGA_blast_lst                 = $prefix."tmp.RGA.blast.lst";
my $RGA_candidates_fasta          = $prefix."RGA.candidates.fasta";
my $RGA_candidates_fasta_nt       = $prefix."RGA.candidates.cdna.fasta";
my $candidate_RGA_pfam_out        = $prefix."candidates_RGA_pfam_out";
my $iprscan_out                   = $prefix."iprscan_out.tsv";
my $iprscan_out_2nd               = $prefix."iprscan_out_further.tsv";

my $nbs_prediction                = $prefix."NBS.res.pfam.txt";   #NBS.res.pfam.txt -lrr LRR.res.pfam.txt -tir TIR.res.pfam.txt 
my $lrr_prediction                = $prefix."LRR.res.pfam.txt";
my $tir_prediction                = $prefix."TIR.res.pfam.txt";

my $cc_prediction                 = $prefix."coils.res.txt";
my $candidate_RGA_lst             = $prefix."candidates_RGA_lst";
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


# ----------------read protein fasta--------------------------
open(ERROR,">$error_report");

Ptime("formatting the input file...");
local $/ = ">";
open(IN, $aa_infile) or die "unalbe to open $aa_infile\n";
while (<IN>) {
    chomp;
    my ($title,$seq) = split/\n/,$_,2;
    next unless ($title and $seq);
    my ($id) = $title =~ /([a-zA-Z0-9\.\-\_]+)/;
    #$title =~ s/\s+//g;
    $seq   =~ s/\s+//g;
    $seq   =~ s/\*//g;
    
    $protein_fasta{$id}= $seq;
}
close IN;

open(OUT,">$aa_formated_infile");
foreach my $id (sort {$a cmp $b} keys %protein_fasta) {
    print OUT ">$id\n$protein_fasta{$id}\n";
}
close OUT;

if ($nt_infile and -s $nt_infile) {
    open(IN,$nt_infile);
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
    open(IN, $g_infile);
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

#-----------------using pfam to scan the pfam-------------------
#--------------supporting broken execution----------------------
#--------remove the $pfam_out if not correctly executed --------

if ($pfam_out and -s $pfam_out) {
    Ptime("$pfam_out detected in current folder, pipeline will jumps to next step - code 001");
}
else {
    pfamscan_parallel($pfam_out, $e_seq, $e_dom, $cpu, $aa_formated_infile, $pfam_index_folder);
}

open(IN, $pfam_out) or die "cant open $pfam_out file\n";  #use @interested_NBS to screen pfam out and acquire all the NBS-encoding genes list
while (<IN>) {
    chomp;
    next if ($_ =~ /^#/ or $_ =~ /^\s/);
    
    my ($geneid,@array) = split/\s+/,$_;
    foreach my $PF (@interested_NBS) {
        if ($_ =~ /$PF\W/) {
            $NBS_pfam_lst{$geneid} = 1;
        }
    }
}
close IN;

#-----------selectively output all pfam scan data for NBS(true) coding gene only
open(IN,"$pfam_out");
open(OUT,">$NBS_pfam_out");
while (<IN>) {
    chomp;
    next unless ($_ =~ /\w/ or $_ =~ /\d/);
    my ($geneid,@array) = split/\s+/,$_;
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
    system("perl -S coils.identification.pl $aa_formated_infile $cpu $cc_prediction");
}

open(IN, $cc_prediction);  #$cc_prediction is output of coils.identification.pl
while (<IN>) {
   chomp;
   my ($id,@res) = split/\t/,$_;
   $coils{$id} = join(" ",@res);
}
close IN;

# -------------------------blastp to RGA database------------------
# -----------this step will get a potential candidates RGA---------
Ptime("Blast with RGA DB...");
if ($RGA_blast and -s $RGA_blast) {
    Ptime("$RGA_blast detected in current folder, pipeline will jumps to next step - code 003");
}
else {
    blastp_parallel($aa_formated_infile, $RGD_index_file, $blast_evalue, $RGA_blast, $cpu);
}

Ptime("Blast is done, now parsing the output...");

open(IN,$RGA_blast);
while (<IN>) {
    chomp;
    my ($geneid,@array) = split/\t/,$_;
    $RGA_blast_lst{$geneid} = 1;
}
close IN;

#--------output candidates of RGA, which include NBS-encoding, RLP, RLK and other potential disease resistance genes analogs.-----------
output_protein_fasta_ram_manner(\%RGA_blast_lst, \%protein_fasta, $RGA_blast_fasta, __LINE__);
push(@deletion,"$RGA_blast_fasta");

#open(OUT,">$RGA_blast_fasta");          
#foreach my $id (sort {$a cmp $b} keys %RGA_blast_lst) {
#    my $seq = $protein_fasta{$id};
#    print OUT ">$id\n$seq\n";
#}
#close OUT;

# -------------------------iprscan---------------------------------
#----using above outputed fasta file to do iprscan for 1st time----
if ($iprscan_out and -s $iprscan_out) {
    Ptime("$iprscan_out detected in current folder, pipeline will jumps to next step - code 004");
}
else {
    Ptime("initializing interproscan...");
    system("interproscan.sh -i $RGA_blast_fasta -appl Pfam -f tsv -iprlookup -o $iprscan_out 1>/dev/null");
}

Ptime("Interproscan is done...");

# -----------extract only RGA related pfam info from $pfam_out-----------
open(IN, $pfam_out) or die "cant open $pfam_out";
open(OUT,">$candidate_RGA_pfam_out");
#push(@deletion,$candidate_RGA_pfam_out);
while (<IN>) {
    chomp;
    next unless ($_);
    my ($geneid,@array) = split/\s+/, $_;
    if (exists $RGA_blast_lst{$geneid}) {
        print OUT join("\t",$geneid,@array);
        print OUT "\n";
        
        $overlap_RGAblast_pfam_lst{$geneid} = 1;  #select pfam_scan output from those which has hits in RGA DB hits
    }
}
close IN;
close OUT;

open(LST,">$candidate_RGA_lst");
foreach my $id (sort {$a cmp $b} keys %overlap_RGAblast_pfam_lst) {
    print LST "$id\n";
}
close LST;

# --------------------------RLK and RLP prediction--------------

if ($RLKorRLP_prediction_output and -s $RLKorRLP_prediction_output) {#$RLKRLP_out_raw
    Ptime("$RLKorRLP_prediction_output detected in current folder, pipeline will jumps to next step - code 005");
}
else {
    system("perl -S RLK.prediction.pl -i $aa_formated_infile -pfx $prefix -pfam $candidate_RGA_pfam_out -iprs $iprscan_out -cpu $cpu -lst $candidate_RGA_lst -o $RLKorRLP_prediction_output");
}

#-----------merge coilsed coil to above output<$RLKorRLP_prediciton_outputs---------------
open(IN,  $RLKorRLP_prediction_output);
open(TMP,     ">$RLKorRLP_merged_domain");
print TMP join("\t", "id","stk","tm","sp","LysM","LRR","CC\n");
#this will add one more column in terms of cc to the $RLKorRLP_prediction_output
while (<IN>) {
   chomp;
   next if($_ =~ /LysM\tLRR/); #ignore first headline prior to processing
   my ($id,@content) = split/\t/,$_;
   my $cc = ($coils{$id}) ? $coils{$id} : '.' ;
   print TMP join("\t", $id, @content, "$cc\n");
}
close IN;
close TMP;


# ----------dissect NBS.pfam.out------------generate NBS.res.pfam.out etc.------------
system("perl -S pfamscan.RGA.summary.pl -i $NBS_pfam_out -pfx $prefix");#output  NBS.res.pfam.txt LRR.res.pfam.txt, TIR.res.pfam.txt and PPR.res.pfam.txt totall 4 files
system("perl -S nbs.domain.result.merge.pl -nbs $nbs_prediction -lrr $lrr_prediction -tir $tir_prediction -cc $cc_prediction -seq $aa_formated_infile >$NBS_merged_domain");  #all $domain_prediction are output of last scirpt
system("perl -S NBS-encoding.amount.summary.pl -i $NBS_merged_domain -o $NBS_pre_candidates_lst -pfx $prefix");

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
push(@deletion,"intermediate.txt") if (-e "intermediate.txt");
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
    system("interproscan.sh -i $tmp_nbsonly_fas -appl pfam,superfamily,panther,coils -f tsv -iprlookup -o $iprscan_out_2nd 1>/dev/null");
}

system("perl -S ipr.specific.id.selection.pl $iprscan_out_2nd $tmp_nbs_ipr $tmp_lrr_ipr $tmp_tir_ipr $tmp_cc_ipr") if ($iprscan_out_2nd and -s $iprscan_out_2nd);#keep the order of output

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


#output lst of RLK, RLP and TMCC
#system("perl -S RLK.prediction.result.parser.v2.pl $RLKorRLP_merged_domain $RLK_candidates_lst $RLP_candidates_lst $TMCC_candidates_lst");
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
if ($gff and $g_infile) {
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
    my ($pfam_out, $e_seq, $e_dom, $cpu, $aa_formated_infile, $pfam_index_folder) = @_;
    my @splitted_out = ();
    my $userID = `echo \$USER`;
    my @fingerprints = ();
    
    my @split_files = fasta_file_split($aa_formated_infile, $cpu);
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
    my ($aa_formated_infile,$RGD_index_file,$blast_evalue,$RGA_blast,$cpu) = @_;
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
    splitted_results_merge($RGA_blast,@splitted_out);

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