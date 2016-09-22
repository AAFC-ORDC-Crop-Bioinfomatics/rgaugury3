package Tool;
use strict;
use warnings;
use Exporter qw(import);
 
our @EXPORT_OK = qw(parallel blastp_parallel files_remove pfamscan_parallel transmembrane iprscan coils_parallel splitted_results_merge);

my $BLASTP = "blastp";
my $PFAMSCAN = "pfamscan";
my $TRANSMEMBRANCE = "transmembrane";
my $COIL = "coil";
my $IPSCAN= "iprscan";

sub parallel {
  my ($task,$input_file,$cpu,$cmd,$output_file) = @_;
  my @split_files = fasta_file_split($input_file, $cpu);
  my @pids=();
  my @splitted_out = ();
  my $pid;
  
  foreach my $file (@split_files) {
      if (-e $file) {
          my $output = ".splited.$task$file.out";
          my $find = "{in}";
          my $splited_cmd = $cmd;
          $find = quotemeta $find; # escape regex metachars if present
          $splited_cmd =~ s/$find/$file/g;
          $find = "{out}";
          $find = quotemeta $find; # escape regex metachars if present
          $splited_cmd =~ s/$find/$output/g;
          if ($pid =fork){
              #parent process here
              push (@pids,$pid);
          } else {
              #child process here
              exec($splited_cmd);
          }
          push(@splitted_out,$output);
      }else {
          print STDERR "unable to find splitted infile $file\n";
      }
  }
  
  waitpid($_,0) for @pids;  
  #splitted result merge to an intact output $tmp_file
  splitted_results_merge($output_file,@splitted_out);
  files_remove(@splitted_out);
  files_remove(@split_files);
} 

sub blastp_parallel {
  my ($aa_formated_infile,$RGD_index_file,$blast_evalue,$RGA_blast_out,$cpu) = @_;
  my $cmd ="blastp -task blastp -query {in} -db $RGD_index_file -evalue $blast_evalue ";
  $cmd .= "-out {out} -num_threads 1 -outfmt '6 qseqid sseqid pident length mismatch ";
  $cmd .= "gapopen qstart qend sstart send evalue bitscore qlen slen stitle'";

  parallel($BLASTP, $aa_formated_infile,$cpu,$cmd,$RGA_blast_out);  
}

sub pfamscan_parallel {
    my ($pfam_out, $e_seq, $e_dom, $cpu, $RGA_blast_fasta, $pfam_index_folder) = @_;
    my $cmd = "perl -S pfam_scan.pl -outfile {out} -e_seq $e_seq -e_dom $e_dom --cpu 1 -as -fasta {in} -dir $pfam_index_folder";            
    parallel($PFAMSCAN, $RGA_blast_fasta, $cpu,$cmd, $pfam_out);
}

sub transmembrane {
    my ($input, $output, $cpu) = @_;
    my $cmd = "phobius.pl -short {in} > {out} 2>/dev/null ";
    parallel($TRANSMEMBRANCE, $input, $cpu,$cmd,$output);
}

sub iprscan{
    my ($input, $output, $cpu, $appl, $format, $logfile) = @_;
    my $cmd = "interproscan.sh -i {in} -appl $appl -f $format -iprlookup -o {out} 1>>./$logfile";
    parallel($IPSCAN, $input, $cpu,$cmd,$output);
}

sub coils_parallel{
    my ($input, $output, $cpu) = @_;
    my $cmd = "scoils-ht -f {in} > {out} 2>/dev/null";
    parallel($COIL, $input, $cpu,$cmd,$output);
}

sub splitted_results_merge {
    my ($output,@splitted_out) = @_;
    local $/ = "\n";
    open(OUT,">$output");
    foreach my $file (@splitted_out) {
        if (-s $file) {
            open(IN,'<'.$file);
            while (<IN>) {
              chomp;
              next if ($_ =~ /^\s*$/ or $_ =~ /^#/);
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


sub fasta_file_split {
    my ($file, $thread) = @_;
    my @splited_files = ();
    my $i = 1;
    my $fileNumber = 1;
    local $/ = ">";
    my $seq_number = `grep ">" $file |wc -l`; 
    $seq_number =~ s/\s+//g;
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

1;

