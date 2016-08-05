#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Log::Log4perl::Level;
use Log::Log4perl qw(:easy);

use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0);
use Tool qw(iprscan);

#owning to iprlookup limitation occasionally, interproscan was invoked by perl script iprscan.pl to run on splitted small data set


GetOptions(my $options = {},
              "-i=s","-appl=s","-f=s","-o=s","-log=s"
);


my $USAGE = <<USAGE;
-------------------------------------------------------------------------------------------------------------
Script: to do iprscan batchly

Usage:   perl $0 <options> <files>
-------------------------------------------------------------------------------------------------------------
 script was coded by Pingchuan Li @ You Lab
Arguments:
        -i          protein input
        -appl       dataset
        -f          output format
        -o          output filename
        -log        log filename

USAGE

die $USAGE unless(defined $options->{i});

my $fasta                 = $options->{i};
my $output                = $options->{o};
my $splitted_files_number = count_split_number($fasta, 4000);  #each splitted file contain 4000 entries
my $logfile               = $options->{log};
my $appl                  = $options->{appl};
my $format                = $options->{f};

iprscan($fasta,  $output, $splitted_files_number, $appl, $format, $logfile);

# -------------------------sub ------------------------------
sub count_split_number{
    my ($file,$NoEach) = @_;
    my $totalSeqNo = `grep ">" $file |wc -l`; $totalSeqNo =~ s/\s+//g;

    my $number = '';
    
    my $value = ($totalSeqNo%$NoEach);
    if ($value == 0) {
        $number = $totalSeqNo/$NoEach;
    }
    else {
        $number = 1 + int($totalSeqNo/$NoEach);
    }
    return $number;
}
=for comment
sub iprscan{
    my ($input, $output)    = @_;
    my @splitted_out        = ();
    my @fingerprints        = ();
    my @renamed_split_files = ();
    my $userID              = `echo \$USER`;

    my @split_files = fasta_file_split($input, $splitted_files_number);

    foreach my $file (@split_files) {
        if (-e $file) {
            my $fingerprint = generate_rand_filename(12); #to replace PID to monitor the status of threads
            push(@fingerprints, $fingerprint);
            
            rename $file, ".splitted.$fingerprint.file.txt"; #because > content won't display in ps, thus extra steps were added here.
            push(@renamed_split_files, ".splitted.$fingerprint.file.txt");
            system("interproscan.sh -i .splitted.$fingerprint.file.txt -appl $appl -f $format -iprlookup -o .splited.res.ipr$file.out 1>>./$logfile");
            push(@splitted_out,".splited.res.ipr$file.out");
        }
        else {
            print STDERR "unable to find splitted infile $file\n";
        }
    }
    
    #below section is to check if all pfam_scan.pl threads are done.
    threads_finish_examination($userID, @fingerprints);

    #splitted result merge to an intact output $tmp_file
    splitted_results_merge($output, @splitted_out);

    files_remove(@splitted_out);
    files_remove(@renamed_split_files);
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
    local $/ = "\n";
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
=cut
