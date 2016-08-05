#!/usr/bin/perl -w
use strict;
use Scalar::Util qw(looks_like_number);

use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0);
use Tool qw(coils_parallel);

die "perl coils.identification.pl <input_fasta> <cpu_number> <output_file>" unless ($#ARGV == 2);

my $input_fasta = $ARGV[0];
my $cpu         = (looks_like_number($ARGV[1])) ? $ARGV[1] : 2 ;
my $out         = ($ARGV[2]) ? $ARGV[2] : "coils.res.txt" ;

local $/ = ">";
open(IN, $input_fasta);
open(OUT,">.temporary1.txt");
while (<IN>) {
    chomp;
    my ($title,$seq) = split/\n/,$_,2;
    next unless ($title and $seq);
    $seq =~ s/\s+//g;
    $seq =~ uc($seq);
    print OUT ">$title\n$seq\n";
}

close IN;
close OUT;

coils_parallel(".temporary1.txt", ".temporary2.txt", $cpu);

open(IN,".temporary2.txt");
my %gene = ();
while (<IN>) {
    chomp;
    my ($title, $seq) = split/\n/,$_,2;
    next unless ($title and $seq);
    
    $seq =~ s/\s+//g;
    my ($id) = $title =~ /(\S+?)\s/;
    my $flag = 0;
    while ($seq =~ /x+/g) {
        my $start = $-[0] + 1;
        my $end   = $+[0];
        push(@{$gene{$id}},join("|","domain_CC","$start-$end"));
        $flag++;
    }
    
    if ($flag == 0) {
        push(@{$gene{$id}},"." );
    }
}
close IN;
open(OUT,">$out");
foreach my $gene (sort {$a cmp $b} keys %gene) {
    print OUT "$gene\t";
    my $value = join("\t",@{$gene{$gene}});
    print OUT "$value\n";
}
close OUT;

unlink".temporary1.txt";
unlink".temporary2.txt";

#-------------------------------sub--------------------------------
sub Ptime{
          my $time = localtime;
          my ($msg)= @_;
          print STDERR "$time: $msg\n";
}
=for comment
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

sub coils_parallel{
    my ($input, $output) = @_;
    my @splitted_out = ();
    my $userID = `echo \$USER`;
    #Ptime("userid = $userID");
    my @fingerprints = ();
    
    my @split_files = fasta_file_split($input, $cpu);
    my @renamed_split_files = ();
    
    foreach my $file (@split_files) {
        if (-e $file) {
            my $fingerprint = generate_rand_filename(12); #to replace PID to monitor the status of threads
            push(@fingerprints,$fingerprint);
            
            rename $file, ".splitted.$fingerprint.file.txt"; #because > content won't display in ps, thus extra steps were added here.
            push(@renamed_split_files, ".splitted.$fingerprint.file.txt");
            system("scoils-ht -f .splitted.$fingerprint.file.txt >.splited.res.coils$file.out 2>/dev/null &");
            push(@splitted_out,".splited.res.coils$file.out");
        }
        else {
            print STDERR "unable to find splitted infile $file\n";
        }
    }
    
    #below section is to check if all pfam_scan.pl threads are done.
    threads_finish_examination($userID,@fingerprints);

    #splitted result merge to an intact output $tmp_file
    splitted_results_merge($output, @splitted_out);

    files_remove(@splitted_out);
    files_remove(@renamed_split_files);
}
=cut


