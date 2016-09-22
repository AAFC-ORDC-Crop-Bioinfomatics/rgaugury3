#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use FindBin;

use File::Basename qw(dirname);
use Cwd  qw(abs_path);
use lib dirname(abs_path $0);
use Tool qw(transmembrane);

#print STDERR "This script is yet to be finished, the RLK and RLP has minor redundant with NBS encoding genes\n";

my $USAGE = <<USAGE;
Scripts: RLK protein prediction, make sure the Phobius can be invoked.

 Written by Pingchuan Li, Frank You Lab

Usage :perl unique.xxx.pl <option> <file>

arguments: 

        -i           fasta input file
        -pfam        pfam result file
        -iprs        interproScan result file
        -lst         if lst file was specified, only the specified gene
                     will be output, one line one gene name
        -o           output filename
        -pfx         prefix for output file
        -cpu         cpu numbers used to analysis

USAGE


#---------- parameters which should be defined in advance -----------
GetOptions(my $options = {},
                    "-i=s","-pfam=s","-iprs=s","-lst=s", "-o=s","-cpu=i","-pfx=s"
                  );
die $USAGE unless defined ($options->{i});

my $input      = $options->{i};
my @types      = ("stk","tm","sp","LysM","LRR"); #tm = transmembrane
my %result     = ();
my @tmp        = ();
my $prefix     = ($options->{pfx}) ? $options->{pfx} : "output." ;
my $outputfile = ($options->{o}) ? $options->{o} : "RLKorP.assay.Alll.genes.lst.txt"; # or RLKorRLP.domain.prediction.txt specified by RGA pipeline

my $cpu            = ($options->{cpu}) ? $options->{cpu} : 1 ;
my $phobius_output = "$prefix"."phobius.txt";
#push(@tmp, "$phobius_output");

my @pfam_rlk       = keys %{file_parser("$FindBin::Bin/configuration_rlk.txt")};

transmembrane($input, $phobius_output, $cpu);

# --------------data processing-----------

my $i = 0;
open(IN,"$phobius_output") || die "cant open the $input file";
while (<IN>) {
    chomp;
    $i++;
    my @array = split/\s+/,$_;
    next if ($array[0] eq 'SEQENCE' and $array[1] eq 'ID');
    
    my $id = $array[0];

    if ($array[1] >= 1 and $array[2] =~ /y/i) {
        $result{$id}->{tm}->{property} = 1;
        $result{$id}->{sp}->{property} = 1;
        
        my ($sp_coor,$tm_coor) = split/\//,$array[3];
        my ($sp_start,$sp_end ) = $sp_coor =~ /(\d+)-(\d+)/;
        push(@{$result{$id}->{sp}->{range}},join("|","phobius","$sp_start-$sp_end"));
        
        while ($tm_coor =~ /(\d+)-(\d+)/g) {
            my ($tm_start,$tm_end) = ($1,$2);
            push(@{$result{$id}->{tm}->{range}},join("|","phobius","$tm_start-$tm_end"));
        }
    }
    elsif ($array[1] >= 1) {
        $result{$id}->{tm}->{property} = 1;
        while ($array[3] =~ /(\d+)-(\d+)/g) {
            my ($tm_start,$tm_end) = ($1,$2);
            push(@{$result{$id}->{tm}->{range}},join("|","phobius","$tm_start-$tm_end"));
        }
    }
}
close IN;
# -------------------------

# merged pfam and iprscan outcome as one %result
if (-e $options->{pfam}) {
    pfam_parse($options->{pfam});
}
else {
    Ptime("not use -pfam file in RLK.prediction.pl");
}

if (-e $options->{iprs}) {
    iprscan_parse($options->{iprs});
}
else {
    Ptime("not use -iprs file in RLK.prediction.pl");
}

# ---------output --------

open(ALL,">$outputfile");
print ALL join("\t","id","stk","tm","sp","LysM","LRR\n");
if ($options->{lst}) {#only processing specified RGA lst
    my $lst = lst_processing($options->{lst});
    
    my %lst = %{$lst};
    
    foreach my $id (sort {$a cmp $b} keys %lst) {
        print ALL "$id";
        foreach my $type (@types) {
            if (looks_like_number($result{$id}->{$type}->{property})) {
                my $value = join(" ",@{$result{$id}->{$type}->{range}});
                print ALL "\t$value";
            }
            else {
                print ALL "\t.";
            }
        }
        print ALL "\n";
    }
}
else {
    foreach my $id (sort {$a cmp $b} keys %result) {
        print ALL "$id";
        foreach my $type (@types) {
            if (looks_like_number($result{$id}->{$type}->{property})) {
                my $value = join(" ",@{$result{$id}->{$type}->{range}});
                print ALL "\t$value";
            }
            else {
                print ALL "\t.";
            }
        }
        print ALL "\n";
    }
}


foreach my $file (@tmp) {
    unlink "$file";
}
close ALL;

# -----------sub--------------
=for comment
sub transmembrane {
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
            
            system("phobius.pl -short \".splitted.$fingerprint.file.txt\">.splited.res.phobious$file.out 2>/dev/null &");
            
            push(@splitted_out,".splited.res.phobious$file.out");
        }else {
            print STDERR "unable to find splitted infile $file\n";
        }
    }
    
    #below section is to check if all pfam_scan.pl threads are done.
    threads_finish_examination($userID,@fingerprints);

    #splitted result merge to an intact output $phobius_output
    splitted_results_merge($output,@splitted_out);

    files_remove(@splitted_out);
    files_remove(@renamed_split_files);
}
=cut

sub pfam_parse{
    my $input = shift;
    open(IN,$input) || die "cant open the $input";
    my $evalue = 0.1;
    my $i = 0;
    while (<IN>) {
        $i++;
        chomp;
        next if (/^#/);
        next if (/^\s*$/);
        
        my @array = split/\s+/,$_;
        next if ($array[12] >$evalue);
        my $id = $array[0];
        my $start = $array[1];
        my $end   = $array[2];  
        
        #if ($array[5] =~ /PF00069/i or $array[5] =~ /PF07714/i or $array[5] =~ /PF12398/i) {
        if (pfam_check($array[5])) {
            push(@{$result{$id}->{stk}->{range}},join("|","pfam","$start-$end"));
                   $result{$id}->{stk}->{property}  = 1;
        }
        elsif($array[5] =~ /PF01476\.?/i) {
            push(@{$result{$id}->{LysM}->{range}},join("|","pfam","$start-$end"));
                   $result{$id}->{LysM}->{property} = 1;
        }
        elsif ($_ =~ /LRR/i) {
            push(@{$result{$id}->{LRR}->{range}},join("|","pfam","$start-$end"));
                   $result{$id}->{LRR}->{property} = 1;
        }
    }
    close IN;
}

sub iprscan_parse {
    my $input = shift;
    open(IN,$input) || die "cant open the $input";
    my $evalue = 1e-10;
    while (<IN>) {
        chomp;
        my ($id,$md5,$len,$database,$hitid,$desc,$start,$end,$evalue2,$true,$date,$ipr,$domain_name) = split/\t/,$_;
        #print "$evalue2\n";next;
        #next if ($evalue2 >$evalue);
        #$len,$database,$hitid,($desc)?$desc:'na',$start,$end,$evalue,$true,($ipr)?$ipr:'na',($domain_name)?$domain_name:'na'
        $desc        = ($desc)? $desc : 'na' ;
        $ipr         = ($ipr) ? $ipr : 'na' ;
        $domain_name = ($domain_name)?$domain_name:'na';
        
        #if ($hitid =~ /PF00069/i or $hitid =~ /PF07714/i or $hitid =~ /PF12398/i) {
        if (pfam_check($hitid)) {# if hitid belong to predefined RLK types    
            push(@{$result{$id}->{stk}->{range}},join("|","iprscan","$start-$end"));
                   $result{$id}->{stk}->{property} = 1;
        }
        elsif ($hitid =~ /PF01476/i) {
            push(@{$result{$id}->{LysM}->{range}},join("|","iprscan","$start-$end"));
                   $result{$id}->{LysM}->{property} = 1;
        }
        elsif ($desc =~ /leucine.*rich/i or $domain_name =~ /leucine.*rich/i) {
            push(@{$result{$id}->{LRR}->{range}},join("|","iprscan","$start-$end"));
                   $result{$id}->{LRR}->{property} = 1;
        }
    }
    close IN;
}

sub lst_processing {
    my $input = shift;
    my %lst = ();
    open(IN, $input) || die "cant open the $input file";
    while (<IN>) {
        chomp;
        my @array = split/\s+/,$_;
        #push(@lst,$array[0]);
        $lst{$array[0]} = 1;
    }
    close IN;
    return \%lst;
}

sub Ptime{
          my $time = localtime;
          my ($msg)= @_;
          print STDERR "$time: $msg\n";
}
=for comment
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
    close OUT;
    local $/ = "\n";
    return(@splited_files);
}

sub files_remove{
    my @files = @_;
    foreach my $file (@files) {
        unlink "$file";
    }
}

sub generate_rand_filename
{
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
    
    open(OUT, ">$output");
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
sub file_parser {
    my $file = shift;
    my %pfam = ();
    
    open(IN,"$file") or die "unable to open $file\n";
    while (<IN>) {
        chomp;
        my ($pf,@other) = split/\s+/,$_;
        next unless ($pf =~ /^pf/i);
        $pfam{uc($pf)} = 1;
    }
    close IN;
    return \%pfam;    
}

sub pfam_check {
    my $array = shift;
    
    my $flag = 0;
    foreach my $pfamid (@pfam_rlk) {
        if ( $array =~ /$pfamid/i) {
            $flag++;
        }
    }
    return $flag;    
}