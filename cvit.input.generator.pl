#!/usr/bin/perl -w
use strict;
use Scalar::Util qw(looks_like_number);
use Getopt::Long;

#this script is better than extract.gff.to.CVIT.kits.by.gene.id.pl

GetOptions(my $options = {},
              "-l=s","-f=s","-t=s","-c=s","-t2=s","-p=s","-dir=s","-pfx=s"
);


my $USAGE = <<USAGE;
-------------------------------------------------------------------------------------------------------------
Script: to extract gff3 for CVIT kits by gene id, specifically for flax, minor modification for other genome

Usage:   perl code.pl -l lst -f gene_gff3 -t gene
-------------------------------------------------------------------------------------------------------------
Coded by Pingchuan Li
Arguments:

        -l       gene list file, first column is the gene id
        -p       protein fasta file name, primary annnotation is recommended
        -f       gff3 file
        -t       gff3 types in 3rd column, like 'gene','mRNA','sRNA', default = 'gene'
        -c       color of annoated gene,'red','blue','green','yellow','orange','black'
        -t2      gene type,like 'RLK', "RLP", "NBS"
        -dir     output directory
        -pfx     prefix for output filename

enjoy it!
USAGE

die $USAGE unless (defined $options->{l});

my %list = ();
my %gene = ();
#my $prefix = $options->{cp};
my %id = ();#all id derived from protein fasta file
my $type = ($options->{t}) ? $options->{t} : 'gene';
my $type2 = $options->{t2};
my $dir = ($options->{dir}) ? $options->{dir}: "./" ;
my $prefix = ($options->{pfx}) ? $options->{pfx} : "input" ;  $prefix =~ s/\.$//g;

local $/ = ">";

open(IN,$options->{p});
while (<IN>) {
    chomp;
    my ($title,$seq) = split/\n/,$_,2;
    next unless ($title and $seq);
    
    my ($id) = $title =~ /([a-zA-Z0-9\.\-\_]+)/;#permit a-z and 0-9 plus dot,sometimes a gene name has dot, ex. Traes_1AL_002FAE6E8.1
    $id{$id} = 1;
}
close IN;
local $/ = "\n";

open(IN,"$options->{l}");
while (<IN>) {
    chomp;
    my @array = split/\t/,$_;
    $list{$array[0]} = 1;
}
close IN;

open(IN,"$options->{f}");
while (<IN>) {
    chomp;
    next unless (/^\w/);
    my @array = split/\t/,$_;
    
    if ($array[2] =~ /mrna/i or $array[2] =~ /transcript/i) {
        # modifiy the below phrase to adapt to other genome gff3
        #my ($id)  = $array[8] =~ /NAME=(.*?)\./;
        
        #below methods is slower than above method
        my $id = '';
        foreach my $key (keys %list) {
            if ($array[8] =~ /$key/) {
                $id = $key;
            }
        }
        
        if ($id) {
            #modify the below contion to adapt to other genome gff3
            my $chr = $array[0];
            my ($start,$end)    = ($array[3],$array[4]);
            $gene{$id}->{chr}   = $chr;
            $gene{$id}->{start} = $start;
            $gene{$id}->{end}   = $end;
            $gene{$id}->{str}   = $array[6];
            $gene{$id}->{color} = $options->{c};
        }
        else {
            #Ptime("No hits found in LIST to match the gene ID within $array[8]") unless ($id);
            next;
        }
    }    
}
close IN;

open(OUT,">$dir/$prefix.CViT.$type2.$options->{c}.txt"); 
foreach my $id (sort {$a cmp $b} keys %gene) {
                    my $chr   =  $gene{$id}->{chr};
                    my $start =  $gene{$id}->{start};
                    my $end   =  $gene{$id}->{end};
                    my $str   =  $gene{$id}->{str};
                    my $color =  $gene{$id}->{color};

                    print OUT join("\t",$chr, 'cvit', $type, $start, $end, ".",$str,".","Name=;color=$color");
                    print OUT "\n";
}
close OUT;

# --------------sub---------------------
sub Ptime{

          my $time = localtime;
          my ($msg)= @_;
          print STDERR "$time: $msg\n";

}























