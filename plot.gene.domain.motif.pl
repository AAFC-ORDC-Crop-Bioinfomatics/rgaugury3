#!/usr/bin/perl -w
use strict;
use GD;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use FindBin;

my $version = 0.1;
# -
GetOptions(my $options = {},
              "-gff3=s","-i=s","-m=s"
);

my $USAGE = <<USAGE;
-------------------------------------------------------------------------------------------------------------
Script: to plot the gene domain and motifs
 
-------------------------------------------------------------------------------------------------------------
Version: $version, coded by Pingchuan Li @ Frank You Lab

Notice: please make folder named as motif where the input file located, in which all the png file have been included. 

Arguments:

        -gff3    gene annotation gff3
        -i       motif and domain input
        -m       output meta data

    enjoy it!
USAGE

die $USAGE unless (defined $options->{i});

#gff3 file, specified on command line
my $gff3   = $options->{gff3};

#motif file,  specified on command line
my $input  = $options->{i};
my $output = ($options->{m}) ? $options->{m} : "RGA.gene.structure.meta.txt" ;

# pay attention to the structure of the %gene
# all the gene gff will be imported into the %gene
# gff3 will be parsed by function <gff3_parser>
my %gene   = %{gff3_parser($gff3)};

#save the coordination mapping info for the amino and nucleotide
my %aa_nn_mapping = (); #the value of aa_nn_mapping will be used by many modules

my @no_found_gff3 = ();

#use to zoom in or zoom out; output in filename as x1..
# so far only 1 works, don't use other factor like bigger than 1 or less than 1.
# more job needs to implemented on motif_tile or domain_draw
my $factor = 1;

#to set a fixed length for intro, pixel unit, if specified n for -s
my $fixed_intron_len = 40;

# figure parameters
my $height     = 100;     #pixel
my $left_space = -1;      #pixel
my $top_space  = 2;       #pixel

my %motif_img = (
    1  => $FindBin::Bin.'/motif/Image001.png',
    2  => $FindBin::Bin.'/motif/Image002.png',
    3  => $FindBin::Bin.'/motif/Image003.png',
    4  => $FindBin::Bin.'/motif/Image004.png',
    5  => $FindBin::Bin.'/motif/Image005.png',
    6  => $FindBin::Bin.'/motif/Image006.png',
    7  => $FindBin::Bin.'/motif/Image007.png',
    8  => $FindBin::Bin.'/motif/Image008.png',
    9  => $FindBin::Bin.'/motif/Image009.png',
    10 => $FindBin::Bin.'/motif/Image010.png',
    11 => $FindBin::Bin.'/motif/Image011.png',
    12 => $FindBin::Bin.'/motif/Image012.png',
    13 => $FindBin::Bin.'/motif/Image013.png',
    14 => $FindBin::Bin.'/motif/Image014.png',
    15 => $FindBin::Bin.'/motif/Image015.png',
    16 => $FindBin::Bin.'/motif/Image016.png',
    17 => $FindBin::Bin.'/motif/Image017.png',
    18 => $FindBin::Bin.'/motif/Image018.png',
    19 => $FindBin::Bin.'/motif/Image019.png',
    20 => $FindBin::Bin.'/motif/Image020.png',
    21 => $FindBin::Bin.'/motif/Image021.png',
    22 => $FindBin::Bin.'/motif/Image022.png',
    23 => $FindBin::Bin.'/motif/Image023.png',
    24 => $FindBin::Bin.'/motif/Image024.png',
    25 => $FindBin::Bin.'/motif/Image025.png',
    26 => $FindBin::Bin.'/motif/Image026.png',
    27 => $FindBin::Bin.'/motif/Image027.png',
    28 => $FindBin::Bin.'/motif/Image028.png',
    29 => $FindBin::Bin.'/motif/Image029.png',
    30 => $FindBin::Bin.'/motif/Image030.png'
);

my %hash = (
    1 => 'NBS',
    2 => 'LRR',
    3 => 'TIR',
    4 => 'CC' ,
    5 => 'stk',
    6 => 'tm' ,
    7 => 'sp' ,
    8 => 'LysM'
);

open(META,">$output");

#processing motif input
open(IN, $input);
while (<IN>) {  #each line represent a gene
    chomp;
    my ($id, @array) = split/\s+/,$_; # one id , one tile
        
    if (looks_like_number($gene{$id}->{len})) {
        #keep reading and processing if existing this gene in gff3
    }
    else {
        push(@no_found_gff3,$id);
        next;
    }

    my $genelen =   $gene{$id}->{len};
    my @exon    = @{$gene{$id}->{exon}};
    my @codon   = @{$gene{$id}->{codon}};
    
    # to generate a canvas according to parameters specified on command line
    my $image = new GD::Image(abs($genelen/$factor)         + 1.2*$left_space, $height + 2*$top_space);    

    my $white  = $image->colorAllocate(255,255,255);
    my $gray   = $image->colorAllocate(222,222,222);
    my $gray1  = $image->colorAllocate(177,177,177);
    my $gray2  = $image->colorAllocate(100,100,100);
    my $black  = $image->colorAllocate(0,0,0);       
    my $red    = $image->colorAllocate(255,0,0);      
    my $blue   = $image->colorAllocate(0,0,255);
    
    #save the nn and aa mapping info to the %aa_nn_mapping
    aa_nn_mapping($genelen,$id,@codon);
    
    # ----show gene id ahead of figure, make sure left side has enough space for text-------------
    #$image->stringFT($black,"/home/lipch/jobs_servers/179.10.113.232/jobs/f.plot.domain.pipeline/DejaVuSansMono.ttf",45,0,10,56,"$id");
    
    $image->transparent($white);
    
    # output file name
    my $outputfile = join(".", $id,"png"); #"$id_x$factor.png";
    
    #the linker line between exon and intron
    my $cat_start = 0 ;
    mkdir("img") unless(-d "img");
    open(OUT,">img/$outputfile") or die "cant write to background\n";


        #draw the gene structure in nn base position
        foreach my $exon (@exon) {
            my ($start,$end) = split/\|/,$exon;
            
            #-----------minus strand---------------
            if ($start>$end) {
                $start = abs($start - $genelen) + 1;
                $end   = abs($end   - $genelen) + 1;
            }
            
            my $x1 = int($start/$factor) + $left_space;
            my $y1 = $top_space;
            my $x2 = int($end/$factor)   + $left_space;
            my $y2 = $top_space + $height;
            
            $image->filledRectangle($x1,$y1,$x2,$y2 + 1, $gray);
            if ($cat_start) {
                # if UTR is directly connected with exon, then there would be no linker between utr and codon region,
                # can be applied to other components besides utr and codon region
                if ($x1  - $cat_start == 1) { 
                }
                else {
                    $image->filledRectangle(abs($cat_start) + 1, abs((2*$top_space+$height)/2) - 2, $start+$left_space -1, abs((2*$top_space+$height)/2) + 2, $gray2);
                }
            }
            $cat_start = $x2;
        }
         
        # draw cds
        foreach my $codon (@codon) {
            my ($start,$end) = split/\|/,$codon;
            
            #-----------minus strand---------------
            if ($start>$end) {
                $start = abs($start - $genelen) + 1;
                $end   = abs($end   - $genelen) + 1;
            }
            
            my $x1 = int($start/$factor) + $left_space;
            my $y1 = $top_space;
            my $x2 = int($end/$factor)   + $left_space;
            my $y2 = $top_space + $height;
            
            $image->filledRectangle($x1,$y1,$x2,$y2 + 1, $gray1);
        }
        
        #------------ coordination convert from aa to nn------------
        #my %domain_motif_mapping = ();
        my %tile_data = (); # to store motif regions in dna manner.
        foreach my $tag (@array) {
            if ($tag =~ /motif/) {
                my ($motif_number,$p_start,$p_end) = $tag =~ /motif_(\d+?)\|(\d+?)\-(\d+)/;
                my @coor = mapping_coor($p_start,$p_end,$id);
                
                #print join("\t",$motif_number, @coor);print "\n";
                
                #print "@coor\n";
                push(@{$tile_data{motif}->{$motif_number}},@coor);
            }
        }
        
        # due to redundancy of motif definition, some of them are merged by coordiation.
        my %tile_data_unique = %{unique_motif_meta(\%tile_data)};
        
        #tile up the image for the motif, sort motif is pretty important.
        foreach my $motif (sort {$a cmp $b} keys %tile_data_unique) {
            #print "$motif\t";
            foreach my $type (sort {$a <=> $b} keys %{$tile_data_unique{$motif}}) {
                print META join("\t", $id, $hash{$type},"color_".$type );
                my @array = @{$tile_data_unique{$motif}->{$type}};
                print META "\t";
                foreach my $region (@array) {
                    my ($start,$end) = @{$region};
                    print META join("|",$start,$end);
                    print META "\t";

                    motif_tile($start,$end,$motif_img{$type},$image);
                }
                print META "\n";
            }
        }
    
    print OUT $image->png;
    close OUT;
}
close META;
ouptut_unfound_gff3_list();

# -------------------- function ------------------------
sub motif_tile{
    #problem yet to be resolved in this subfunction, but it works well now. the only potential issue is that if the motif_image is not sorted while invoking 'foreach', then the plot figure will randomly change the color pattern.
    
    my ($start, $end, $motif_image, $im) = @_;
    my $xpos = $start;
    my $ypos = $top_space + 2;
    #print "$motif_image\n";
    while ($xpos <$end) {
        my $width = (($end - $xpos + 1)>=100) ? 100 : $end - $xpos + 1;
        my $motif_image_new = GD::Image->newFromPng("$motif_image",0);
        $im->copy($motif_image_new, $xpos + $left_space, $ypos - 1 , 0, 0, $width, 100); #$top_space + $height - 5
        $xpos += 100;
    }
}

sub domain_draw {
    my ($start,$end,$color,$im) = @_;
    my $x1pos = $start;
    my $x2pos = $end;
    my $y1pos = $top_space + 2;
    my $y2pos = $height;
    
    $im->filledRectangle($x1pos+ $left_space,$y1pos,$x2pos + $left_space,$y2pos,$color);
}

sub gff3_parser {
    #All the data will saved in %gene, especially the @exon
    my ($file) = @_;
    my $offset = 0;
    my $len    = 0;
    
    # use to store original gff3
    my %gff3   = ();
    my %gene   = ();
    
    # use to store sorted gff3 by column 8, 3 and column 4
    my @gff3_sorted = ();
    
    # sorting gff3 by geneid (column 8)
    open(IN,$file) or die "unable to open $file";
    while (<IN>) {
        chomp;
        next if (/^#/ or /^\s*$/);
        my @array = split/\t/,$_;
        
        # only Feature : mRNA, CDS, UTR and gene are acceptable for parsing.
        next unless($array[2]=~ /mRNA/i or $array[2]=~ /exon/i or $array[2]=~ /CDS/i  or $array[2]=~ /UTR/i or $array[2]=~ /gene/i);
        
        my ($geneid) = $array[8] =~ /ID=(.*?)\;/i;
        
        ($array[3],$array[4]) = ($array[4],$array[3]) if ($array[3]>$array[4]);
        
        if ($array[2] =~ /mRNA/i or $array[2] =~ /gene/i) {
            $gff3{$geneid}->{mRNA} = $_;
        }
        elsif ($array[2] =~ /CDS/i or $array[2] =~ /exon/i or $array[2] =~ /UTR/i) {
            push(@{$gff3{$geneid}->{exon}},[@array]);
        }
    }
    close IN;

    #sorting each gene by column 3 and column 4
    foreach my $id (sort {$a cmp $b} keys %gff3) {
        my $header = $gff3{$id}->{mRNA};
        
        my @body = @{$gff3{$id}->{exon}};
        my @sorted = sort {$a->[3] <=> $b->[3] or $a->[4] <=> $b->[4]} @body;
        
        push(@gff3_sorted,[my @tmp = split/\t/,$header]);
        foreach my $key (@sorted) {
            push(@gff3_sorted,$key);
        }
    }
    
    # processing gff3 parsing
    foreach my $ref (@gff3_sorted) {
        my @array = @{$ref};

        if ($array[2] eq 'mRNA') {
            $offset = minone($array[3],$array[4]) - 1;
            my ($geneid)  = $array[8] =~ /ID=(.*?)\;/i;
            
            $gene{$geneid}->{len} = abs($array[4] - $array[3]) + 1;
        }
        elsif ($array[2] =~ /CDS/i  or $array[2] =~ /exon/i or $array[2] =~ /UTR/i) {
            my $start = minone($array[3],$array[4]) - $offset;
            my $end   = maxone($array[3],$array[4]) - $offset;

            my ($geneid) = $array[8] =~ /ID=(.*?)\;/i;

            # ---------consider the strand info-----------
            if ($array[6] eq '+') {
                push(@{$gene{$geneid}->{exon}},  join("|",$start, $end));
                push(@{$gene{$geneid}->{codon}}, join("|",$start, $end)) if ($array[2] =~ /CDS/i);
            }
            else {
                #because flax gff3 always put smaller number at the first, thus in minus strand the smaller coordination is actually the 3' end
                unshift(@{$gene{$geneid}->{exon}}, join("|", $end, $start));
                unshift(@{$gene{$geneid}->{codon}},join("|", $end, $start))  if ($array[2] =~ /CDS/i); 
            }
        }
        else {
            next;
        }        
    }
    
    close IN;
    return (\%gene);
}

sub aa_nn_mapping{
    # this function will generate a mapping table for each aa<=>nn-represented by pixel position
    my ($genelen,$geneid,@codon) = @_;
    
    # put into all DNA exon coordination into @tmp array
    my @tmp = ();
    
    #to record the position coordination prior to the start of coding bp
    #my $coor = 0;
    
    foreach my $cds(@codon){
        my ($start,$end) = split/\|/, $cds;

        if ($start>$end) {#when gene located in minus strand
            $start = abs($start - $genelen) + 1;
            $end   = abs($end   - $genelen) + 1;
        }
        
        for my $j ($start..$end) {
           push(@tmp,$j);
        }
    }
    my $ab = ($#tmp + 1)%3;
    #die "cds has problem" unless ($ab == 0);

    my $aa_number = 0;

    my @tmp3 = ();

    foreach my $number (sort {$a <=> $b}  @tmp){
        
        push(@tmp3, $number);
        
        if ($#tmp3 == 2) {  #until temporary tmp3 array has 3 scalor.

            $aa_number++;
            
            # the mapping info is aa => (nn,nn,nn);
            push(@{$aa_nn_mapping{$geneid}->{$aa_number}},@tmp3);
            
            #print "@tmp3\n";
            @tmp3 = ();
        }
    }
    
    # --------------debug test purpose ------------------
    #foreach my $gene (keys %aa_nn_mapping) {
    #    foreach my $aa_p (sort {$a <=> $b} keys %{$aa_nn_mapping{$gene}}) {
    #        my @array = @{$aa_nn_mapping{$gene}->{$aa_p}};
    #        #print "$aa_p => (@array)\n";
    #    }
    #}
}

#to conver the motif aa coor to the nn coor
sub mapping_coor {
   my($p_start, $p_end, $geneid) = @_;
   my @tmp = ();
   for my $i ($p_start..$p_end) {
        push(@tmp,@{$aa_nn_mapping{$geneid}->{$i}});
   }

   my @result;
   
    my $index1 = 0;
    my $index2 = 0;
    
    foreach my $i (sort {$a <=> $b} @tmp ) {
        
        $index1++;
        
        if ($index1 == 1) {
            $result[$index2]->[0] = $i;
            $result[$index2]->[1] = $i;
        }
        else {
            if ($i - $result[$index2]->[1] <= 1) {
                $result[$index2]->[1] = $i;
            }
            else {
                $index2++;
                $result[$index2]->[0] = $i;
                $result[$index2]->[1] = $i;
            }
        }
    }

    my @regions = ();
    foreach my $i (@result) {
        
        my $start = $i->[0];
        my $end   = $i->[1];
        
        #each of exon block will be push to regions array.
        push(@regions,join("|",$start,$end));
    }
    return @regions;
}

sub ouptut_unfound_gff3_list {
    foreach my $id (sort {$a cmp $b} @no_found_gff3) {
        print "$id wasn\'t found in the gff3 file\n";
    }
}

sub minone {
    my @array = @_;
    my @a = sort(@array);
    return $a[0];
}
sub maxone {
    my @array = @_;
    my @a = sort(@array);
    return $a[-1];
}

sub unique_motif_meta {
    my $motif_regions_ref = shift;
    my %motif_regions = %{$motif_regions_ref};
    my %motif_regions_unique = ();
    
    foreach my $class (sort {$a cmp $b} keys %motif_regions) {
        foreach my $type (sort {$a <=> $b} keys %{$motif_regions{$class}}) {
            my @regions1 = @{$motif_regions{$class}->{$type}};
            my @regions2 = ();
            foreach my $region (@regions1) {
                #print join("\t",$class,$type,$region,"\n");
                my ($start,$end) = split/\|/,$region;
                ($start,$end) = ($end,$start) if ($start>$end);
                push(@regions2, [$start, $end]);
            }
                
            my @sorted = sort {$a->[0] <=> $b->[0]} @regions2;
            my @merged = ();
            
            my $tmp_start = "";
            my $tmp_end   = "";
            my $i = 0;
            
            foreach my $ref (@sorted) {
                $i++;
                my ($start,$end) = @{$ref};
                
                if ($i == 1) {
                    $tmp_start = $start;
                    $tmp_end   = $end  ;
                }
                else {
                    if ($start>$tmp_end) {
                        push(@{$motif_regions_unique{$class}->{$type}},[$tmp_start, $tmp_end]);
                        $tmp_start = $start;
                        $tmp_end   = $end  ;
                    }
                    elsif ($start>=$tmp_start and $start<=$tmp_end) {
                        if ($end >= $tmp_end) {
                            $tmp_end  = $end;
                        }
                        elsif ($end<$tmp_end) {

                        }
                    }
                    elsif ($start<$tmp_start) {
                        print STDERR "meta data wasn't successfully sorted prior to merge\n";
                    }
                }
            }
            push(@{$motif_regions_unique{$class}->{$type}},[$tmp_start, $tmp_end]);
        }
    }
    return(\%motif_regions_unique);
}

