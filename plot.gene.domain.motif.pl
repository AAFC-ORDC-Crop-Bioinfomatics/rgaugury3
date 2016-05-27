#!/usr/bin/perl -w
use strict;
use GD;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use FindBin;

my $version = 0.1;
# -
GetOptions(my $options = {},
              "-gff3=s","-i=s","-s=s"
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
        -s       plot the intron lenght in ratio to extron[y] or no[n], default = 'n'
    enjoy it!
USAGE

die $USAGE unless (defined $options->{i});

#gff3 file
my $gff3   = $options->{gff3};

#motif file
my $input  = $options->{i};

#intron scale true or false
my $scale  = ($options->{s})?($options->{s}):'n';

#pay attention to the structure of the %gene
my %gene   = ();
my %gene_alt = ();

#save the coordination mapping info for the amino and nucleotide
my %aa_nn_mapping = (); #the value of aa_nn_mapping will be used by many modules

my @no_found_gff3 = ();

#use to zoom in or zoom out; output in filename as x1..
my $factor = 1;

#to set a fixed length for intro, pixel
my $fixed_intron_len = 50;

my $height = 100;
my $left_space = 20; #pixel
my $top_space  = 2;   #pixel

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

# all the gene gff will be imported into the %gene

# below function will populate the gene hash with standard method
gff3_parser($gff3); 

#processing motif input
open(IN, $input);
while (<IN>) {
    chomp;
    my ($id,@array) = split/\s+/,$_; # one id , one tile
    
    if (looks_like_number($gene{$id}->{len})) {
        #keep reading and processing if existing this gene in gff3
    }
    else {
        push(@no_found_gff3,$id);
        #print "id is $id\n";
        next;
    }

    my $genelen = $gene{$id}->{len};
    my @cds = @{$gene{$id}->{cds}};
    
    #initialize the modifed gene length, because the intron length has been modifed as fixed legnth for each of them
    my $epi_gene_length = 0;
    my $intron_No = -1;
    
    foreach my $cds (@cds) {
        my ($start,$end) = split/\|/,$cds;
        $epi_gene_length += abs($start - $end + 1);
        $intron_No++;
    }
    
    my $image = '';
    if ($scale eq 'y') {
        $image = new GD::Image(abs($genelen/$factor)+1.2*$left_space, $height+2*$top_space);    
    }
    elsif ($scale eq 'n') {
        #the final length for a modified gene
        $epi_gene_length += $intron_No*$fixed_intron_len;
        
        $image = new GD::Image(abs($epi_gene_length/$factor)+1.2*$left_space, $height+2*$top_space);    
    }
    else {
        die "scale can be only y or n\n";
    }

    
    #save the nn and aa mapping info to the %aa_nn_mapping
    aa_nn_mapping($genelen,$epi_gene_length,$id,$scale,@cds);
    
    my $white = $image->colorAllocate(255,255,255);
    my $gray  = $image->colorAllocate(222,222,222);
    my $gray1  = $image->colorAllocate(155,155,155);
    my $black = $image->colorAllocate(0,0,0);       
    my $red   = $image->colorAllocate(255,0,0);      
    my $blue  = $image->colorAllocate(0,0,255);
    
    #$image->stringFT($black,"/home/lipch/jobs_servers/179.10.113.232/jobs/f.plot.domain.pipeline/DejaVuSansMono.ttf",45,0,10,56,"$id");
    
    my  %domain_color = (
        'nbs' => $image->colorAllocate(155,155,155),
        'NBS' => $image->colorAllocate(155,155,155),
        #'lrr' => $image->colorAllocate(55,55,55),
        #'LRR' => $image->colorAllocate(55,55,55),
        'tir' => $image->colorAllocate(90,90,90),
        'TIR' => $image->colorAllocate(90,90,90),
        'cc'  => $image->colorAllocate(120,120,120),
        'CC'  => $image->colorAllocate(120,120,120)
    );
    
    $image->transparent($white);
    
    my $outputfile = join(".", $id,"png"); #"$id_x$factor.png";
    
    #the linker line between exon and intron
    my $cat_start = "";
    mkdir("img") unless(-d "img");
    open(OUT,">img/$outputfile") or die "cant write to background\n";
    

    if ($scale eq 'n') {
        # -------------------draw the gene $factorcuture, including exon and intron(linker)-----------------------
        my $coor = 0;
        my $flag = -1;
        foreach my $cds (@cds) {
            $flag++; # to get the number of intron
            my ($start,$end) = split/\|/,$cds;
            
            my $exon_len = int(abs($start - $end)/$factor) + 1;
            
            my $x1 = $coor + $left_space + 1;  #pixel coordination
            my $y1 = $top_space;
            
            my $x2 = $x1 + $exon_len - 1;   #pixel coordination
            my $y2 = $top_space + $height;
                
            $coor += $exon_len + $fixed_intron_len;  #pixel coordination
            
            #-------------draw exon structure-----------------
            $image->filledRectangle($x1,$y1,$x2,$y2, $gray);
            
            #-------------draw intron linker------------------
            if ($cat_start) {
                my $cat_end = $cat_start + $fixed_intron_len + 1;
                $image->filledRectangle($cat_start + 1 ,abs((2*$top_space+$height)/2) - 2, $cat_end - 1, abs((2*$top_space+$height)/2) + 2,$gray1);
            }
            $cat_start = $x2;
        }
        
        #------------ coordination convert from aa to nn------------
        #my %domain_motif_mapping = ();
        my %tile_data = ();
        foreach my $tag (@array) {
            if ($tag =~ /motif/) {
                my ($motif_number,$p_start,$p_end) = $tag =~ /motif_(\d+?)\|(\d+?)\-(\d+)/;
                my @coor = mapping_coor($p_start,$p_end,$id);
                #print "@coor\n";
                push(@{$tile_data{motif}->{$motif_number}},@coor);
            }
            elsif ($tag =~ /domain/) {
                my ($domain_name,$p_start,$p_end) = $tag =~ /domain_(\S+?)\|(\d+?)\-(\d+)/;
                my @coor = mapping_coor($p_start,$p_end,$id);
                #print "@coor\n";
                push(@{$tile_data{domain}->{$domain_name}},@coor);
            }
        }
        
        #plot the domain
        foreach my $name (sort {$a cmp $b} keys %{$tile_data{domain}}) {
            my @array = @{$tile_data{domain}->{$name}};
            #print join("\t",@array,"\n");
            foreach my $cat_coor (@array) {
                my ($start,$end) = split/\|/,$cat_coor;
                domain_draw($start,$end,$domain_color{$name},$image);
            }
        }
        
        #tile up the image for the motif
        foreach my $motif(sort {$a <=> $b} keys %{$tile_data{motif}}) {
            my @array = @{$tile_data{motif}->{$motif}};
            foreach my $cat_coor (@array) {
                my ($start,$end) = split/\|/,$cat_coor;
                motif_tile($start,$end,$motif_img{$motif},$image);
            }
        }                
    }
    elsif ($scale eq 'y') {
        foreach my $cds (@cds) {
            my ($start,$end) = split/\|/,$cds;
            
            #-----------minus strand---------------
            if ($start>$end) {
                $start = abs($start - $genelen) + 1;
                $end   = abs($end - $genelen) + 1;
            }
            
            my $x1 = int($start/$factor) + $left_space;
            my $y1 = $top_space;
            my $x2 = int($end/$factor) + $left_space;
            my $y2 = $top_space + $height;
            
            $image->filledRectangle($x1,$y1,$x2,$y2, $gray);
            if ($cat_start) {
                $image->filledRectangle(abs($cat_start) + 1, abs((2*$top_space+$height)/2) - 2, $start+$left_space -1, abs((2*$top_space+$height)/2) + 2,$gray1);
            }
            $cat_start = $x2;
        }
        
        #------------ coordination convert from aa to nn------------
        #my %domain_motif_mapping = ();
        my %tile_data = ();
        foreach my $tag (@array) {
            #print "$tag\n";
            if ($tag =~ /motif/) {
                my ($motif_number,$p_start,$p_end) = $tag =~ /motif_(\d+?)\|(\d+?)\-(\d+)/;
                my @coor = mapping_coor($p_start,$p_end,$id);
                #print "@coor\n";
                push(@{$tile_data{motif}->{$motif_number}},@coor);
            }
            elsif ($tag =~ /domain/) {
                
                my ($domain_name,$p_start,$p_end) = $tag =~ /domain_(\S+?)\|(\d+?)\-(\d+)/;
                my @coor = mapping_coor($p_start,$p_end,$id);
                #print "@coor\n";
                push(@{$tile_data{domain}->{$domain_name}},@coor);
            }
        }
        
        #plot the domain
        foreach my $name (keys %{$tile_data{domain}}) {
            my @array = @{$tile_data{domain}->{$name}};
            #print join("\t",@array,"\n");
            foreach my $cat_coor (@array) {
                my ($start,$end) = split/\|/,$cat_coor;
                domain_draw($start,$end,$domain_color{$name},$image);
            }
        }
        
        #tile up the image for the motif
        foreach my $motif(sort {$a <=> $b} keys %{$tile_data{motif}}) {
            my @array = @{$tile_data{motif}->{$motif}};
            foreach my $cat_coor (@array) {
                my ($start,$end) = split/\|/,$cat_coor;
                motif_tile($start,$end,$motif_img{$motif},$image);
            }
        }                
    }
    else {
        die "scale can be only y or n\n";
    }
    
    print OUT $image->png;
    close OUT;
}

ouptut_unfound_gff3_list();

# -------------------- function ------------------------
sub motif_tile{
    #problem yet to be resolved in this subfunction, but it works well know. the only potential issue is that if the motif_image is not sorted for its invoking 'foreach', then the plot figure will randomly change the color pattern.
    
    my ($start,$end,$motif_image,$im) = @_;
    my $xpos = $start;
    my $ypos = $top_space + 2;
    #print "$motif_image\n";
    while ($xpos <$end) {
        my $width = (($end - $xpos + 1)>=100) ? 100 : $end - $xpos + 1;
        my $motif_image_new = GD::Image->newFromPng("$motif_image",0);
        $im->copy($motif_image_new,$xpos + $left_space,$ypos - 1 , 0, 0, $width, 100); #$top_space + $height - 5
        
        $xpos += 100;
        last if ($xpos >=100);
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
    #All the data will saved in %gene, especially the @cds
    my ($file) = @_;
    my $offset = 0;
    
    my $len = 0;
    open(IN,$file);
    while (<IN>) {
        chomp;
        
        my @array = split/\t/,$_;
        if ($array[2] eq 'gene') {
            $len = abs($array[4] - $array[3]) + 1;
            $offset = $array[3] - 1;
        }
        elsif ($array[2] eq 'mRNA') {
            my $geneid  = $array[8];
            
            $gene{$geneid}->{len} = $len;
        }
        elsif ($array[2] =~ /CDS/i or $array[2] =~ /UTR/i) {
            my $start = $array[3] - $offset;
            my $end   = $array[4] - $offset;

            my $geneid = $array[8];

            # consider the strand info-----------
            if ($array[6] eq '+') {
                push(@{$gene{$geneid}->{cds}},join("|",$start,$end));
            }
            else {
                #because flax gff3 always put smaller number at the first, thus in minus strand the smaller coordination is actually the 3' end
                unshift(@{$gene{$geneid}->{cds}},join("|",$end,$start));
            }
        }
        else {
            next;
        }
    }
    close IN;
}

sub aa_nn_mapping{
    my ($genelen,$epi_gene_len,$geneid,$scale_y_n,@cds) = @_;
    my $i = 0;
    
    # put into all DNA exon coordination into @tmp array
    my @tmp = ();
    
    #to record the coordination of last intron
    my $coor = 0;
    
    foreach my $exon(@cds){
        my ($start,$end) = split/\|/,$exon;
        
        if ($start>$end) {
            $start = abs($start - $genelen) + 1;
            $end   = abs($end - $genelen) + 1;
        }
        
        if ($scale_y_n eq 'n') {
            
            my $start_new = $coor + 1;
            my $end_new   = $start_new + abs($start - $end);
            
            $start = $start_new;
            $end   = $end_new;
            $coor += (abs($start_new - $end_new) + 1) + $fixed_intron_len;
        }
        
        #print join("\t",$geneid,$start,$end);
        #print "\n";
        
        
        for my $j ($start..$end) {
           push(@tmp,$j);
        }
    }
    my $ab = ($#tmp + 1)%3;
    #die "cds has problem" unless ($ab == 0);
    
    my $nn_number = 0;
    my $aa_number = 0;
    
    my @tmp3 = ();
    
    
    foreach my $number (@tmp){
        $nn_number++;
        if ($nn_number%3 == 0) {
            $aa_number++;

            push(@tmp3,$number);
            # the mapping info is aa => (nn,nn,nn);
            push(@{$aa_nn_mapping{$geneid}->{$aa_number}},@tmp3);
            
            #print "@tmp3\n";
            @tmp3 = ();
        }
        else {
            push(@tmp3, $number);
        }
    }
}

#to conver the motif aa coor to the nn coor
sub mapping_coor {
   my($p_start,$p_end,$geneid) = @_;
   my @tmp = ();
   for my $i ($p_start..$p_end) {
        push(@tmp,@{$aa_nn_mapping{$geneid}->{$i}});
   }

   my $flag = 0;
   my $tmp;
   
   my $start;
   my @result;
   my $total = $#tmp + 1;
   
   foreach my $coor (sort {$a <=> $b} @tmp) {
        $flag++;
        if ($flag == 1) {
            $tmp = $coor;
            $start = $coor;
        }
        elsif ($flag == $total){
            if ($coor - $tmp == 1) {
                push(@result,join("|",$start,$coor));
            }
            
            else{
                push(@result,join("|",$start,$tmp));
                push(@result,join("|",$coor,$coor));
            }
        }
        else {
            if ($coor - $tmp == 1)   {#if motif was broken by intron, this sub will help to tell
                $tmp = $coor;
            }
            else {
                push(@result,join("|",$start,$tmp));
                $start = $coor;
                $tmp = $coor;
            }
        }
   }
   return @result;
}

sub ouptut_unfound_gff3_list {
    foreach my $id (sort {$a cmp $b} @no_found_gff3) {
        print "$id wasn\'t found in the gff3 file\n";
    }
}