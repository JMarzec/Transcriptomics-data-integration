#!/usr/bin/perl -w
################################################################################
#
#   File name: ArrayAnnot.pl
#
#   Authors: Jacek Marzec ( j.marzec@qmul.ac.uk )
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

################################################################################
#
#	Description: Script retrieving Ensembl microarray probe mapping results for selected arrays. Script adopted from http://cvs.sanger.ac.uk/cgi-bin/viewvc.cgi/ensembl-functgenomics/scripts/examples/microarray_annotation_example.pl?revision=1.1&root=ensembl&view=markup . It requiries Bioperl, Perl Application Programme Interface (API), 'ensembl' and 'ensembl-functgenomics' packages to be installed a priori.
#
#	Command line use example: ./ArrayAnnot.pl -i /scratch/jack/data/PhD/Transcriptomics_project/GenExpression_InputFiles.txt
#
#	-i info_file:       File name with full path that lists the datasets to be analysed. The file is expected to include the following columns: (1) dataset name; (2) target file name with full path; (3) data type (raw or processed) and (4) platform used
#
################################################################################

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Data::Dumper;

sub usage ();

#===============================================================================
#    Functions
#===============================================================================

#===============================================================================
#    Main
#===============================================================================

our @ARGV;
my $arg;
my $inFile;
my @dir;
my $dir;
my %Arrays = ();
my $annotFile;
my $registry;
my $slice_adaptor;
my $pfa;
my $aa;

while ($arg = shift) {
    
    if ($arg =~ /^-i$/) {
        $inFile = shift;
    } elsif ($arg =~ /help/) {
        usage ();
    }
}

if ($inFile) {
    @dir = split('/', $inFile);
    $dir = join('/', @dir[0 .. $#dir-1]);
} else {
    usage ();
}

##### Get information about datasets to be analysed
open (INFILE, $inFile) or die $!;

while (my $record = <INFILE>) {

    chomp $record;

    my @info = split(/\t/, $record);

    ##### Create hash with arrays to be analysed
    $Arrays{ $info[3] } = $info[3];
}
close(INFILE);

##### Use the Registry, load database details and then get the slice/array adaptors
$registry = "Bio::EnsEMBL::Registry";

$registry->load_registry_from_db(
				-host => 'ensembldb.ensembl.org',
				-user => 'anonymous',
				-db_version => 87,
				 );

print "Using api: ",$registry->software_version,"\n";

$slice_adaptor = $registry->get_adaptor( "human", "core", "Slice" );
$pfa           = $registry->get_adaptor( "human", "funcgen", "ProbeFeature" );
$aa            = $registry->get_adaptor( "human", "funcgen", "Array" );

################################################################################
#
#			AFFYMETRIX HuEx-1_0-st-v2                                          
#
################################################################################

if ( exists  $Arrays{ 'Affy_HuEx1ST' } ) {

    $annotFile = $dir."/Affy_HuEx1ST_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HuEx-1_0-st-v2', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
    
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			AFFYMETRIX HuGene-2_0-st-v1
#
################################################################################

if ( exists  $Arrays{ 'Affy_HuGene2ST' } ) {
    
    $annotFile = $dir."/Affy_HuGene2ST_annot.txt";
    
    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );
    
    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
        
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
        
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HuGene-2_0-st-v1', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";
        
        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
            
            $probename =~ s/;//g;
            
            my @probe_region = split(  ":", $slice_name );
            
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
            
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
            
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
            
            foreach my $gene  (@genes) {
                $gene_no++;
                
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
            
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			AFFYMETRIX HuGene-1_0-st-v1                                          
#
################################################################################

if ( exists  $Arrays{ 'Affy_HuGene1ST' } ) {

    $annotFile = $dir."/Affy_HuGene1ST_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HuGene-1_0-st-v1', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
    
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			AFFYMETRIX HG-U133_Plus_2                                          
#
################################################################################

if ( exists  $Arrays{ 'Affy_U133Plus2' } ) {

    $annotFile = $dir."/Affy_U133Plus2_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HG-U133_Plus_2', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
    
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			AFFYMETRIX HG-U133A_2                                          
#
################################################################################

if ( exists  $Arrays{ 'Affy_U133A2' } ) {

    $annotFile = $dir."/Affy_U133A2_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HG-U133A_2', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
    
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			AFFYMETRIX HG-U133A                                          
#
################################################################################

if ( exists  $Arrays{ 'Affy_U133A' } ) {

    $annotFile = $dir."/Affy_U133A_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HG-U133A', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
    
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			AFFYMETRIX HG-U133B                                          
#
################################################################################

if ( exists  $Arrays{ 'Affy_U133B' } ) {

    $annotFile = $dir."/Affy_U133B_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HG-U133B', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
    
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			AFFYMETRIX HG-U95Av2                                          
#
################################################################################

if ( exists  $Arrays{ 'Affy_U95Av2' } ) {

    $annotFile = $dir."/Affy_U95Av2_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HG-U95Av2', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
    
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			AFFYMETRIX HG-U95B                                          
#
################################################################################

if ( exists  $Arrays{ 'Affy_U95B' } ) {

    $annotFile = $dir."/Affy_U95B_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HG-U95B', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
    
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			AFFYMETRIX HG-U95C                                          
#
################################################################################

if ( exists  $Arrays{ 'Affy_U95C' } ) {

    $annotFile = $dir."/Affy_U95C_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" ) or die $!;
    print ( ANNOT_FILE "probe\tcell\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor( 'HG-U95C', 'AFFY' );
        
        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probeset   = $pf->probeset->name;
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probeset\t$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
    
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#			ILLUMINA HumanHT-12
#
################################################################################

if ( exists  $Arrays{ 'Illum_HT_12_V3' } ) {

    $annotFile = $dir."/Illum_HT_12_V3_annot.txt";

    open ( ANNOT_FILE, ">$annotFile" );
    print ( ANNOT_FILE "probe\tlocus\tgene_id\tgene_short_name\tcigar\tgene_no\n" );

    ##### Get the slice for each chromosome and get oligo array object
    for (my $i=1; $i<25; $i++) {
    
        my $chr = $i;
        $chr =~ s/23/X/;
        $chr =~ s/24/Y/;
    
        my $slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr );
        my $array = $aa->fetch_by_name_vendor('HumanHT-12', 'ILLUMINA');

        ##### Get mapped probes for that slice
        my $oligo_features = $pfa->fetch_all_by_Slice_Array($slice,$array);
        print "\nProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n";

        foreach my $pf(@$oligo_features){
            my $probename  = $pf->probe->get_probename($array->name);
            my $slice_name = $pf->feature_Slice->name;
    
            $probename =~ s/;//g;
    
            my @probe_region = split(  ":", $slice_name );
    
            ##### Report probe name and genomic position
            print ( ANNOT_FILE "$probename\t"."chr".$probe_region[2].":".$probe_region[3]."-".$probe_region[4]."\t" );
            
            ##### Get the  ensembl gene ID and short name
            my $probe_slice = $slice_adaptor->fetch_by_region( 'chromosome', $probe_region[2], $probe_region[3], $probe_region[4] );
            my @genes = @{ $probe_slice->get_all_Genes() };
    
            ##### Consider only the first gene record (probably the more certain). Extended cigar_string as defined by SAMTools group
            my $gene_no = 0;
    
            foreach my $gene  (@genes) {
                $gene_no++;
        
                if ( $gene_no == 1 ) {
                    print( ANNOT_FILE $gene->display_id()."\t".$gene->external_name."\t".$pf->cigar_string."\t" );
                }
            }
        
            if ( $gene_no == 0 ) {
                print( ANNOT_FILE "-\t-\t-\t".$gene_no."\n" )
            } else {
                print( ANNOT_FILE $gene_no."\n" )
            }
        }
    }
}
close(ANNOT_FILE);

################################################################################
#
#		Export all microarray platforms and associated information
#       ( http://www.ensembl.org/info/docs/api/funcgen/regulation_tutorial.html )
#
#
##### Grab the adaptors
#my $array_adaptor = $registry->get_adaptor('Human','funcgen','array');
#
##### Grab all the arrays
#my @array = @{$array_adaptor->fetch_all};
#
##### Print some array info
#foreach my $array ( @array ){
#    print "\nArray:\t".$array->name."\n";
#    print "Type:\t".$array->type."\n";
#    print "Vendor:\t".$array->vendor."\n";
#    ##### Grab the ArrayChips from the array design
#    my @array_chips   = @{$array->get_ArrayChips};
#
#    ##### Print some ArrayChip info
#    foreach my $ac ( @array_chips ){
#        print "ArrayChip:".$ac->name."\tDesignID:".$ac->design_id."\n";
#    }
#}
################################################################################

#===============================================================================
#    End of main
#===============================================================================

exit;

#===============================================================================
#    Subroutines
#===============================================================================

sub usage () {
    print <<"EOS";
    
    usage: $0 -i info_file
    
    Input data
    -i info_file:  File name with full path that lists the datasets to be analysed.
                   The file is expected to include the following columns:
                   (1) dataset name; (2) target file name with full path; (3) data
                   type (raw or processed) and (4) platform used
   
EOS
exit;
}
