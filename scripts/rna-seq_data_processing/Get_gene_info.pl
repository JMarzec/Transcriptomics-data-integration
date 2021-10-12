#!/usr/bin/perl -w
################################################################################
#
#   File name: Get_gene_info.pl
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
#	Description: Script retrieving gene symbol, lenght and GC content using Ensembl Perl Application Programme Interface (API). GFF/GTF file format is expected as an  input. It requiries Bioperl, Perl API and 'ensembl' package to be installed a priori
#
#	Command line use example: ./Get_gene_info.pl -i /scratch/jack/genome_annotation/ensembl_Homo_sapiens.GRCh37.74.gtf
#
#	-i gene_ref:        The name of the reference annotation GFF/GTF file with full path
#
################################################################################

use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Number::Range;

sub usage ();

#===============================================================================
#    Functions
#===============================================================================

#===============================================================================
#    Main
#===============================================================================

my $debug = 0;
our @ARGV;
my $arg;
my $geneFile;
my $geneInfoFile;
my $registry;
my $gene_adaptor;
my $slice_adaptor;
my $range = Number::Range->new("1..22");
my $geneSymbol;
my %geneList = ();

while ($arg = shift) {
    
    if ($arg =~ /^-i$/) {
        $geneFile = shift;
    } elsif ($arg =~ /help/) {
        usage ();
    }
}

if ( $geneFile ) {
    
    print("\nProcess genes in $geneFile...\n");
    
} else {
    usage ();
}

##### Get information from Ensembl about all genes !!!very time consuming step depending on number of genes to be annotated!!!
print("\nRetrieving infrmation from Ensembl about captured genes ...\n");

##### Use the Registry, load database details and then get the slice/array adaptors
$registry = "Bio::EnsEMBL::Registry";

$registry->load_registry_from_db(
-host => 'ensembldb.ensembl.org',
-user => 'anonymous',
-db_version => 83,
);

print "\nUsing api: ",$registry->software_version,"\n";

$gene_adaptor = $registry->get_adaptor( "human", "core", "gene" );
$slice_adaptor = $registry->get_adaptor( "human", "core", "slice" );

##### Get genes lenght and GC content
open (GENE_FILE, $geneFile) or die $!;

$geneInfoFile = $geneFile.".gene_info.txt";
open (GENE_INFO, ">$geneInfoFile");

print ( GENE_INFO "ensembl_id\tgene_symbol\tlength\tGC_content\n" );

my $count = 0;

while (my $record = <GENE_FILE>) {

    chomp $record;
    
    if ( $record =~ m/gene_id \"(.*?)\"/ ) {
        
        my $ensemblID = $1;
        
        if ( !exists $geneList{ $ensemblID } ) {
            
            $count++;
            
            if ($count % 1000 == 0) {
                print "$count...\n";
            }
            
            my $gene = $gene_adaptor->fetch_by_stable_id($ensemblID);
        
            ##### Get gene symbol
            if ( $record =~ m/gene_name \"(.*?)\"/ ) {
            
                $geneSymbol = $1;
            } else {
                $geneSymbol = "";
            }
        
            if ( $gene ) {
            
                ##### Store all captured genes in a hash
                $geneList{ $ensemblID } = $geneSymbol;
            
                ##### Get canonical transcript assocaited with a gene ( http://lists.ensembl.org/pipermail/dev/2012-April/007435.html )
                my $canon_tran = $gene->canonical_transcript();
        
                my $tranID = $canon_tran->stable_id;
                my $tranLength = $canon_tran->length;
        
                ##### Get transcript GC content
                my $chr = $canon_tran->slice->seq_region_name;
                my $start = $canon_tran->start;
                my $end = $canon_tran->end;
            
                if ( $range->inrange($chr) or $chr eq "X" or $chr eq "Y" ) {
                
                    my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr, $start, $end);
                    my $GC = $slice->get_base_count->{'%gc'};
                    $GC = $GC/100;
                
                    ##### Write the Ensembl gene IDs along with read-counts into a file
                    print ( GENE_INFO $ensemblID."\t".$geneSymbol."\t".$tranLength."\t".$GC."\n" );
            
                } else {
                    print("$chr\n") if $debug;
                }
            }
        } else {
            next;
        }
    } else {
        print ( "\nGene listed in $count row does not exist in Ensembl database\n\n" );
    }
}
close(GENE_FILE);
close(GENE_INFO);


#===============================================================================
#    End of main
#===============================================================================

exit;

#===============================================================================
#    Subroutines
#===============================================================================

sub usage () {
    print <<"EOS";
    
    usage: $0 -i gene_ref
    
    Input data
    -i gene_ref:     The name of the reference annotation GTF file with full path

EOS
exit;
}
