#!/usr/bin/perl -w
################################################################################
#
#   File name: Convert_RSEMprocessed.pl
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
#	Description: Script preparing gene expression files processed with transcript quantification software RSEM ( http://deweylab.biostat.wisc.edu/rsem/ ) to be merged with other RNA-seq datasets. It reads in [GTF]_gene_info.txt file to map gene symbols with Ensembl gene IDs retrieved with Perl script 'Get_gene_info.pl'
#
#	Command line use example: ./Convert_RSEMprocessed.pl -t /scratch/jack/data/PhD/processed/Illum_HiSeq2000_RNAseq/TCGA/target_TCGA.txt -i /scratch/jack/genome_annotation/ensembl_Homo_sapiens.GRCh37.74.gtf.gene_info.txt
#
#	-t target_file:     Full path with name of the target file listing data to be analysed. The file is expected to include the following columns: (1) sample name; (2) file name and (3) target
#	-i gene_info:       The name of the [GTF]_gene_info.txt file with full path
#
################################################################################

use strict;
use warnings;

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
my $inFile;
my $geneFile;
my @dir;
my $dir;
my %geneInfo = ();

while ($arg = shift) {
    
    if ($arg =~ /^-t$/) {
        $inFile = shift;
    } elsif ($arg =~ /^-i$/) {
        $geneFile = shift;
    } elsif ($arg =~ /help/) {
        usage ();
    }
}

if ( $inFile and $geneFile ) {
    
    @dir = split('/', $inFile);
    $dir = join('/', @dir[0 .. $#dir-1]);
    
} else {
    usage ();
}


##### Get Ensembl gene IDs from [GTF]_gene_info.txt file
open (GENE_FILE, $geneFile) or die $!;

while (my $record = <GENE_FILE>) {
    
    chomp $record;
    
    my @info = split(/\t/, $record);
    $geneInfo{ $info[1] } = $info[0];
}
close(GENE_FILE);


##### Get information about samples to be analysed
open (INFILE, $inFile) or die $!;

##### Skip the first line
<INFILE>;

while (my $record = <INFILE>) {
    
    chomp $record;
    
    my @info = split(/\t/, $record);
    
    ##### Open RNA-seq expression file with read counts
    my $targetFile = join('/', $dir, $info[1]);
    
    open (TARGET_FILE, $targetFile) or die $!;
    print("\nProcessing $info[1] ...\n\n");
    
    ##### Create RNA-seq expression file for FPKM values
    my $countFile = $targetFile;
    $countFile =~ s/\.results/\_count\.txt/;
    
    open ( COUNT_FILE, ">$countFile" );

    ##### Skip the first line
    <TARGET_FILE>;
    
    while (my $record = <TARGET_FILE>) {
        
        chomp $record;
        
        my @sampleInfo = split(/\t/, $record);
        
        my @geneSymbol = split(/\|/, $sampleInfo[0]);
        
        ##### Retrieve gene information
        if ( exists  $geneInfo{ $geneSymbol[0] } ) {
            
            ##### Write the Ensembl gene IDs along with read-counts into a file
            print ( COUNT_FILE $geneInfo{ $geneSymbol[0] }."\t".$sampleInfo[1]."\n" );
        }
    }
    close(TARGET_FILE);
    close(COUNT_FILE);
}
close(INFILE);


#===============================================================================
#    End of main
#===============================================================================

exit;

#===============================================================================
#    Subroutines
#===============================================================================

sub usage () {
    print <<"EOS";
    
    usage: $0 -t target_file -i gene_info
    
    Input data
    -t target_file:  Target file name with full path to data to be analysed.
                     The file is expected to include the following columns:
                     (1) sample name; (2) file name and (3) target
    
    -i gene_info:    The name of the [GTF]_gene_info.txt file with full path
   
EOS
exit;
}
