#!/usr/bin/perl
################################################################################
#
#   File name: FastQC_RNAseq.pl
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
#	Description: Pipeline for assessing technical quality of a RNA-seq data using FastQC tool ( http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ ). As an input the raw fastq files are used. It produces a PDF file with useful plots to assess the technical quality of a run. It creates an HTML report with embedded graphs, but also a zip file containing individual graph files and additional data files containing the raw data from which plots were drawn. User has to specify the location of the FastQC tool.
#
#	Command line use example: ./FastQC_RNAseq.pl -l /scratch/jack/data/PhD/raw/Illum_HiSeq2000_RNAseq/E-MTAB-567/E-MTAB-567_fastq_list.txt -m 1 -d /data/home/hfw456/applications/FastQC.app
#
#	-l fastq_list:      Full path with name of a file listing fastq files to be analysed. The file is expected to include one fastq file per row
#	-m threads:         Number of files which can be processed simultaneously.  Each thread will be allocated 250MB of memory (default: 1)
#	-d FastQC:          Full path to the FastQC tool
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

our @ARGV;
my $arg;
my $inFile;
my @dir;
my $dir;
my $FastQCdir;
my $studyID;
my $threads = 1;
my $QCdir;
my $FastqFile;
my $fileCount = 0;

while ($arg = shift) {
    if ($arg =~ /^-l$/) {
        $inFile = shift;
    } elsif ($arg =~ /^-m$/) {
        $threads = shift;
    } elsif ($arg =~ /^-d$/) {
        $FastQCdir = shift;
    }
}

if ($inFile and $FastQCdir) {
    
    @dir = split('/', $inFile);
    $dir = join('/', @dir[0 .. $#dir-1]);
    $studyID = $dir[$#dir-1];
    
} else {
    usage ();
}

if ($inFile) {
    
    print( "\nProcessing dataset $studyID in $dir\n\n" );
    
    ##### Set/create a directory for QC plots
    $QCdir = join('/', $dir, join('_', "FastQC", $studyID) );
    
    unless(-d $QCdir){
        mkdir $QCdir or die;
    }
    
    open (INFILE, $inFile) or die $!;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $FastqFile = $info[ 0 ];
        
        
        system("$FastQCdir/Contents/Resources/Java/fastqc $dir/$FastqFile --outdir $QCdir");
        
        if ( $FastqFile =~ m/fastq/ and $FastqFile !~ m/\_1\.fastq/ ) {
            
            $fileCount++;
            
            print( "\nQC check for file nr $fileCount has finished at ");
            system( "date" );
            print( "\n\n" );
        }
    }
    close(INFILE);
    
} else {
    usage ();
}

#===============================================================================
#    End of main
#===============================================================================

exit;

#===============================================================================
#    Subroutines
#===============================================================================

sub usage () {
    print <<"EOS";
    
usage: $0 -l fastq_list -i ebwt_base -o my_project
    
    Input data
    -l fastq_list:   Full path with name of a file listing fastq files to be
                     analysed. The file is expected to include one fastq file
                     per row
    
    Options
    -m :             Number of files which can be processed simultaneously.
                     Each thread will be allocated 250MB of memory (default: 1)
    -d FastQC:       Full path to the FastQC tool
    
EOS
exit;
}
