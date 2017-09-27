#!/usr/bin/perl
################################################################################
#
#   File name: QC_RNAseq.pl
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
#	Description: Pipeline for assessing technical quality of a RNA-seq data using htseq-qa tool within HTSeq Python package ( http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html ). As an input the TopHat2 alignment results are required. It produces a PDF file with useful plots to assess the technical quality of a run.
#
#	Command line use example: ./QC_RNAseq.pl -t /scratch/jack/data/PhD/raw/Illum_HiSeq2000_RNAseq/E-MTAB-567/target_E-MTAB-567.txt -m 40
#
#	-t target_file:     Full path with name of the target file listing data to be analysed. The file is expected to include the following columns: (1) sample name; (2) file name, (3) target and (4) technical replicates indicated by same number
#	-m max_qual:        The maximum quality score that appears in the data (default: 40)
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
my $studyID;
my $alignDir;
my $maxQual = 40;
my $sampleName;
my $alignOut;
my $QCdir;
my $QCout;
my $FastqFile;
my @FastqFileRoot;
my $sampleCount = 0;

while ($arg = shift) {
    if ($arg =~ /^-t$/) {
        $inFile = shift;
    } elsif ($arg =~ /^-m$/) {
        $maxQual = shift;
    }
}

if ($inFile) {
    
    @dir = split('/', $inFile);
    $dir = join('/', @dir[0 .. $#dir-1]);
    $studyID = $dir[$#dir-1];
    $alignDir = join('/', $dir, "alignment");
    
} else {
    usage ();
}

if ($inFile) {
    
    print( "\nProcessing dataset $studyID in $dir\n\n" );
    
    ##### Set a directory with all alignment results
    unless(-d $alignDir){
        print( "alignment results directory $alignDir does not exist!\n") or die;
    }
    
    ##### Set/create a directory for QC plots
    $QCdir = join('/', $dir, join('_', "QC", $studyID) );
    
    unless(-d $QCdir){
        mkdir $QCdir or die;
    }
    
    open (INFILE, $inFile) or die $!;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
			
        $FastqFile = $info[ 1 ];
        $sampleName = $info[ 0 ];
        
        if ( $FastqFile =~ m/fastq/ and $FastqFile !~ m/\_2\.fastq/ ) {
            
            $sampleCount++;
            
            @FastqFileRoot = split(/\_1\.fastq/, $FastqFile);
            $FastqFileRoot[ 0 ] =~ s/\.fastq.*?\.gz/\_/g;
            $FastqFileRoot[ 0 ] =~ s/\_$//g;
            $FastqFileRoot[ 0 ] =~ s/\,//g;
            
            $alignOut = join('/', $alignDir, $FastqFileRoot[ 0 ]);
            
            ##### Set a directory with sample alignment results
            unless(-d $alignOut){
                print( "alignment results directory for $FastqFileRoot[ 0 ] does not exist!\n") or die;
            }
            
            chdir( $alignOut );
            
            ##### Convert mapped and unmapped BAM files to SAM format and produce QC plots using htseq-qa
            print( "Converting $FastqFileRoot[ 0 ] BAM files to SAM format\n\n" );
            system( "samtools view -h $FastqFileRoot[ 0 ]_accepted_hits.bam > $FastqFileRoot[ 0 ]_accepted_hits.sam");
            system( "samtools view -h $FastqFileRoot[ 0 ]_unmapped.bam > $FastqFileRoot[ 0 ]_unmapped.sam");
            
            print( "Merging SAM files with mapped and unmapped reads\n\n" );
            system( "cat $FastqFileRoot[ 0 ]_accepted_hits.sam $FastqFileRoot[ 0 ]_unmapped.sam > $FastqFileRoot[ 0 ].sam");
            
            ##### Remove SAM files
            system( "rm $FastqFileRoot[ 0 ]_accepted_hits.sam" );
            system( "rm $FastqFileRoot[ 0 ]_unmapped.sam" );
            
            $QCout = join('/', $dir, join('_', "QC", $studyID), join('.', $FastqFileRoot[ 0 ], "pdf") );
            
            print( "Producing QC plots for $FastqFileRoot[ 0 ]\n\nWrite HTseq-qa output to:\n$QCdir\n\n" );
            system( "python -m HTSeq.scripts.qa -m $maxQual -o $QCout $FastqFileRoot[ 0 ].sam" );
            
            print( "QC check for sample nr $sampleCount has finished at ");
            system( "date" );
            print( "\n\n" );
            
            ##### Remove SAM file
            system( "rm $FastqFileRoot[ 0 ].sam" );
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
    
usage: $0 -t target_file -i ebwt_base -o my_project [options]
    
    Input data
    -t target_file:  Target file name with full path to data to be analysed.
                     The file is expected to include the following columns:
                     (1) sample name; (2) file name and (3) target
    
    Options
    -m max_qual:     The maximum quality score that appears in the data (default: 40).
    
EOS
exit;
}
