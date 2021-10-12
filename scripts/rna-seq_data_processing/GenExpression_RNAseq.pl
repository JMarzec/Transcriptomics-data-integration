#!/usr/bin/perl
################################################################################
#
#   File name: GenExpression_RNAseq.pl
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
#	Description: Pipeline for calculating read counts for genes in GFF/GTF format using htseq-count tool within HTSeq Python package ( http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html ). As an input the TopHat2 alignment results are required.
#
#	Command line use example: ./GenExpression_RNAseq.pl -t /scratch/jack/data/PhD/raw/Illum_HiSeq2000_RNAseq/ICGC/target_ICGC.txt -i /scratch/jack/genome_annotation/ensembl_Homo_sapiens.GRCh37.74.gtf -o /scratch/jack/data/PhD/Transcriptomics_project
#
#	-t target_file:     Full path with name of the target file listing data to be analysed. The file is expected to include the following columns: (1) sample name; (2) file name, (3) target and (4) technical replicates indicated by same number
#	-i gene_ref:        The name of the reference annotation GFF/GTF file with full path
#	-o my_project:      Project workspace. This is the directory to which all the results will be written
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
my $sampleName;
my $alignOut;
my $Project;
my %datasetsList = ();
my $samples2excludeFile;
my %samples2exclude;
my $refAnnot;
my $FastqFile;
my @FastqFileRoot;
my $sampleCount = 0;

while ($arg = shift) {
    if ($arg =~ /^-t$/) {
        $inFile = shift;
    } elsif ($arg =~ /^-i$/) {
        $refAnnot = shift;
    } elsif ($arg =~ /^-o$/) {
        $Project = shift;
    }
}

if ($inFile and $refAnnot and $Project) {
    
    @dir = split('/', $inFile);
    $dir = join('/', @dir[0 .. $#dir-1]);
    $studyID = $dir[$#dir-1];
    $alignDir = join('/', $dir, "alignment");
    
} else {
    usage ();
}

##### Store sample to be excluded in a hash
$samples2excludeFile = join('', $Project, "/outliers_", $studyID, ".txt");

open (INFILE, $samples2excludeFile) or die $!;

while (my $record = <INFILE>) {
    
    chomp $record;
    
    my @info = split(/\t/, $record);
    
    ##### Create hash with samples to be excluded
    if ( $info[2] ) {
        
        my @samples = split(/,/, $info[2]);
        
        foreach (@samples) {
            $samples2exclude{ $_ } = $_;
        }
    }
}
close(INFILE);


##### Add dataset to the lists of datasets to be integrated or create such file if it does not exist yet
unless(-e $Project."/GenExpression_InputFiles.txt"){
    
    open ( DATASETS_LIST, ">$Project/GenExpression_InputFiles.txt");
    print( DATASETS_LIST "DatasetName\tTargetFile\tType\tPlatform\n".$studyID."\t".$inFile."\traw\tRNAseq\n" );
    
} else {
    open (DATASETS_LIST, $Project."/GenExpression_InputFiles.txt") or die $!;
    
    ##### Check if the dataset is already listed in the file
    while (my $record = <DATASETS_LIST>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $datasetsList{ $info[0] } = $info[0];
    }
    close( DATASETS_LIST);
    
    if ( !exists $datasetsList{ $studyID } ) {
        
        open (DATASETS_LIST, ">>$Project/GenExpression_InputFiles.txt");
        print( DATASETS_LIST $studyID."\t".$inFile."\traw\tRNAseq\n" );
        
    } else {
        print( "\nDataset ".$studyID." is already in the datasets list.\n\n" );
    }
}
close( DATASETS_LIST);


if ($inFile) {
    
    print( "\nProcessing dataset $studyID in $dir\n\n" );
    
    ##### Set a directory with all alignment results
    unless(-d $alignDir){
        print( "alignment results directory $alignDir does not exist!\n") or die;
    }
    
    open (INFILE, $inFile) or die $!;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
			
        $FastqFile = $info[ 1 ];
        $sampleName = $info[ 0 ];
        
        if ( exists $samples2exclude{ $sampleName }  ) {
        
            print("Sample $sampleName is skipped!\n\n");
            
        } elsif ( $FastqFile =~ m/fastq/ and $FastqFile !~ m/\_2\.fastq/ ) {
            
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
            
            ##### Sort BAM files and calculate read count for each genes using htseq-count
            print( "Sorting $FastqFileRoot[ 0 ] BAM file\n\n" );
            system( "samtools sort -n $FastqFileRoot[ 0 ]_accepted_hits.bam $FastqFileRoot[ 0 ]_accepted_hits_sorted");
            
            print( "Calculating read counts for $FastqFileRoot[ 0 ]\n\nWrite HTseq-count output to:\n$alignOut\n\n" );
            system( "python -m HTSeq.scripts.count --stranded=no --format=bam --quiet $FastqFileRoot[ 0 ]_accepted_hits_sorted.bam $refAnnot >> $FastqFileRoot[ 0 ]_accepted_hits_count.txt" );
            
            print( "Read count for sample nr $sampleCount has finished at ");
            system( "date" );
            print( "\n\n" );
            
            ##### Remove sorted BAM file
            system( "rm $FastqFileRoot[ 0 ]_accepted_hits_sorted.bam" );
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
    
usage: $0 -t target_file -i gene_ref -o my_project
    
    Input data
    -t target_file:  Target file name with full path to data to be analysed.
                     The file is expected to include the following columns:
                     (1) sample name; (2) file name and (3) target
    -i gene_ref:     The name of the reference annotation GTF file with full path
    
    Output data
    -o my_project:   Project workspace. This is the directory to which all the
                     results will be written
    
EOS
exit;
}
