#!/usr/bin/perl
################################################################################
#
#   File name: RNAseq_fastqAlign.pl
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
#	Description: Pipeline for RNA sequencing data alignment. It performs read alignment using Bowtie 2 read aligner (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) within TopHat splice junction mapper (http://tophat.cbcb.umd.edu/) and converts the ouptput to sorted bam files
#
#	Command line use example: ./RNAseq_fastqAlign.pl -t /scratch/jack/data/PhD/raw/Illum_HiSeq2000_RNAseq/E-MTAB-567/target_E-MTAB-567.txt -l PE -i /scratch/jack/genome_annotation/hg19 -o /scratch/jack/data/PhD/Transcriptomics_project -p 10
#
#	-t target_file:     Full path with name of the target file listing data to be analysed. The file is expected to include the following columns: (1) sample name; (2) file(s) name(s) (separated by comma, single-end experiments only), (3) target and (4) technical replicates indicated by same number
#	-l library_type:    Single- (SE) or paired-end (PE) library
#	-i ebwt_base:       The basename of the Bowtie index with full path
#	-o my_project:      Project workspace. This is the directory to which all the results will be written
#	-p threads_num:     Use this many threads to align reads. The default is 1
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
my $libType;
my $threadsNum = 1;
my $alignOut;
my $bowtieRefIndex;
my %datasetsList = ();
my $samples2excludeFile;
my $Project;
my $FastqFile;
my @FastqFileRoot;
my $FastqFile1;
my $FastqFile2;
my @FastqFiles;
my $FastqSEfiles;
my $sampleCount = 0;
my $AlignSummaryFile;
my $summaryFile;

while ($arg = shift) {
    if ($arg =~ /^-t$/) {
        $inFile = shift;
    } elsif ($arg =~ /^-l$/) {
        $libType = shift;
    } elsif ($arg =~ /^-i$/) {
        $bowtieRefIndex = shift;
    } elsif ($arg =~ /^-o$/) {
        $Project = shift;
    } elsif ($arg =~ /^-p$/) {
        $threadsNum = shift;
    }
}

if ($inFile and $bowtieRefIndex and $Project and $libType) {
    
    @dir = split('/', $inFile);
    $dir = join('/', @dir[0 .. $#dir-1]);
    $studyID = $dir[$#dir-1];
    $alignDir = join('/', $dir, "alignment");
    
    ##### Set/create a directory for the project
    unless(-d $Project){
        mkdir $Project or die;
    }
    
    
    ##### Add dataset to the lists of datasets to be integrated or create such file if it does not exist yet
    unless(-e $Project."/GenExpression_InputFiles.txt"){
        
        open ( DATASETS_LIST, ">$Project/GenExpression_InputFiles.txt");
        print( DATASETS_LIST "DatasetName\tTargetFile\tType\tPlatform\n".$studyID."\t".$inFile."\t".$libType."\tRNAseq\n" );
        
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
            print( DATASETS_LIST $studyID."\t".$inFile."\t".$libType."\tRNAseq\n" );
            
        } else {
            print( "\nDataset ".$studyID." is already in the datasets list.\n" );
        }
    }
    close( DATASETS_LIST);

    
    ##### Report detected outliers in the project workspace (THIS FEATURE WILL BE ADDED)
    $samples2excludeFile = join('', $Project, "/outliers_", $studyID, ".txt");
    open (EXCLUDE_FILE, ">$samples2excludeFile");
    
    print( EXCLUDE_FILE "DatasetName\tDataDir\tSamples2exclude\n$studyID\t$dir\t\n" );
    
    print( "\nProcessing dataset $studyID in $dir using $threadsNum threads\n\n" );
    
    ##### Set/create a directory for files to be generated
    unless(-d $alignDir){
        mkdir $alignDir or die;
    }
    
    ##### Create file for alignment summary
    $summaryFile = join('', $alignDir, '/', $studyID, "_summary.txt");
    
    if ( -e $summaryFile ) {
        
        $summaryFile = join('', $alignDir, '/', $studyID, "_summary_", int(rand(1000)), ".txt");
    }
    
    open (SUMMARY_FILE, ">$summaryFile");
    
    if ( $libType eq "PE" ) {
        print( SUMMARY_FILE "Sample\tForward_reads\tForward_mapped\tForward_mapped (%)\tReverse_reads\tReverse_mapped\tReverse_mapped (%)\tAligned_pairs\tDiscordant_pairs\tConcordant_pairs (%)\n" );
    } elsif ( $libType eq "SE" ) {
        print( SUMMARY_FILE "Sample\tTotal_reads\tTotal_mapped\tTotal_mapped (%)\n" );
    }
    
    open (INFILE, $inFile) or die $!;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
			
        $FastqFile = $info[ 1 ];

        if ( $FastqFile =~ m/fastq/ and $FastqFile !~ m/\_2\.fastq/ ) {
            
            $sampleCount++;
            
            @FastqFileRoot = split(/\_1\.fastq/, $FastqFile);
            $FastqFileRoot[ 0 ] =~ s/\.fastq.*?\.gz/\_/g;
            $FastqFileRoot[ 0 ] =~ s/\_$//g;
            $FastqFileRoot[ 0 ] =~ s/\,//g;
            
            $alignOut = join('/', $alignDir, $FastqFileRoot[ 0 ]);
            
            ##### Set/create a directory for files to be generated
            unless(-d $alignOut){
                mkdir $alignOut or die;
            }
            chdir( $alignOut );
            
            ##### Determine sequencing file names and perform alignment for paired-end...
            if ( $libType eq "PE" ) {
            	    
            	    $FastqFile1 = join('', $dir, '/', $FastqFileRoot[ 0 ], "_1.fastq", $FastqFileRoot[ 1 ]);
            	    $FastqFile2 = join('', $dir, '/', $FastqFileRoot[ 0 ], "_2.fastq", $FastqFileRoot[ 1 ]);
            
            	    ##### Align reads with Bowtie2 using TopHat (turned off the coverage-search algorithm, which is turned on as default for reads < 75bp, as this step takes too much tiem and memory):
                    print( "Mapping reads:\n$FastqFile1\n$FastqFile2\n\nWrite TopHat output to:\n$alignOut\n\n" );
            	    #system( "tophat --no-coverage-search -o $alignOut -p $threadsNum $bowtieRefIndex $FastqFile1 $FastqFile2" );
                    system( "tophat -o $alignOut -p $threadsNum $bowtieRefIndex $FastqFile1 $FastqFile2" );
            
            ##### ... or single-end experinet
            } elsif ( $libType eq "SE" ) {
            
            	    @FastqFiles = split(/\,/, $info[ 1 ]);
                
            	    if ( scalar(@FastqFiles) == 1 ) {
            	    	    
            	    	    $FastqSEfiles = $info[ 1 ];

            	    ##### Process multiple files for single sample
                    } else {
                            print( "Mapping reads:\n" );
            	    	    for (my $i = 0; $i < scalar(@FastqFiles); $i++) {
            	    	    
            	    	    	    $FastqFiles[$i] = join('', $dir, '/', $FastqFiles[$i]);
            	    	    	    print( "$FastqFiles[$i]\n" );
            	    	    }
            	    	    $FastqSEfiles = join(',', @FastqFiles);
                    }
 		    
            	    ##### Align reads with Bowtie2 using TopHat (turned off the coverage-search algorithm, which is turned on as default for reads < 75bp, as this step takes too much tiem and memory):
            	    print( "\nWrite TopHat output to:\n$alignOut\n\n" );
            	    #system( "tophat --no-coverage-search -o $alignOut -p $threadsNum $bowtieRefIndex $FastqSEfiles" );
                    system( "tophat -o $alignOut -p $threadsNum $bowtieRefIndex $FastqSEfiles" );
            }
            
            ##### Change Tophat2 output names:
            rename( "accepted_hits.bam", "$FastqFileRoot[ 0 ]_accepted_hits.bam" );
            rename( "deletions.bed", "$FastqFileRoot[ 0 ]_deletions.bed" );
            rename( "insertions.bed", "$FastqFileRoot[ 0 ]_insertions.bed" );
            rename( "junctions.bed", "$FastqFileRoot[ 0 ]_junctions.bed" );
            rename( "unmapped.bam", "$FastqFileRoot[ 0 ]_unmapped.bam" );
            rename( "prep_reads.info", "$FastqFileRoot[ 0 ]_prep_reads.info" );
            
            print( "TopHat2 has finished data processing for run nr $sampleCount at ");
            system( "date" );
            print( "\n" );
            
            ##### Index alignment:
            print( "Indexing the alignment results...\n\n" );
            system( "samtools index $FastqFileRoot[ 0 ]_accepted_hits.bam" );
            
            ##### Get simple stats from TopHat2 output file 'align_summary.txt'
            $AlignSummaryFile = $alignOut . "/align_summary.txt";

            open (ALIGN_SUMMARY_FILE, $AlignSummaryFile) or die $!;
            
            print( SUMMARY_FILE "$FastqFileRoot[ 0 ]\t" );
            
            while (my $record = <ALIGN_SUMMARY_FILE>) {
                 
                if ( $record =~ m/Input:\s+(\d+)/ or $record =~ m/Mapped:\s+(\d+)/ or $record =~ m/Aligned pairs:\s+(\d+)/ or $record =~ m/and:\s+(\d+)/ or $record =~ m/^(\d+\.\d+%) concordant/ ) {
                    
                    print( SUMMARY_FILE "$1\t" );

                    if ( $record =~ m/\((\d+\.\d+)/ and $record !~ m/and:\s+(\d+)/) {
                        
                        print( SUMMARY_FILE "$1%\t" );
                    }
                }
            }
            print( SUMMARY_FILE "\n" );
            close (ALIGN_SUMMARY_FILE);

        }
    }
    close (INFILE);
    close (EXCLUDE_FILE);
    close (SUMMARY_FILE);
    
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
    
usage: $0 -t target_file -l library_type -i ebwt_base -o my_project [options] 
    
    Input data
    -t target_file:  Target file name with full path to data to be analysed.
                     The file is expected to include the following columns:
                     (1) sample name; (2) file(s) name(s) (separated by comma, 
                     single-end experiments only) and (3) target
    -l library_type: Single- (SE) or paired-end (PE) library
    -i ebwt_base:    The basename of the Bowtie index with full path
    
    Output data
    -o my_project:   Project workspace. This is the directory to which all the
                     results will be written
    
    Options
    -p threads_num:  Use this many threads to align reads. The default is 1
    
EOS
exit;
}
