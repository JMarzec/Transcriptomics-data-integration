#!/usr/bin/perl -w
################################################################################
#
#   File name: GenExpressionComb_RNAseq.pl
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
#	Description: Script combining RNA-seq gene expression level data (read counts) for specified samples into a single gene-by-sample matrix. !!!CURRENTLY WORKS ON PAIRED-END RNA-SEQ DATA ONLY!!!
#
#	Command line use example: ./GenExpressionComb_RNAseq.pl -o /scratch/jack/data/PhD/Transcriptomics_project
#
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
my $Project;
my @dir;
my $dir;
my $studyID;
my $sampleName;
my $FastqFile;
my %geneList=();
my $sampleNo=0;
my $sampleList='';
my $countMatrix='';
my $samples2excludeFile;
my %samples2exclude;

while ($arg = shift) {
    
    if ($arg =~ /^-o$/) {
        $Project = shift;
    } elsif ($arg =~ /help/) {
        usage ();
    }
}

if ( $Project ) {
    
    ##### Get list of datasets to be analysed
    if (-e $Project."/GenExpression_InputFiles.txt") {
        
        open (DATASETS_LIST, $Project."/GenExpression_InputFiles.txt") or die $!;

        ##### Skip the first line
        <DATASETS_LIST>;

        while (my $dataset = <DATASETS_LIST>) {
    
            chomp $dataset;
    
            my @datasetInfo = split(/\s/, $dataset);
            
            ##### Consider RNA-seq datasets only
            if ( $datasetInfo[3] eq "RNAseq" ) {

                @dir = split('/', $datasetInfo[1]);
                $dir = join('/', @dir[0 .. $#dir-1]);
                $studyID = $dir[$#dir-1];
                
                print( "\nProcessing dataset $studyID...\n");
                
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
                
                open (INFILE, $datasetInfo[1]) or die $!;
                
                ##### Skip the first line
                <INFILE>;
                
                while (my $record = <INFILE>) {
                    
                    chomp $record;
                    
                    my @info = split(/\t/, $record);
                    
                    ##### Merge read-count data from processed files
                    if ( $datasetInfo[2] eq "processed" ) {
                        
                        $sampleName = $info[ 0 ];
	
                        if ( exists $samples2exclude{ $sampleName }  ) {
                            
                            print("Sample $sampleName is skipped!\n\n");
                            
                        } else {

                        ##### Capture sample names to be written in the output file header
                        $sampleList = $sampleList."\t".$info[0];
                        
                        my $countFile = $dir."/".$info[1];
                        $countFile =~ s/\.results/\_count\.txt/;
                        
                        open (COUNT_FILE, $countFile) or die $!;
                        
                        while (my $countInfo = <COUNT_FILE>) {
                            
                            chomp $countInfo;
                                
                            my @sampleInfo = split(/\t/, $countInfo);
                            my $readCount = $sampleInfo[1];
                            $readCount =~ s/\s$//;
                                
                            ##### Report all captured genes and store their expression (read counts) per sample in separate hashes
                            $geneList{ $sampleInfo[0] }[$sampleNo] = $readCount;
                        }
                        close(COUNT_FILE);
                        
                        $sampleNo++;
                        }
                        
                    ##### Merge read-count data generated from raw fastq files
                    } elsif ($datasetInfo[2] eq "raw") {
                
                        $FastqFile = $info[ 1 ];
                        $sampleName = $info[ 0 ];
                
                        if ( exists $samples2exclude{ $sampleName }  ) {
                            
                            print("Sample $sampleName is skipped!\n\n");
                            
                        } elsif ( $FastqFile =~ m/fastq/ and $FastqFile !~ m/\_2\.fastq/ ) {
                    
                            my @FastqFileRoot = split(/\_1\.fastq/, $FastqFile);
                            $FastqFileRoot[ 0 ] =~ s/\.fastq.*?\.gz/\_/g;
                            $FastqFileRoot[ 0 ] =~ s/\_$//g;
                            $FastqFileRoot[ 0 ] =~ s/\,//g;
                    
                            my $alignOut = join('/', $dir, "alignment", $FastqFileRoot[ 0 ]);
                            
                            ##### Set/create a directory for files to be generated
                            unless(-d $alignOut){
                                mkdir $alignOut or die;
                            }
                                
                            ##### Capture sample names to be written in the output file header
                            $sampleList = $sampleList."\t".$info[0];

                            my $countFile = $alignOut."/$FastqFileRoot[ 0 ]_accepted_hits_count.txt";
                            
                            open (COUNT_FILE, $countFile) or die $!;
    
                            while (my $countInfo = <COUNT_FILE>) {
        
                                chomp $countInfo;
        
                                ##### Omit lines with special counters listed by HTSeq
                                if ( $countInfo !~ m/\_\_/ ) {
                                    
                                    my @sampleInfo = split(/\t/, $countInfo);
                                    my $readCount = $sampleInfo[1];
                                    $readCount =~ s/\s$//;
            
                                    ##### Report all captured genes and store their expression (read counts) per sample in separate hashes
                                    $geneList{ $sampleInfo[0] }[$sampleNo] = $readCount;
                                }
                            }
                            $sampleNo++;
                            close(COUNT_FILE);
                        }
                    }
                }
            close(INFILE);
            }
        }
        close(DATASETS_LIST);
        
    } else {
        print( "\n".$Project."/GenExpression_InputFiles.txt does not exist!\n\n");
        exit;
    }
} else {
    usage ();
}

print( "\nMerging RNA-seq datasets...\n\n");
$countMatrix = join('/', $Project, "RNAseq_count_matrix.txt");

open ( COUNT_FILE, ">$countMatrix" );

##### Write the header with sample names
print ( COUNT_FILE "gene_id".$sampleList."\n" );

##### Create single gene-by-sample read-count matrix and save it into a file
for my $geneID ( keys %geneList ) {
    
    print ( COUNT_FILE $geneID );
    
    for ( my $i=0; $i<$sampleNo; $i++ ) {
        
        if ( ${ $geneList{$geneID} }[$i] ) {
            print ( COUNT_FILE "\t".${ $geneList{$geneID} }[$i] );
        } else {
            print ( COUNT_FILE "\t0" );
        }
    }
    print ( COUNT_FILE "\n" );
}

close(COUNT_FILE);

#===============================================================================
#    End of main
#===============================================================================

exit;

#===============================================================================
#    Subroutines
#===============================================================================

sub usage () {
    print <<"EOS";
    
    usage: $0 -o my_project
    
    Output data
    -o my_project:   Project workspace. This is the directory to which all the results
                     will be written
   
EOS
exit;
}
