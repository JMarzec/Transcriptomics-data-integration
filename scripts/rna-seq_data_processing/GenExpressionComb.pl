#!/usr/bin/perl
################################################################################
#
#   File name: genExpression.pl
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
#	Description: Pipeline combining cross-platform datasets
#
#	Command line use example: ./GenExpressionComb.pl -i /scratch/jack/data/PhD/GenExpressionComb_InputFiles.txt -o /scratch/jack/data/PhD/Transcriptomics_project
#	
#	-i info_file:       File name with full path that lists the datasets to be analysed. The file is expected to include the following columns: (1) dataset name; (2) target file name with full path; (3) data type (raw or processed) and (4) platform used
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
my $projectName;
my $targeFile;
my %arrays = ();
my %samples = ();

if (defined $arg)  {###
    print $arg, "\n";####
}#####

while ($arg = shift) {
    
    if ($arg =~ /^-i$/) {
        $inFile = shift;
    } elsif ($arg =~ /^-o$/) {
        $projectName = shift;
    } elsif ($arg =~ /help/) {
        usage ();
    }
}

if ($inFile and $projectName) {
    
} else {
    usage ();
}

##### Get information about datasets to be analysed
open (INFILE, $inFile) or die $!;

##### Create the target file for all datasets
$targeFile = join('', $projectName, "/target.txt");
open (TARGET_OUTFILE, ">$targeFile");
print("$targeFile\n\n");

##### Skip the first line
<INFILE>;

while (my $record = <INFILE>) {
    
    chomp $record;
    
    my @info = split(/\t/, $record);
    
    ##### Create hash with platforms used
    $arrays{ $info[3] } = $info[3];
    print("$info[3]\n");
    
    ##### Read the target files to get the information about samples
    open (TARGET_FILE, $info[1]) or die $!;
    ##### Skip the first line
    <TARGET_FILE>;
    
    while (my $record = <TARGET_FILE>) {
        
        chomp $record;
        
        my @targetInfo = split(/\t/, $record);
        
        ##### Create hash with samples to be analysed
        $samples{ $targetInfo[1] } = $targetInfo[1];
        #print("$targetInfo[1]\n");
    }
    close(TARGET_FILE);
    
    ##### Create hash with arrays to be analysed
    $arrays{ $info[3] } = $info[3];
}

close(INFILE);
close(TARGET_OUTFILE);

#===============================================================================
#    End of main
#===============================================================================

exit;

#===============================================================================
#    Subroutines
#===============================================================================

sub usage () {
    print <<"EOS";
    
    usage: $0 -i info_file -o my_project
    
    Input data
    -i info_file:  File name with full path that lists the datasets to be analysed.
                   The file is expected to include the following columns:
                   (1) dataset name; (2) target file name with full path; (3) data
                   type (raw or processed) and (4) platform used
    
    Output data
    -o my_project: Project workspace. This is the directory to which all the results
                   will be written
   
EOS
exit;
}


