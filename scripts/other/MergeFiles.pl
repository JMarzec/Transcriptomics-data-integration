#!/usr/bin/perl -w
################################################################################
#
#   File name: MergeFiles.pl
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
#	Description: Script merging content of two files based on common content in specified column. The entire content of the primary file is recorded and only the rows with common content of specified column in the secondary file is merged
#
#	Command line use example: ./MergeFiles.pl -f1 /Users/marzec01/data/PhD/Transcriptomics_project/Comb_E-MEXP-1243_GSE17951_GSE45016_GSE55945_HGPINvsTumour_topTable.txt -f2 /Users/marzec01/data/PhD/Transcriptomics_project/PCa_assocaited_genes_DDPC.txt -c 1 -o /Users/marzec01/data/PhD/Transcriptomics_project/Comb_HGPINvsTumour_topTable.txt
#
#	-f1 file_1: Primary file
#	-f2 file_2: Secondary file to merged to the primary file based on the common content of specified column
#	-c column:  The number of column on which the merge should be based
#	-o out_file: The output file name with full path
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
my $inFile1;
my $inFile2;
my $inFile1_length;
my $colNo;
my $outFile;
my $rowsFile1 = 0;
my $rowsFile2 = 0;
my $rowsOverlap = 0;
my %file2Index = ();

while ($arg = shift) {
    
    if ($arg =~ /^-f1$/) {
        $inFile1 = shift;
    } elsif ($arg =~ /^-f2$/) {
        $inFile2 = shift;
    } elsif ($arg =~ /^-c$/) {
        $colNo = shift;
    } elsif ($arg =~ /^-o$/) {
        $outFile = shift;
    } elsif ($arg =~ /help/) {
        usage ();
    }
}

if (!$inFile1 or !$inFile2 or !$colNo or !$outFile) {
    usage ();
} else {
    $colNo--;
}

##### Store the content of the secondary file in a hash with specified column content as keys
open (FILE_2, $inFile2) or die $!;

#if ( $inFile1 && $inFile2 && $colNo && $outFile )
while (my $record = <FILE_2>) {
	
    $rowsFile2++;
    my @info = split(/\t/, $record);
    $inFile1_length = scalar(@info) - 1;
    
	if ( $info[$colNo] ) {
        ##### Remove the content of the specified column to avoid redundancy when the files are combined
        $record =~ s/$info[$colNo]\t//;
        $record =~ s/\s$//g;
        $info[$colNo] =~ s/\s$//g;
        $file2Index{ $info[$colNo] } = $record;
    }
}

close(FILE_2);

##### Compare files and report entire content of the primary file and only rows with common content in specified column of the secondary file
open (FILE_1, $inFile1) or die $!;

open (OUTFILE, ">$outFile") or die $!;

while (my $record = <FILE_1>) {
    
    $rowsFile1++;
    my @info = split(/\t/, $record);
    
    $record =~ s/\s$//g;
    $info[$colNo] =~ s/\s$//g;
    
	if ( exists $file2Index{ $info[$colNo] } ) {
        
        print ( OUTFILE $record."\t".$file2Index{ $info[$colNo] }."\n" );
        $rowsOverlap++;
    } else {
		
        print ( OUTFILE $record."\tNA" x $inFile1_length, "\n" );
	}		
}
close(FILE_1);
close(OUTFILE);

print ("\nNumber of rows in the primary file ".$rowsFile1."\nNumber of rows in the secondary file: ".$rowsFile2."\nNumber of common rows: ".$rowsOverlap."\n\n" );

#===============================================================================
#    End of main
#===============================================================================

exit;

#===============================================================================
#    Subroutines
#===============================================================================

sub usage () {
    print <<"EOS";
    
usage: $0 -f1 file_1 -f2 file_2 -c column -o out_file
    
    Input data
    -f1 file_1: Primary file
    -f2 file_2: Secondary file to merged to the primary file based on the common content of specified column
    -c column:  The number of column on which the merge should be based
        
    Output data
    -o out_file: The output name with full path

EOS
exit;
}

