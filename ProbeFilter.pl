#!/usr/bin/perl
################################################################################
#
#   File name: ProbeFilter.pl
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
#	Description: Script filtering out unreliable probesets based on cigar strings derived from Ensembl microarray probe mapping. It uses ArrayAnnot.pl annotation files as input files. For Affymetrix 3' IVT arrays, only probesets with < 3 probes with mismatch/insertion/deletion or mapping to none or more than one gene are retained.  For Affymetrix exon arrays, only probesets with < 2 probes mismatch/insertion/deletion and mapping to one gene are retained. For Illumina array, only probes with perfect match and mapping to exactly one gene are retained.
#
#	Command line use example: ./ProbeFilter.pl -i /scratch/jack/data/PhD/Transcriptomics_project/GenExpression_InputFiles.txt
#
#	-i info_file:       File name with full path that lists the datasets to be analysed. The file is expected to include the following columns: (1) dataset name; (2) target file name with full path; (3) data type (raw or processed) and (4) platform used
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
my %Arrays =();
my $annotFile;
my $flatFile;
my $summaryFile;

while ($arg = shift) {
    
    if ($arg =~ /^-i$/) {
        $inFile = shift;
        @dir = split('/', $inFile);
        $dir = join('/', @dir[0 .. $#dir-1]);
    } elsif ($arg =~ /help/) {
        usage ();
    }
}

if (!$inFile) {
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

################################################################################
#
#			AFFYMETRIX HuEx-1ST v2
#
################################################################################

my %Probes =();
my %Probes2rm =();
my $probeUniqID = '';
my %probesUniq = ();
my %probesPerSet = ();
my %nProbesPerSet = ();
my %nXprobesPerSet = ();

if ( exists  $Arrays{ 'Affy_HuEx1ST' } ) {
    
    $annotFile = $dir."/Affy_HuEx1ST_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;    
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $probeUniqID = $info[0].$info[1];
        
        ##### Count number of probes in each probeset
        if ( !exists $probesUniq{ $probeUniqID } ) {
            
            $probesPerSet{ $info[0] }[0]++;
            
            ##### Count probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            if ( $info[5] ne "25=" or $info[6] != 1 ) {
                
                $probesPerSet{ $info[0] }[1]++;
                
            ##### Count probes mapping to other gene than remaining probes within the probeset
            } elsif ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                
                $probesPerSet{ $info[0] }[1]++;
            }
        }
        
        if ( $info[3] ne "-" ) {
            ##### Report probesets containing probes that map to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                    
                $Probes2rm{ $info[0] }=10;
            
            ##### Count probes, within a probeset, with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[5] ne "25=" or $info[6] != 1 ) {
            
                $Probes2rm{ $info[0] }++;
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[3];
            $Probes{ $info[0] }[2] = $info[4];
        }
        $probesUniq{ $probeUniqID } = $probeUniqID;
    }
    
    close(INFILE);
    
    $flatFile = $dir."/Affy_HuEx1ST.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probesets
    while ( (my $key, my $value) = each %Probes ) {
        
        ##### Report only probesets with < 2 probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
        if ( exists  $Probes2rm{ $key } and $Probes2rm{ $key } > 1 ) {
            
            delete $Probes{ $key };
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);

################################################################################
#			Mapping summary for AFFYMETRIX HuEx-1ST v2
################################################################################
    
    $summaryFile = $dir."/Affy_HuEx1ST.perProbeset_summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    print ( SUMMARY_FILE "probeset\tprobes_number\tunreliable_probes_number\n" );
    
    ##### Write mapping summary per each probeset
    while ( (my $key, my $value) = each %probesPerSet ) {
        
        ##### Count probesets consisting of n probes
        $nProbesPerSet{ $probesPerSet{ $key }[0] }++;
        
        ##### Report probesets with > 1 unrealiable probes
        if ($probesPerSet{ $key }[1] ) {
        
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t".$probesPerSet{ $key }[1]."\n" );
            
            ##### Count probesets with n unrealiable probes
            $nXprobesPerSet{ $probesPerSet{ $key }[1] }++;
            
        ##### Report probesets with no unrealiable probes
        } else {
            
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t0\n" );
            
            ##### Count probesets with no unrealiable probes
            $nXprobesPerSet{ 0 }++;
        }
    }
    close(SUMMARY_FILE);
    
    
    $summaryFile = $dir."/Affy_HuEx1ST.summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    
    ##### Report probesets with n probes
    print ( SUMMARY_FILE "probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nProbesPerSet ) {
            
            print( SUMMARY_FILE $key."\t".$nProbesPerSet{ $key }."\n" )
    }
    
    ##### Report probesets with n unrealiable probes
    print ( SUMMARY_FILE "\nunreliable_probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nXprobesPerSet ) {
        
        print( SUMMARY_FILE $key."\t".$nXprobesPerSet{ $key }."\n" )
    }
    close(SUMMARY_FILE);
}             

################################################################################
#
#			AFFYMETRIX HuGene-1ST v1
#
################################################################################

%Probes =();
%Probes2rm =();
$probeUniqID = '';
%probesUniq = ();
%probesPerSet = ();
%nProbesPerSet = ();
%nXprobesPerSet = ();

if ( exists  $Arrays{ 'Affy_HuGene1ST' } ) {
    
    $annotFile = $dir."/Affy_HuGene1ST_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;    
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $probeUniqID = $info[0].$info[1];
        
        ##### Count number of probes in each probeset
        if ( !exists $probesUniq{ $probeUniqID } ) {
            
            $probesPerSet{ $info[0] }[0]++;
            
            ##### Count probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            if ( $info[5] ne "25=" or $info[6] != 1 ) {
                
                $probesPerSet{ $info[0] }[1]++;
                
            ##### Count probes mapping to other gene than remaining probes within the probeset
            } elsif ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                
                $probesPerSet{ $info[0] }[1]++;
            }
        }
        
        if ( $info[3] ne "-" ) {
            ##### Report probesets containing probes that map to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                    
                $Probes2rm{ $info[0] }=10;
            
            ##### Count probes, within a probeset, with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[5] ne "25=" or $info[6] != 1 ) {
            
                $Probes2rm{ $info[0] }++;
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[3];
            $Probes{ $info[0] }[2] = $info[4];
        }
        $probesUniq{ $probeUniqID } = $probeUniqID;
    }
    
    close(INFILE);
    
    $flatFile = $dir."/Affy_HuGene1ST.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probesets
    while ( (my $key, my $value) = each %Probes ) {
        
        ##### Report only probesets with < 2 probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
        if ( exists  $Probes2rm{ $key } and $Probes2rm{ $key } > 1 ) {
            
            delete $Probes{ $key };
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);

################################################################################
#			Mapping summary for AFFYMETRIX HuGene-1ST v1
################################################################################
    
    $summaryFile = $dir."/Affy_HuGene1ST.perProbeset_summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    print ( SUMMARY_FILE "probeset\tprobes_number\tunreliable_probes_number\n" );
    
    ##### Write mapping summary per each probeset
    while ( (my $key, my $value) = each %probesPerSet ) {
        
        ##### Count probesets consisting of n probes
        $nProbesPerSet{ $probesPerSet{ $key }[0] }++;
        
        ##### Report probesets with > 1 unrealiable probes
        if ($probesPerSet{ $key }[1] ) {
        
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t".$probesPerSet{ $key }[1]."\n" );
            
            ##### Count probesets with n unrealiable probes
            $nXprobesPerSet{ $probesPerSet{ $key }[1] }++;
            
        ##### Report probesets with no unrealiable probes
        } else {
            
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t0\n" );
            
            ##### Count probesets with no unrealiable probes
            $nXprobesPerSet{ 0 }++;
        }
    }
    close(SUMMARY_FILE);
    
    
    $summaryFile = $dir."/Affy_HuGene1ST.summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    
    ##### Report probesets with n probes
    print ( SUMMARY_FILE "probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nProbesPerSet ) {
            
            print( SUMMARY_FILE $key."\t".$nProbesPerSet{ $key }."\n" )
    }
    
    ##### Report probesets with n unrealiable probes
    print ( SUMMARY_FILE "\nunreliable_probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nXprobesPerSet ) {
        
        print( SUMMARY_FILE $key."\t".$nXprobesPerSet{ $key }."\n" )
    }
    close(SUMMARY_FILE);
}         

################################################################################
#
#			AFFYMETRIX HG-U133 Plus 2
#
################################################################################

%Probes =();
%Probes2rm =();
$probeUniqID = '';
%probesUniq = ();
%probesPerSet = ();
%nProbesPerSet = ();
%nXprobesPerSet = ();

if ( exists  $Arrays{ 'Affy_U133Plus2' } ) {
    
    $annotFile = $dir."/Affy_U133Plus2_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;    
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $probeUniqID = $info[0].$info[1];
        
        ##### Count number of probes in each probeset
        if ( !exists $probesUniq{ $probeUniqID } ) {
            
            $probesPerSet{ $info[0] }[0]++;
            
            ##### Count probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            if ( $info[5] ne "25=" or $info[6] != 1 ) {
                
                $probesPerSet{ $info[0] }[1]++;
                
            ##### Count probes mapping to other gene than remaining probes within the probeset
            } elsif ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                
                $probesPerSet{ $info[0] }[1]++;
            }
        }
        
        if ( $info[3] ne "-" ) {
            ##### Report probesets containing probes that map to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                    
                $Probes2rm{ $info[0] }=10;
            
            ##### Count probes, within a probeset, with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[5] ne "25=" or $info[6] != 1 ) {
            
                $Probes2rm{ $info[0] }++;
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[3];
            $Probes{ $info[0] }[2] = $info[4];
        }
        $probesUniq{ $probeUniqID } = $probeUniqID;
    }
    
    close(INFILE);
    
    $flatFile = $dir."/Affy_U133Plus2.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probesets
    while ( (my $key, my $value) = each %Probes ) {
        
        ##### Report only probesets with < 3 probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
        if ( exists  $Probes2rm{ $key } and $Probes2rm{ $key } > 2 ) {
            
            delete $Probes{ $key };
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);

################################################################################
#			Mapping summary for AFFYMETRIX HG-U133 Plus 2
################################################################################
    
    $summaryFile = $dir."/Affy_U133Plus2.perProbeset_summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    print ( SUMMARY_FILE "probeset\tprobes_number\tunreliable_probes_number\n" );
    
    ##### Write mapping summary per each probeset
    while ( (my $key, my $value) = each %probesPerSet ) {
        
        ##### Count probesets consisting of n probes
        $nProbesPerSet{ $probesPerSet{ $key }[0] }++;
        
        ##### Report probesets with > 1 unrealiable probes
        if ($probesPerSet{ $key }[1] ) {
        
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t".$probesPerSet{ $key }[1]."\n" );
            
            ##### Count probesets with n unrealiable probes
            $nXprobesPerSet{ $probesPerSet{ $key }[1] }++;
            
        ##### Report probesets with no unrealiable probes
        } else {
            
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t0\n" );
            
            ##### Count probesets with no unrealiable probes
            $nXprobesPerSet{ 0 }++;
        }
    }
    close(SUMMARY_FILE);
    
    
    $summaryFile = $dir."/Affy_U133Plus2.summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    
    ##### Report probesets with n probes
    print ( SUMMARY_FILE "probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nProbesPerSet ) {
            
            print( SUMMARY_FILE $key."\t".$nProbesPerSet{ $key }."\n" )
    }
    
    ##### Report probesets with n unrealiable probes
    print ( SUMMARY_FILE "\nunreliable_probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nXprobesPerSet ) {
        
        print( SUMMARY_FILE $key."\t".$nXprobesPerSet{ $key }."\n" )
    }
    close(SUMMARY_FILE);
}     

################################################################################
#
#			AFFYMETRIX HG-U133A 2
#
################################################################################

%Probes =();
%Probes2rm =();
$probeUniqID = '';
%probesUniq = ();
%probesPerSet = ();
%nProbesPerSet = ();
%nXprobesPerSet = ();

if ( exists  $Arrays{ 'Affy_U133A2' } ) {
    
    $annotFile = $dir."/Affy_U133A2_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;    
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $probeUniqID = $info[0].$info[1];
        
        ##### Count number of probes in each probeset
        if ( !exists $probesUniq{ $probeUniqID } ) {
            
            $probesPerSet{ $info[0] }[0]++;
            
            ##### Count probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            if ( $info[5] ne "25=" or $info[6] != 1 ) {
                
                $probesPerSet{ $info[0] }[1]++;
                
            ##### Count probes mapping to other gene than remaining probes within the probeset
            } elsif ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                
                $probesPerSet{ $info[0] }[1]++;
            }
        }
        
        if ( $info[3] ne "-" ) {
            ##### Report probesets containing probes that map to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                    
                $Probes2rm{ $info[0] }=10;
            
            ##### Count probes, within a probeset, with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[5] ne "25=" or $info[6] != 1 ) {
            
                $Probes2rm{ $info[0] }++;
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[3];
            $Probes{ $info[0] }[2] = $info[4];
        }
        $probesUniq{ $probeUniqID } = $probeUniqID;
    }
    
    close(INFILE);
    
    $flatFile = $dir."/Affy_U133A2.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probesets
    while ( (my $key, my $value) = each %Probes ) {
        
        ##### Report only probesets with < 3 probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
        if ( exists  $Probes2rm{ $key } and $Probes2rm{ $key } > 2 ) {
            
            delete $Probes{ $key };
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);

################################################################################
#			Mapping summary for AFFYMETRIX HG-U133A 2
################################################################################
    
    $summaryFile = $dir."/Affy_U133A2.perProbeset_summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    print ( SUMMARY_FILE "probeset\tprobes_number\tunreliable_probes_number\n" );
    
    ##### Write mapping summary per each probeset
    while ( (my $key, my $value) = each %probesPerSet ) {
        
        ##### Count probesets consisting of n probes
        $nProbesPerSet{ $probesPerSet{ $key }[0] }++;
        
        ##### Report probesets with > 1 unrealiable probes
        if ($probesPerSet{ $key }[1] ) {
        
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t".$probesPerSet{ $key }[1]."\n" );
            
            ##### Count probesets with n unrealiable probes
            $nXprobesPerSet{ $probesPerSet{ $key }[1] }++;
            
        ##### Report probesets with no unrealiable probes
        } else {
            
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t0\n" );
            
            ##### Count probesets with no unrealiable probes
            $nXprobesPerSet{ 0 }++;
        }
    }
    close(SUMMARY_FILE);
    
    
    $summaryFile = $dir."/Affy_U133A2.summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    
    ##### Report probesets with n probes
    print ( SUMMARY_FILE "probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nProbesPerSet ) {
            
            print( SUMMARY_FILE $key."\t".$nProbesPerSet{ $key }."\n" )
    }
    
    ##### Report probesets with n unrealiable probes
    print ( SUMMARY_FILE "\nunreliable_probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nXprobesPerSet ) {
        
        print( SUMMARY_FILE $key."\t".$nXprobesPerSet{ $key }."\n" )
    }
    close(SUMMARY_FILE);
}     

################################################################################
#
#			AFFYMETRIX HG-U133A
#
################################################################################

%Probes =();
%Probes2rm =();
$probeUniqID = '';
%probesUniq = ();
%probesPerSet = ();
%nProbesPerSet = ();
%nXprobesPerSet = ();

if ( exists  $Arrays{ 'Affy_U133A' } ) {
    
    $annotFile = $dir."/Affy_U133A_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;    
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $probeUniqID = $info[0].$info[1];
        
        ##### Count number of probes in each probeset
        if ( !exists $probesUniq{ $probeUniqID } ) {
            
            $probesPerSet{ $info[0] }[0]++;
            
            ##### Count probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            if ( $info[5] ne "25=" or $info[6] != 1 ) {
                
                $probesPerSet{ $info[0] }[1]++;
                
            ##### Count probes mapping to other gene than remaining probes within the probeset
            } elsif ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                
                $probesPerSet{ $info[0] }[1]++;
            }
        }
        
        if ( $info[3] ne "-" ) {
            ##### Report probesets containing probes that map to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                    
                $Probes2rm{ $info[0] }=10;
            
            ##### Count probes, within a probeset, with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[5] ne "25=" or $info[6] != 1 ) {
            
                $Probes2rm{ $info[0] }++;
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[3];
            $Probes{ $info[0] }[2] = $info[4];
        }
        $probesUniq{ $probeUniqID } = $probeUniqID;
    }
    
    close(INFILE);
    
    $flatFile = $dir."/Affy_U133A.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probesets
    while ( (my $key, my $value) = each %Probes ) {
        
        ##### Report only probesets with < 3 probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
        if ( exists  $Probes2rm{ $key } and $Probes2rm{ $key } > 2 ) {
            
            delete $Probes{ $key };
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);

################################################################################
#			Mapping summary for AFFYMETRIX HG-U133A
################################################################################
    
    $summaryFile = $dir."/Affy_U133A.perProbeset_summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    print ( SUMMARY_FILE "probeset\tprobes_number\tunreliable_probes_number\n" );
    
    ##### Write mapping summary per each probeset
    while ( (my $key, my $value) = each %probesPerSet ) {
        
        ##### Count probesets consisting of n probes
        $nProbesPerSet{ $probesPerSet{ $key }[0] }++;
        
        ##### Report probesets with > 1 unrealiable probes
        if ($probesPerSet{ $key }[1] ) {
        
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t".$probesPerSet{ $key }[1]."\n" );
            
            ##### Count probesets with n unrealiable probes
            $nXprobesPerSet{ $probesPerSet{ $key }[1] }++;
            
        ##### Report probesets with no unrealiable probes
        } else {
            
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t0\n" );
            
            ##### Count probesets with no unrealiable probes
            $nXprobesPerSet{ 0 }++;
        }
    }
    close(SUMMARY_FILE);
    
    
    $summaryFile = $dir."/Affy_U133A.summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    
    ##### Report probesets with n probes
    print ( SUMMARY_FILE "probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nProbesPerSet ) {
            
            print( SUMMARY_FILE $key."\t".$nProbesPerSet{ $key }."\n" )
    }
    
    ##### Report probesets with n unrealiable probes
    print ( SUMMARY_FILE "\nunreliable_probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nXprobesPerSet ) {
        
        print( SUMMARY_FILE $key."\t".$nXprobesPerSet{ $key }."\n" )
    }
    close(SUMMARY_FILE);
}     

################################################################################
#
#			AFFYMETRIX HG-U133B
#
################################################################################

%Probes =();
%Probes2rm =();
$probeUniqID = '';
%probesUniq = ();
%probesPerSet = ();
%nProbesPerSet = ();
%nXprobesPerSet = ();

if ( exists  $Arrays{ 'Affy_U133B' } ) {
    
    $annotFile = $dir."/Affy_U133B_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;    
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $probeUniqID = $info[0].$info[1];
        
        ##### Count number of probes in each probeset
        if ( !exists $probesUniq{ $probeUniqID } ) {
            
            $probesPerSet{ $info[0] }[0]++;
            
            ##### Count probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            if ( $info[5] ne "25=" or $info[6] != 1 ) {
                
                $probesPerSet{ $info[0] }[1]++;
                
            ##### Count probes mapping to other gene than remaining probes within the probeset
            } elsif ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                
                $probesPerSet{ $info[0] }[1]++;
            }
        }
        
        if ( $info[3] ne "-" ) {
            ##### Report probesets containing probes that map to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                    
                $Probes2rm{ $info[0] }=10;
            
            ##### Count probes, within a probeset, with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[5] ne "25=" or $info[6] != 1 ) {
            
                $Probes2rm{ $info[0] }++;
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[3];
            $Probes{ $info[0] }[2] = $info[4];
        }
        $probesUniq{ $probeUniqID } = $probeUniqID;
    }
    
    close(INFILE);
    
    $flatFile = $dir."/Affy_U133B.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probesets
    while ( (my $key, my $value) = each %Probes ) {
        
        ##### Report only probesets with < 3 probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
        if ( exists  $Probes2rm{ $key } and $Probes2rm{ $key } > 2 ) {
            
            delete $Probes{ $key };
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);

################################################################################
#			Mapping summary for AFFYMETRIX HG-U133B
################################################################################
    
    $summaryFile = $dir."/Affy_U133B.perProbeset_summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    print ( SUMMARY_FILE "probeset\tprobes_number\tunreliable_probes_number\n" );
    
    ##### Write mapping summary per each probeset
    while ( (my $key, my $value) = each %probesPerSet ) {
        
        ##### Count probesets consisting of n probes
        $nProbesPerSet{ $probesPerSet{ $key }[0] }++;
        
        ##### Report probesets with > 1 unrealiable probes
        if ($probesPerSet{ $key }[1] ) {
        
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t".$probesPerSet{ $key }[1]."\n" );
            
            ##### Count probesets with n unrealiable probes
            $nXprobesPerSet{ $probesPerSet{ $key }[1] }++;
            
        ##### Report probesets with no unrealiable probes
        } else {
            
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t0\n" );
            
            ##### Count probesets with no unrealiable probes
            $nXprobesPerSet{ 0 }++;
        }
    }
    close(SUMMARY_FILE);
    
    
    $summaryFile = $dir."/Affy_U133B.summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    
    ##### Report probesets with n probes
    print ( SUMMARY_FILE "probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nProbesPerSet ) {
            
            print( SUMMARY_FILE $key."\t".$nProbesPerSet{ $key }."\n" )
    }
    
    ##### Report probesets with n unrealiable probes
    print ( SUMMARY_FILE "\nunreliable_probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nXprobesPerSet ) {
        
        print( SUMMARY_FILE $key."\t".$nXprobesPerSet{ $key }."\n" )
    }
    close(SUMMARY_FILE);
}     

################################################################################
#
#			AFFYMETRIX HG-U95Av2
#
################################################################################

%Probes =();
%Probes2rm =();
$probeUniqID = '';
%probesUniq = ();
%probesPerSet = ();
%nProbesPerSet = ();
%nXprobesPerSet = ();

if ( exists  $Arrays{ 'Affy_U95Av2' } ) {
    
    $annotFile = $dir."/Affy_U95Av2_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;    
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $probeUniqID = $info[0].$info[1];
        
        ##### Count number of probes in each probeset
        if ( !exists $probesUniq{ $probeUniqID } ) {
            
            $probesPerSet{ $info[0] }[0]++;
            
            ##### Count probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            if ( $info[5] ne "25=" or $info[6] != 1 ) {
                
                $probesPerSet{ $info[0] }[1]++;
                
            ##### Count probes mapping to other gene than remaining probes within the probeset
            } elsif ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                
                $probesPerSet{ $info[0] }[1]++;
            }
        }
        
        if ( $info[3] ne "-" ) {
            ##### Report probesets containing probes that map to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                    
                $Probes2rm{ $info[0] }=10;
            
            ##### Count probes, within a probeset, with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[5] ne "25=" or $info[6] != 1 ) {
            
                $Probes2rm{ $info[0] }++;
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[3];
            $Probes{ $info[0] }[2] = $info[4];
        }
        $probesUniq{ $probeUniqID } = $probeUniqID;
    }
    
    close(INFILE);
    
    $flatFile = $dir."/Affy_U95Av2.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probesets
    while ( (my $key, my $value) = each %Probes ) {
        
        ##### Report only probesets with < 3 probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
        if ( exists  $Probes2rm{ $key } and $Probes2rm{ $key } > 2 ) {
            
            delete $Probes{ $key };
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);

################################################################################
#			Mapping summary for AFFYMETRIX HG-U95Av2
################################################################################
    
    $summaryFile = $dir."/Affy_U95Av2.perProbeset_summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    print ( SUMMARY_FILE "probeset\tprobes_number\tunreliable_probes_number\n" );
    
    ##### Write mapping summary per each probeset
    while ( (my $key, my $value) = each %probesPerSet ) {
        
        ##### Count probesets consisting of n probes
        $nProbesPerSet{ $probesPerSet{ $key }[0] }++;
        
        ##### Report probesets with > 1 unrealiable probes
        if ($probesPerSet{ $key }[1] ) {
        
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t".$probesPerSet{ $key }[1]."\n" );
            
            ##### Count probesets with n unrealiable probes
            $nXprobesPerSet{ $probesPerSet{ $key }[1] }++;
            
        ##### Report probesets with no unrealiable probes
        } else {
            
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t0\n" );
            
            ##### Count probesets with no unrealiable probes
            $nXprobesPerSet{ 0 }++;
        }
    }
    close(SUMMARY_FILE);
    
    
    $summaryFile = $dir."/Affy_U95Av2.summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    
    ##### Report probesets with n probes
    print ( SUMMARY_FILE "probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nProbesPerSet ) {
            
            print( SUMMARY_FILE $key."\t".$nProbesPerSet{ $key }."\n" )
    }
    
    ##### Report probesets with n unrealiable probes
    print ( SUMMARY_FILE "\nunreliable_probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nXprobesPerSet ) {
        
        print( SUMMARY_FILE $key."\t".$nXprobesPerSet{ $key }."\n" )
    }
    close(SUMMARY_FILE);
}     

################################################################################
#
#			AFFYMETRIX HG-U95B
#
################################################################################

%Probes =();
%Probes2rm =();
$probeUniqID = '';
%probesUniq = ();
%probesPerSet = ();
%nProbesPerSet = ();
%nXprobesPerSet = ();

if ( exists  $Arrays{ 'Affy_U95B' } ) {
    
    $annotFile = $dir."/Affy_U95B_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;    
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $probeUniqID = $info[0].$info[1];
        
        ##### Count number of probes in each probeset
        if ( !exists $probesUniq{ $probeUniqID } ) {
            
            $probesPerSet{ $info[0] }[0]++;
            
            ##### Count probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            if ( $info[5] ne "25=" or $info[6] != 1 ) {
                
                $probesPerSet{ $info[0] }[1]++;
                
            ##### Count probes mapping to other gene than remaining probes within the probeset
            } elsif ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                
                $probesPerSet{ $info[0] }[1]++;
            }
        }
        
        if ( $info[3] ne "-" ) {
            ##### Report probesets containing probes that map to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                    
                $Probes2rm{ $info[0] }=10;
            
            ##### Count probes, within a probeset, with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[5] ne "25=" or $info[6] != 1 ) {
            
                $Probes2rm{ $info[0] }++;
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[3];
            $Probes{ $info[0] }[2] = $info[4];
        }
        $probesUniq{ $probeUniqID } = $probeUniqID;
    }
    
    close(INFILE);
    
    $flatFile = $dir."/Affy_U95B.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probesets
    while ( (my $key, my $value) = each %Probes ) {
        
        ##### Report only probesets with < 3 probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
        if ( exists  $Probes2rm{ $key } and $Probes2rm{ $key } > 2 ) {
            
            delete $Probes{ $key };
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);

################################################################################
#			Mapping summary for AFFYMETRIX HG-U95B
################################################################################
    
    $summaryFile = $dir."/Affy_U95B.perProbeset_summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    print ( SUMMARY_FILE "probeset\tprobes_number\tunreliable_probes_number\n" );
    
    ##### Write mapping summary per each probeset
    while ( (my $key, my $value) = each %probesPerSet ) {
        
        ##### Count probesets consisting of n probes
        $nProbesPerSet{ $probesPerSet{ $key }[0] }++;
        
        ##### Report probesets with > 1 unrealiable probes
        if ($probesPerSet{ $key }[1] ) {
        
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t".$probesPerSet{ $key }[1]."\n" );
            
            ##### Count probesets with n unrealiable probes
            $nXprobesPerSet{ $probesPerSet{ $key }[1] }++;
            
        ##### Report probesets with no unrealiable probes
        } else {
            
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t0\n" );
            
            ##### Count probesets with no unrealiable probes
            $nXprobesPerSet{ 0 }++;
        }
    }
    close(SUMMARY_FILE);
    
    
    $summaryFile = $dir."/Affy_U95B.summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    
    ##### Report probesets with n probes
    print ( SUMMARY_FILE "probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nProbesPerSet ) {
            
            print( SUMMARY_FILE $key."\t".$nProbesPerSet{ $key }."\n" )
    }
    
    ##### Report probesets with n unrealiable probes
    print ( SUMMARY_FILE "\nunreliable_probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nXprobesPerSet ) {
        
        print( SUMMARY_FILE $key."\t".$nXprobesPerSet{ $key }."\n" )
    }
    close(SUMMARY_FILE);
}     

################################################################################
#
#			AFFYMETRIX HG-U95C
#
################################################################################

%Probes =();
%Probes2rm =();
$probeUniqID = '';
%probesUniq = ();
%probesPerSet = ();
%nProbesPerSet = ();
%nXprobesPerSet = ();

if ( exists  $Arrays{ 'Affy_U95C' } ) {
    
    $annotFile = $dir."/Affy_U95C_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;    
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);
        
        $probeUniqID = $info[0].$info[1];
        
        ##### Count number of probes in each probeset
        if ( !exists $probesUniq{ $probeUniqID } ) {
            
            $probesPerSet{ $info[0] }[0]++;
            
            ##### Count probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            if ( $info[5] ne "25=" or $info[6] != 1 ) {
                
                $probesPerSet{ $info[0] }[1]++;
                
            ##### Count probes mapping to other gene than remaining probes within the probeset
            } elsif ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                
                $probesPerSet{ $info[0] }[1]++;
            }
        }
        
        if ( $info[3] ne "-" ) {
            ##### Report probesets containing probes that map to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[3] ) {
                    
                $Probes2rm{ $info[0] }=10;
            
            ##### Count probes, within a probeset, with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[5] ne "25=" or $info[6] != 1 ) {
            
                $Probes2rm{ $info[0] }++;
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[3];
            $Probes{ $info[0] }[2] = $info[4];
        }
        $probesUniq{ $probeUniqID } = $probeUniqID;
    }
    
    close(INFILE);
    
    $flatFile = $dir."/Affy_U95C.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probesets
    while ( (my $key, my $value) = each %Probes ) {
        
        ##### Report only probesets with < 3 probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
        if ( exists  $Probes2rm{ $key } and $Probes2rm{ $key } > 2 ) {
            
            delete $Probes{ $key };
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);

################################################################################
#			Mapping summary for AFFYMETRIX HG-U95C
################################################################################
    
    $summaryFile = $dir."/Affy_U95C.perProbeset_summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    print ( SUMMARY_FILE "probeset\tprobes_number\tunreliable_probes_number\n" );
    
    ##### Write mapping summary per each probeset
    while ( (my $key, my $value) = each %probesPerSet ) {
        
        ##### Count probesets consisting of n probes
        $nProbesPerSet{ $probesPerSet{ $key }[0] }++;
        
        ##### Report probesets with > 1 unrealiable probes
        if ($probesPerSet{ $key }[1] ) {
        
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t".$probesPerSet{ $key }[1]."\n" );
            
            ##### Count probesets with n unrealiable probes
            $nXprobesPerSet{ $probesPerSet{ $key }[1] }++;
            
        ##### Report probesets with no unrealiable probes
        } else {
            
            print( SUMMARY_FILE $key."\t".$probesPerSet{ $key }[0]."\t0\n" );
            
            ##### Count probesets with no unrealiable probes
            $nXprobesPerSet{ 0 }++;
        }
    }
    close(SUMMARY_FILE);
    
    
    $summaryFile = $dir."/Affy_U95C.summary.txt";
    
    open ( SUMMARY_FILE, ">$summaryFile" );
    
    ##### Report probesets with n probes
    print ( SUMMARY_FILE "probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nProbesPerSet ) {
            
            print( SUMMARY_FILE $key."\t".$nProbesPerSet{ $key }."\n" )
    }
    
    ##### Report probesets with n unrealiable probes
    print ( SUMMARY_FILE "\nunreliable_probes_number\tfrequency\n" );
    
    while ( (my $key, my $value) = each %nXprobesPerSet ) {
        
        print( SUMMARY_FILE $key."\t".$nXprobesPerSet{ $key }."\n" )
    }
    close(SUMMARY_FILE);
}     

################################################################################
#
#			ILLUMINA HumanHT-12
#
################################################################################

%Probes =();
%Probes2rm =();

if ( exists  $Arrays{ 'Illum_HT_12_V3' } ) {
    
    $annotFile = $dir."/Illum_HT_12_V3_annot.txt";
    
    ##### Read the annotation file
    open (INFILE, $annotFile) or die $!;
    ##### Skip the first line
    <INFILE>;
    
    while (my $record = <INFILE>) {
        
        chomp $record;
        
        my @info = split(/\t/, $record);

        if ( $info[2] ne "-" ) {
            
            ##### Report probes mapping to various genes
            if ( exists $Probes{ $info[0] } and $Probes{ $info[0] }[1] ne $info[2] ) {
                
                $Probes2rm{ $info[0] } = $info[0];
                
            ##### Report probes with mismatch/insterion/deletion or mapping to a region annotated with none or more than one gene
            } elsif ( $info[4] ne "50=" or $info[5] != 1 ) {
                
                $Probes2rm{ $info[0] } = $info[0];
            }
            $Probes{ $info[0] }[0] = $info[0];
            $Probes{ $info[0] }[1] = $info[2];
            $Probes{ $info[0] }[2] = $info[3];
        }
    }
    close(INFILE);
    
    $flatFile = $dir."/Illum_HT_12_V3.flat";
    
    open ( FLAT_FILE, ">$flatFile" );
    print ( FLAT_FILE "probe\tgene_idD\tgene_short_name\n" );
    
    ##### Go through all probes
    while ( (my $key, my $value) = each %Probes ) {
        
        my @info = split(/\t/, $value);
        
        if ( exists  $Probes2rm{ $key } ) {
            
            delete $Probes{ $key };
            
            ##### Report only probes with perfect match
        } else {
            print( FLAT_FILE $Probes{ $key }[0]."\t".$Probes{ $key }[1]."\t".$Probes{ $key }[2]."\n" )
        }
    }
    close(FLAT_FILE);
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
        
    usage: $0 -i info_file
        
    Input data
    -i info_file:  File name with full path that lists the datasets to be analysed.
                   The file is expected to include the following columns:
                   (1) dataset name; (2) target file name with full path; (3) data
                   type (raw or processed) and (4) platform used

EOS
exit;
}
