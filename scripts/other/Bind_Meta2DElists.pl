#!/usr/bin/perl -w
################################################################################
#
#   File name: Bind_Meta2DElists.pl
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
#	Description: Script combining meta-analysis results with differential expression analyses results performed on individual studies.
#
#	Command line use example: ./Bind_Meta2DElists.pl -o /scratch/jack/data/PhD/Transcriptomics_project -g /scratch/jack/data/PhD/Transcriptomics_project/PCa_assocaited_genes_DDPC-Schoenborn-COSMIC.txt -e /scratch/jack/data/PhD/Tissue_components_analysis/stroma_genes_to_exclude.txt -m /scratch/jack/data/PhD/Tissue_components_analysis/Stuart_tumour_genes.txt,/scratch/jack/data/PhD/Tissue_components_analysis/Stuart_T_greater_S_genes.txt,/scratch/jack/data/PhD/Tissue_components_analysis/Stuart_B_and_T_genes.txt,/scratch/jack/data/PhD/Tissue_components_analysis/ESTIMATE_stromal_genes.txt -c Tumour_Stuart_T,Tumour_Stuart_T_greater_S,Tumour_Stuart_B_and_T,Tumour_GSE20758
#
#	-o my_project:      Project workspace. This is the directory to which all the results will be written
#	-g known_genes (optional):  List of known genes associated with studied phenotype. The first column is expected to list the Ensembl gene IDs
#	-e exclude_genes (optional):   List of genes to be excluded. The first column is expected to list the Ensembl gene IDs. These genes will be reported in separete files
#	-m mark_genes (optional):   List(s) of genes to be marked in the output files. The first column is expected to list the Ensembl gene IDs. Multiple gene lists should be separated by comma (spaces are not allowed)
#	-c mark_char (optional):    Character to be used to mark uder-defined genes. Used if the '-c' option is used. Multiple characters for corresponding gene lists should be separated by comma
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
my $studyNo;
my $datasetsList;
my $metaResults;
my $comparison;
my $DElist;
my %geneList=();
my $metaHeader;

my $knownGenes=0;
my %knownGenesList=();
my $genes2excl;
my %genes2exclList=();
my $genes2mark;
my %genes2markList=();
my $char;
my @chars;

while ($arg = shift) {
    
    if ($arg =~ /^-o$/) {
        $Project = shift;
    } elsif ($arg =~ /^-g$/) {
        $knownGenes = shift;
    } elsif ($arg =~ /^-e$/) {
        $genes2excl = shift;
    } elsif ($arg =~ /^-m$/) {
        $genes2mark = shift;
    } elsif ($arg =~ /^-c$/) {
        $char = shift;
    } elsif ($arg =~ /help/) {
        usage ();
    }
}

if ( $Project ) {
    
    ##### Get list of known genes associated with studied phenotype (if provided)
    if ( $knownGenes ) {
    	    
    	     open (KNOWN_GENE_LIST, $knownGenes) or die $!;
    	     
    	     while (my $record = <KNOWN_GENE_LIST>) {
                    
    	     	     chomp $record;
    	     	     
    	     	     my @genesInfo = split(/\s/, $record);
            	    	    	    	    	    
    	     	     $knownGenesList{ $genesInfo[0] } = $genesInfo[0];
    	     }
    	     print( "\n\nList of ".keys(%knownGenesList)." known genes provided\n\n");
    }
    
    ##### Get list of genes to be excluded (if provided). Duplicated gene IDs will be excluded
    if ( $genes2excl ) {
        
        open (ENE2EXCL_LIST, $genes2excl) or die $!;
        
        while (my $record = <ENE2EXCL_LIST>) {
            
            chomp $record;
            
            my @genesInfo = split(/\s/, $record);
            
            $genes2exclList{ $genesInfo[0] } = $genesInfo[0];
        }
        print( "\n\nList of ".keys(%genes2exclList)." genes to be excluded\n\n");
    }

    ##### Get list of genes to be marked (if provided). Duplicated gene IDs will be excluded
    if ( $genes2mark ) {
        
        ##### Go through all files with genes to mark
        my @filesInfo = split(/\,/, $genes2mark);
        
        ##### Get characters to be used for genes marking (if provided)
        if ( $char ) {
            
            @chars = split(/\,/, $char);
            
        ##### Use numeric marks for genes in case user didn't provided the list of marks
        } else {
            
            @chars = ( 1..scalar(@filesInfo) );
        }

        my $marks=0;
    
        foreach my $file (@filesInfo) {
            
            open (GENE2MARK_LIST, $file) or die $!;

            while (my $record = <GENE2MARK_LIST>) {
            
                chomp $record;
            
                my @genesInfo = split(/\s/, $record);
            
                $genes2markList{ $genesInfo[0] } = $chars[$marks];
            }
            
            $marks++;
        }
        
        print( "\n\nList of ".keys(%genes2markList)." genes to be marked\n\n");
    }
    

    ##### Open a directory to project workspace
    opendir (PROJECT_DIR, $Project) or die $!;
    
    ##### Look for files with meta-analysis results
    while (my $file = readdir(PROJECT_DIR)) {
        
        %geneList=();
        
    	if ( $knownGenes ) {
    		
    		$datasetsList="\tKnown";
    	} else {
    		$datasetsList='';
    	}
    	    
        if ( $file =~ m/Meta_(.*?vs.*?)\.txt/ and $file !~ m/Meta_.*?vs.*?\_GO.*?\.txt/ and $file !~ m/Meta_.*?vs.*?\_DEbound.*?\.txt/ ) {
        
            $metaResults = join('/', $Project, $file );
            
            $comparison = $1;
            
            print( "\nProcessing meta-analsis results for $comparison\n\n");

            ##### Get list of datasets to be analysed
            if (-e $Project."/GenExpression_InputFiles.txt") {
        
            	    open (DATASETS_LIST, $Project."/GenExpression_InputFiles.txt") or die $!;

            	    ##### Skip the first line
            	    <DATASETS_LIST>;
            	    
            	    $studyNo=0;

            	    while (my $dataset = <DATASETS_LIST>) {
    
            	    	    chomp $dataset;
    
            	    	    my @datasetInfo = split(/\s/, $dataset);

            	    	    @dir = split('/', $datasetInfo[1]);
            	    	    $dir = join('/', @dir[0 .. $#dir-1]);
            	    	    $studyID = $datasetInfo[0];
            	    	   
            	    	    opendir (PROJECT, $Project) or die $!;
            	    	    
            	    	    while (my $Comb_file = readdir(PROJECT)) {
            	    	    	    	    
            	    	    	    if ( $Comb_file =~ m/Comb\_$studyID\_$1\_topTable\.txt/ ) {
                                        
                                            print("Reading statistics for gene list $Comb_file\n\n");
            	    	    	    	    
            	    	    	    	    $datasetsList .= "\t".$studyID." (log2FC)"."\t".$studyID." (adj p-value)";
            	    	    	    	    
            	    	    	    	    ##### Store DE gene list in a hash
            	    	    	    	    $DElist = join('/', $Project, $Comb_file);
                
            	    	    	    	    open (INFILE, $DElist) or die $!;
                
            	    	    	    	    while (my $record = <INFILE>) {
                    
            	    	    	    	    	    chomp $record;
                                                
            	    	    	    	    	    my @info = split(/\t/, $record);
                                                
            	    	    	    	    	    $geneList{ $info[0] }[$studyNo] = $info[8]."\t".$info[11];
            	    	    	    	    }
            	    	    	    	    $studyNo++;
            	    	    	    }
            	    	    }
            	    }
            	    close(DATASETS_LIST);
        
            } else {
            	    print( "\n".$Project."/GenExpression_InputFiles.txt does not exist!\n\n");
            	    exit;
            }

            print( "Binding results...\n\n");
            
            open ( META_FILE, $metaResults) or die $!;
            
            $metaResults =~ s/\.txt/_DEbound\.txt/;
            open ( BOUND_FILE, ">$metaResults" );
                
            ##### Write the header
            $metaHeader = scalar <META_FILE>;
            chomp $metaHeader;

            print( BOUND_FILE $metaHeader."\t".$datasetsList."\n");
            
            
            ##### Report results for gene to be exluded (if provided)
            if ( $genes2excl ) {
                
                $metaResults =~ s/\.txt/_excluded\.txt/;
                open ( GENE2EXCL_BOUND_FILE, ">$metaResults" );
                
                ##### Write the header
                print( GENE2EXCL_BOUND_FILE $metaHeader."\t".$datasetsList."\n");
            }

            while (my $metaFile = <META_FILE>) {
            	    
            	    chomp $metaFile;
                
            	    ##### Report statistics for meta-analysis
            	    my @metaInfo = split("\t", $metaFile);
                
                    my $gene = $metaInfo[0];

                    ##### Check which genes need to be MARKED
                    if ( exists $genes2markList{ $gene } ) {
                    
                        $metaInfo[0] = $genes2markList{ $metaInfo[0] }."_".$metaInfo[0];
                        
                        $metaFile = join("\t", @metaInfo);
                    }
                
                    ##### Check which genes need to be EXCLUDED
                    if ( exists $genes2exclList{ $gene } ) {
                        
                        print ( GENE2EXCL_BOUND_FILE $metaFile."\t" );
                        
                        ##### Tag known genes
                        if ( $knownGenes ) {
            	    	    
                            if ( exists $knownGenesList{ $gene } ) {
                                print ( GENE2EXCL_BOUND_FILE "\tyes" );
            	    	    } else {
                                print ( GENE2EXCL_BOUND_FILE "\tNo" );
            	    	    }
                        }
                        
                        ##### Report statistics for individual studies
                        for (my $i=0; $i<$studyNo; $i++) {
            	    	    
            	    	    if ( ${ $geneList{$gene} }[$i] ) {
                                print ( GENE2EXCL_BOUND_FILE "\t".${ $geneList{ $gene } }[$i] );
            	    	    } else {
                                print ( GENE2EXCL_BOUND_FILE "\tNA\tNA" );
            	    	    }
                        }
                        print ( GENE2EXCL_BOUND_FILE "\n" );
                        
                    } else {
                        
                        print ( BOUND_FILE $metaFile."\t" );
                        
                        ##### Tag known genes
                        if ( $knownGenes ) {
            	    	    
                            if ( exists $knownGenesList{ $gene } ) {
                                print ( BOUND_FILE "\tyes" );
            	    	    } else {
                                print ( BOUND_FILE "\tNo" );
            	    	    }
                        }
                        
                        ##### Report statistics for individual studies
                        for (my $i=0; $i<$studyNo; $i++) {
            	    	    
            	    	    if ( ${ $geneList{$gene} }[$i] ) {
                                print ( BOUND_FILE "\t".${ $geneList{ $gene } }[$i] );
            	    	    } else {
                                print ( BOUND_FILE "\tNA\tNA" );
            	    	    }
                        }
                        print ( BOUND_FILE "\n" );
                        
                    }
            }
            close(BOUND_FILE);
            close(META_FILE);
            
            if ( $genes2excl ) {
                close(GENE2EXCL_BOUND_FILE);
            }
        }
    }
    
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
    
    usage: $0 -o my_project [optional] -g known_genes -e exclude_genes -m mark_genes -c mark_char
    
    Output data
    -o my_project:   Project workspace. This is the directory to which all the results
                     will be written
                     
    -g known_genes:  List of known genes associated with studied phenotype.
    		     The first column is expected to list the Ensembl gene IDs
    
    -e exclude_genes:  List of genes to be excluded. The first column is
                 expected to list the Ensembl gene IDs. These genes will be
                 reported in separete files
    
    -m mark_genes:  List(s) of genes to be marked in the output files. The
                 first column is expected to list the Ensembl gene IDs.
                 Multiple gene lists should be separated by comma (spaces
                 are not allowed)
    
    -c mark_char:  Character to be used to mark uder-defined genes. Used if
                 the '-c' option is used. Multiple characters for
                 corresponding gene lists should be separated by comma
   
EOS
exit;
}
