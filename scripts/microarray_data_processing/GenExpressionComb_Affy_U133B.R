################################################################################
#
#   File name: GenExpressionComb_Affy_U133B.R
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
#	Description: Pipeline combining datasets produced using Affymetrix Human Genome U133B Array. The script reads in three a priori generated files: (1) list of datasets to be analysed; (2) file produced with QC_Affy_U133B.R which lists outliers identified with ArrayOutliers function and (3) '.flat' file produced with ArrayAnnot.pl which lists reliable probes annotated with Ensembl gene IDs and gene symbol. Pipeline executes 'multiGene2ProbeFilter.R' to remove multi gene to probe mappings
#
#	Command line use example: R --file=./GenExpressionComb_Affy_U133B.R --args "/scratch/jack/data/PhD/Transcriptomics_project"
#
#	First arg:     Project workspace. This is the directory to which cross-platform analysis results will be written
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Functions
#===============================================================================

##### Prepare object to write into a file
prepare2write <- function (x) {
	
	x2write <- cbind(rownames(x), x)
    colnames(x2write) <- c("",colnames(x))
	return(x2write)
}

##### Assign colours to analysed datasets
getDatasetsColours <- function(datasets) {
    
    ##### Predefined selection of colours for datasets
    datasets.colours <- c("bisque","orange","firebrick","lightslategrey","darkseagreen","darkcyan","dodgerblue")
    
    f.datasets <- factor(datasets)
    vec.datasets <- datasets.colours[1:length(levels(f.datasets))]
    datasets.colour <- rep(0,length(f.datasets))
    for(i in 1:length(f.datasets))
    datasets.colour[i] <- vec.datasets[ f.datasets[i]==levels(f.datasets)]

    return( list(vec.datasets, datasets.colour) )
}

#===============================================================================
#    Load libraries
#===============================================================================

library(affy)
library(gcrma)
library(simpleaffy)
library(sva)

source("MultiGene2ProbeFilter.R")

#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

ProjectDir = args[4]
probesFlatFile = paste(paste(ProjectDir, collapse = "/"),"Affy_U133B.flat",sep="/")

##### Read file with datasets information
DatasetInput=read.table(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE,row.names=1)
DatasetInput=DatasetInput[DatasetInput[,"Platform"] %in% "Affy_U133B",]

##### Change working directory to the project workspace
setwd(ProjectDir)

#===============================================================================
#    Load data from affymetrix '.CEL' files
#===============================================================================

fileInfo = strsplit(DatasetInput[,"TargetFile"], split='/', fixed=TRUE)

targetFile <- read.table(DatasetInput[1,"TargetFile"],sep="\t",as.is=TRUE,header=TRUE,row.names=1)[,c(1:3)]
targetFile <- cbind(targetFile,rownames(DatasetInput[1,]))
colnames(targetFile)[ncol(targetFile)] <- "Dataset"
targetFile[,"FileName"] <- paste(paste(fileInfo[[1]][1:length(fileInfo[[1]])-1], collapse="/"), targetFile[,"FileName"], sep="/")

if ( nrow(DatasetInput) > 1 ) {
    for ( i in 2:nrow(DatasetInput) ) {
    
        targetFileTmp <- read.table(DatasetInput[i,"TargetFile"],sep="\t",as.is=TRUE,header=TRUE,row.names=1)[,c(1:3)]
        targetFileTmp <- cbind(targetFileTmp,rownames(DatasetInput[i,]))
        colnames(targetFileTmp)[ncol(targetFileTmp)] <- "Dataset"
        targetFileTmp[,"FileName"] <- paste(paste(fileInfo[[i]][1:length(fileInfo[[1]])-1], collapse="/"), targetFileTmp[,"FileName"], sep="/")
    
        ##### Deal with replicates
        if ( any(!is.na(targetFile[,"Replicate"])) ) {
            maxRep <- max(targetFile[!is.na(targetFile[,"Replicate"]),"Replicate"])
            targetFileTmp[,"Replicate"] <- targetFileTmp[,"Replicate"] + maxRep
        }
        targetFile <- rbind(targetFile, targetFileTmp)
    }
}

cat("Loading datasets produced using Affymetrix GeneChip Human Genome U133B Array into an AffyBatch object...\n")

##### Remove samples identified as outliers
studyIDs <- rownames(DatasetInput)

outlierList <- NULL
j<-1
for (i in 1:length(studyIDs)) {
	outliersInfo = read.table(paste(ProjectDir,"/outliers_", studyIDs[i],".txt",sep=""),sep="\t",as.is=TRUE,header=TRUE,row.names=1)
	
	if ( !is.na(outliersInfo[,2]) ) {
		outliers =unlist(strsplit(outliersInfo[,2], split=',', fixed=TRUE))
		for (k in 1:length(outliers)) {
			cat(paste("Removing sample:", outliers[k], "\n", sep=" "))
			outlierList[j] <- outliers[k]
			j<-j+1
		}
	}
}
targetFile=targetFile[setdiff(rownames(targetFile),as.vector(as.list(outlierList))),]

##### Keep only reliable probes based on cigar strings derived from Ensembl microarray probe mapping
Affy_U133B_probes <- read.table(probesFlatFile,sep="\t",as.is=TRUE,header=TRUE,row.names=1)

##### Select an affybatch probeset subset prior to normalization, based on the following post https://stat.ethz.ch/pipermail/bioconductor/2009-January/026046.html
library(hgu133bcdf)

##### Report probesets to remove
difference=setdiff(ls(hgu133bcdf),rownames(Affy_U133B_probes))
##### Remove these genes from the cdf environment -> the cdf will have only the interesting genes
rm(list=difference, envir=hgu133bcdf)

dat<-ReadAffy(filenames = targetFile[,"FileName"], sampleNames = rownames(targetFile), phenoData=as(as.data.frame(targetFile[,"Target"]), "AnnotatedDataFrame") )
colnames(pData(dat)) <- "Target"

#===============================================================================
#    Normalisation with GC-RMA
#===============================================================================

cat("Robust Multiarray Average normalisation with GC-content background correction.\n")
Affy_U133B_gcrma <- NULL
Affy_U133B_gcrma <- gcrma(dat)

#save.image("Affy_U133B_gcrma.R")
#load("Affy_U133B_gcrma.R")

Affy_U133B_data <- exprs(Affy_U133B_gcrma)

dim(Affy_U133B_data)

#===============================================================================
#    Boxplot on normalised data
#===============================================================================

Affy_U133B_datasets <- targetFile$Dataset

Affy_U133B_datasets_No <- max(as.numeric(factor(Affy_U133B_datasets)))

Affy_U133B_datasets.colour <- getDatasetsColours(Affy_U133B_datasets)


pdf("Affy_U133B_boxplot.pdf", pointsize = 8 ,width = 0.2*length(targetFile$FileName), height = 6)
par(mar=c(13, 4, 3, 2))
boxplot(dat,col=Affy_U133B_datasets.colour[[2]], las = 2) # Generates boxplot of non-normalized intensity values.
boxplot(data.frame(Affy_U133B_data),col=Affy_U133B_datasets.colour[[2]], main="Normalised GC-RMA data", las = 2) # Generates boxplot of normalized log intensity values.
dev.off()

#===============================================================================
#    Remove multi gene to probeset mappings
#===============================================================================

#####  Dataset is presumed to contain non-unique gene names, so write the gene names into the first column rather then into the rownames
Affy_U133B_data <- cbind(Affy_U133B_probes[intersect(rownames(Affy_U133B_data), rownames(Affy_U133B_probes)),1], Affy_U133B_data)

Affy_U133B_data  <- multiGene2ProbeFilter(Affy_U133B_data)

dim(Affy_U133B_data)

cat("Writing normalised expression data to Affy_U133B_gcrma.exp\n")
write.table(prepare2write(Affy_U133B_data), file="Affy_U133B_gcrma.exp",sep="\t", row.names=FALSE)


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
