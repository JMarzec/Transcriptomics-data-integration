################################################################################
#
#   File name: GenExpressionComb_Illum_HT_12_V3.R
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
#	Description: Pipeline combining datasets produced using Illumina Human HT-12 v3 expression beadchip. The script reads in three a priori generated files: (1) list of datasets to be analysed; (2) file produced with QC_Illum_HT_12_V3.R which lists outliers identified with ArrayOutliers function and (3) '.flat' file produced with ArrayAnnot.pl which lists reliable probes annotated with Ensembl gene IDs and gene symbol. Pipeline executes 'multiGene2ProbeFilter.R' to remove multi gene to probe mappings
#
#	Command line use example: R --file=./GenExpressionComb_Illum_HT_12_V3.R --args "/scratch/jack/data/PhD/Transcriptomics_project"
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

#####  Bind datasets retaining only the common list of genes between all studies
bindDatasets <- function(datasetList) {
    
    cat(paste("Loading expression data", datasetList[[1]], "...\n"), sep="")
    dataBound <- read.table(datasetList[[1]],sep="\t",as.is=TRUE,header=FALSE,row.names=1)
    
	for (i in 2:length(datasetList)) {
		
        cat(paste("Loading expression data", datasetList[[i]], "...\n"), sep="")
        data2bind <- read.table(datasetList[[i]],sep="\t",as.is=TRUE,header=FALSE,row.names=1)
        
        dataBound <- cbind(dataBound[ intersect(rownames(dataBound), rownames(data2bind)),], data2bind[ intersect(rownames(dataBound), rownames(data2bind)),])
	}
    return( dataBound )
}

##### Prepare object to write into a file
prepare2write <- function (x) {
	
	x2write <- cbind(rownames(x), x)
    colnames(x2write) <- c("",colnames(x))
	return(x2write)
}

##### Replace signal intensities of '0' with minimal value detected on the array
replace0s <- function (lumiBatch) {
	
	for (i in 1:ncol(exprs(lumiBatch))) {
		
		exprs.ordered <- sort(exprs(lumiBatch)[,i])
		exprs.ordered <- exprs.ordered[ exprs.ordered != 0 ]
        
		exprs(lumiBatch)[exprs(lumiBatch)[,i]==0] <- exprs.ordered[1]
	}
	return(lumiBatch)
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

library(lumi)
library(simpleaffy)
library(sva)

source("MultiGene2ProbeFilter.R")

#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

ProjectDir = args[4]
probesFlatFile = paste(paste(ProjectDir, collapse = "/"),"Illum_HT_12_V3.flat",sep="/")

##### Read file with datasets information
DatasetInput=read.table(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE,row.names=1)
DatasetInput=DatasetInput[DatasetInput[,"Platform"] %in% "Illum_HT_12_V3",]

##### Change working directory to the project workspace
setwd(ProjectDir)

#===============================================================================
#    Load data
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

datasetList <- as.list(unique(targetFile[,"FileName"]))

if ( length(datasetList) > 1 ) {

    datasets <- bindDatasets(datasetList)
    cat("Combining datasets produced using Illumina HumanHT-12 v3.0 Expression Beadchip ...\n")
    cat("Writing bound expression data into a file.\n")
    write.table(prepare2write(datasets), file="Illum_HT_12_V3_exprs_matrix.txt",sep="\t", row.names=FALSE, col.names = FALSE)
    
    cat("Loading bound expression data into LumiBatch object\n")
    dat <- lumiR("Illum_HT_12_V3_exprs_matrix.txt", sep = NULL, detectionTh = 0.01, na.rm = TRUE, lib = NULL)
    
} else {
    cat("Loading dataset produced using Illumina HumanHT-12 V3.0 Expression Beadchip into LumiBatch object...\n")
    
    dat <- lumiR(datasetList[[1]], sep = NULL, detectionTh = 0.01, na.rm = TRUE, lib = NULL)
}

##### Remove samples identified as outliers
studyIDs <- rownames(DatasetInput)

sampleNames(dat) <- rownames(targetFile)
dat@phenoData$sampleID <- rownames(targetFile)


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
dat <- get.array.subset(dat,"sampleID",setdiff(sampleNames(dat),as.vector(as.list(outlierList))))

targetFile <- targetFile[setdiff(rownames(targetFile),as.vector(as.list(outlierList))),]


##### Keep only reliable probes based on cigar strings derived from Ensembl microarray probe mapping
Illum_HT_12_V3_probes <- read.table(probesFlatFile,sep="\t",as.is=TRUE,header=TRUE,row.names=1)

dat <- dat[featureNames(dat) %in% rownames(Illum_HT_12_V3_probes), ]

#===============================================================================
#    Normalisation with Robust Spline Normalisation (RSN)
#===============================================================================
##### RSN method combines features from loess and quantile normalization methods. It was recommended by Schmid et al ( http://www.biomedcentral.com/1471-2164/11/349 ) as it  combines the positive effects of quantile normalization, i.e. preservation of the rank order, and spline interpolation, i.e. continuous mapping of the values, and at the same time circumvents their drawbacks, i.e. discontinuous mapping of intensity values and no rank preservation, respectively

##### Replace signal intensities of '0', which introduce errors in RSN, with minimal value detected on the array
dat  <- replace0s(dat)
Illum_HT_12_V3_rsn <- lumiN(dat, method='rsn')

Illum_HT_12_V3_data <- exprs(Illum_HT_12_V3_rsn)

dim(Illum_HT_12_V3_data)

##### Add '1' to all values to avaid issues with expression of 0 in log2-transformation. Expression of '1' will become '0' after log2-transformation
Illum_HT_12_V3_data <- Illum_HT_12_V3_data+1
Illum_HT_12_V3_data <- log2(Illum_HT_12_V3_data)

#===============================================================================
#    Boxplot on normalised data
#===============================================================================

Illum_HT_12_V3_datasets <- targetFile$Dataset

Illum_HT_12_V3_datasets_No <- max(as.numeric(factor(Illum_HT_12_V3_datasets)))

Illum_HT_12_V3_datasets.colour <- getDatasetsColours(Illum_HT_12_V3_datasets)


pdf("Illum_HT_12_V3_boxplot.pdf", pointsize = 8 ,width = 0.2*length(sampleNames(dat)), height = 6)
par(mar=c(13, 4, 3, 2))
boxplot(dat,col=Illum_HT_12_V3_datasets.colour[[2]], las = 2, ylab="") # Generates boxplot of non-normalized intensity values.
boxplot(data.frame(Illum_HT_12_V3_data),col=Illum_HT_12_V3_datasets.colour[[2]], main="Normalised RSN data", las = 2) # Generates boxplot of normalized log intensity values.
dev.off()

#===============================================================================
#    Remove multi gene to probe mappings
#===============================================================================

#####  Dataset is presumed to contain non-unique gene names, so write the gene names into the first column rather then into the rownames
Illum_HT_12_V3_data <- cbind(Illum_HT_12_V3_probes[intersect(rownames(Illum_HT_12_V3_data), rownames(Illum_HT_12_V3_probes)),1], Illum_HT_12_V3_data)

Illum_HT_12_V3_data  <- multiGene2ProbeFilter(Illum_HT_12_V3_data)

dim(Illum_HT_12_V3_data)
colnames(Illum_HT_12_V3_data) <- sampleNames(dat)

cat("Writing normalised expression data to Illum_HT_12_V3_rsn.exp\n")
write.table(prepare2write(Illum_HT_12_V3_data), file="Illum_HT_12_V3_rsn.exp",sep="\t", row.names=FALSE)


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
