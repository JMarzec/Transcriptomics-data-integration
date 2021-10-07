################################################################################
#
#   File name: GenExpressionComb_Affy_U133Plus2.R
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
#	Description: Pipeline combining datasets produced using Affymetrix Human Genome U133 Plus 2 Array. The script reads in three a priori generated files: (1) list of datasets to be analysed; (2) file produced with QC_Affy_U133Plus2.R which lists outliers identified with ArrayOutliers function and (3) '.flat' file produced with ArrayAnnot.pl which lists reliable probes annotated with Ensembl gene IDs and gene symbol. Pipeline executes 'multiGene2ProbeFilter.R' to remove multi gene to probe mappings
#
#	Command line use example: R --file=./GenExpressionComb_Affy_U133Plus2.R --args "/scratch/jack/data/PhD/Transcriptomics_project"
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
library(biomaRt)

source("MultiGene2ProbeFilter.R")

#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

ProjectDir = args[4]

##### Read file with datasets information
DatasetInput=read.table(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE,row.names=1)
DatasetInput=DatasetInput[DatasetInput[,"Platform"] %in% "Affy_U133Plus2",]

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

cat("Loading datasets produced using Affymetrix GeneChip Human Genome U133 Plus 2.0 Array into an AffyBatch object...\n")

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

##### Select an affybatch probeset subset prior to normalization, based on the following post https://stat.ethz.ch/pipermail/bioconductor/2009-January/026046.html
library(hgu133plus2cdf)

dat<-ReadAffy(filenames = targetFile[,"FileName"], sampleNames = rownames(targetFile), phenoData=as(as.data.frame(targetFile[,"Target"]), "AnnotatedDataFrame") )
colnames(pData(dat)) <- "Target"

#===============================================================================
#    Normalisation with GC-RMA
#===============================================================================

cat("Robust Multiarray Average normalisation with GC-content background correction.\n")
Affy_U133Plus2_gcrma <- NULL
Affy_U133Plus2_gcrma <- gcrma(dat)

#save.image("Affy_U133Plus2_gcrma.R")
#load("Affy_U133Plus2_gcrma.R")

Affy_U133Plus2_data <- exprs(Affy_U133Plus2_gcrma)

dim(Affy_U133Plus2_data)

#===============================================================================
#    Boxplot on normalised data
#===============================================================================

Affy_U133Plus2_datasets <- targetFile$Dataset

Affy_U133Plus2_datasets_No <- max(as.numeric(factor(Affy_U133Plus2_datasets)))

Affy_U133Plus2_datasets.colour <- getDatasetsColours(Affy_U133Plus2_datasets)


pdf("Affy_U133Plus2_boxplot.pdf", pointsize = 8 ,width = 0.2*length(targetFile$FileName), height = 6)
par(mar=c(13, 4, 3, 2))
boxplot(dat,col=Affy_U133Plus2_datasets.colour[[2]], las = 2) # Generates boxplot of non-normalized intensity values.
boxplot(data.frame(Affy_U133Plus2_data),col=Affy_U133Plus2_datasets.colour[[2]], main="Normalised GC-RMA data", las = 2) # Generates boxplot of normalized log intensity values.
dev.off()















#===============================================================================
#    Convert probes to gene symbols
#===============================================================================

##### Retrieve gene annotation information
annotGenes <- function (topGenes) {


    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    #genesAnnot <- getGene( id=rownames(Affy_U133Plus2_data), type="ensembl_gene_id", mart = mart)  ##### Currently does not work
    #genesAnnot <- getBM(c("ensembl_gene_id","hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position"), "ensembl_gene_id", rownames(Affy_U133Plus2_data), mart = mart)

    ensembl = useMart(biomart= "ensembl",dataset="hsapiens_gene_ensembl")
    affy_ensembl= c("affy_hg_u133_plus_2", "ensembl_gene_id")
    genesAnnot <- getBM(attributes= affy_ensembl, mart= ensembl, values = "*", uniqueRows=T)

    genesAnnotUnique <- unique(genesAnnot$ensembl_gene_id)
    annotated <- NULL

    for (j in 1:length(genesAnnotUnique)) {

        annotated <- rbind(annotated, c(genesAnnot[genesAnnot$ensembl_gene_id %in% genesAnnotUnique[j], 1:8], Affy_U133Plus2_data[genesAnnotUnique[j],]) )

    }
    rownames(annotated) <- annotated[,"ensembl_gene_id"]

    ##### Add not annotated genes
    notAnnotated <- topGenes[rownames(topGenes) %!in% genesAnnot$ensembl_gene_id,]
    notAnnotated <- cbind(as.data.frame(matrix(NA,nrow(notAnnotated),8)), notAnnotated)
    colnames(notAnnotated) <- colnames(annotated)

    annotated <- as.matrix(rbind(annotated, notAnnotated))

	return(annotated[,2:ncol(annotated)])
}




ensembl = useMart(biomart= "ensembl",dataset="hsapiens_gene_ensembl")
affy_ensembl= c("affy_hg_u133_plus_2", "ensembl_gene_id")
getBM(attributes= affy_ensembl, mart= ensembl, values = "*", uniqueRows=T)


#===============================================================================
#    Remove multi gene to probeset mappings
#===============================================================================

#####  Dataset is presumed to contain non-unique gene names, so write the gene names into the first column rather then into the rownames
Affy_U133Plus2_data <- cbind(Affy_U133Plus2_probes[intersect(rownames(Affy_U133Plus2_data), rownames(Affy_U133Plus2_probes)),1], Affy_U133Plus2_data)

Affy_U133Plus2_data  <- multiGene2ProbeFilter(Affy_U133Plus2_data)

dim(Affy_U133Plus2_data)




#===============================================================================
#    Write data into a file
#===============================================================================

cat("Writing normalised expression data to Affy_U133Plus2_gcrma.exp\n")
write.table(prepare2write(Affy_U133Plus2_data), file="Affy_U133Plus2_gcrma.exp",sep="\t", row.names=FALSE)


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
