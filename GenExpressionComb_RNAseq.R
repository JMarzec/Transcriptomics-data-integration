################################################################################
#
#   File name: GenExpressionComb_RNAseq.R
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
#	Description: Pipeline merging datasets produced using RNA-seq technology. The script reads a priori generated file containing a single gene-by-sample expression (read counts) matrix produced with GenExpressionComb_RNAseq.pl. It reads in [GTF]_gene_info.txt file, generated with Perl script 'Get_gene_info.pl', containing genes' length and GC content required for conditional quantile normalisation (CQN)
#
#	Command line use example: R --file=./GenExpressionComb_RNAseq.R --args "/scratch/jack/genome_annotation/ensembl_Homo_sapiens.GRCh37.74.gtf.gene_info.txt" "/scratch/jack/data/PhD/Transcriptomics_project"
#
#	First arg:      [GTF]_gene_info.txt file
#	Second arg:     Project workspace. This is the directory to which cross-platform analysis results will be written
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

library(preprocessCore)
library(cqn)

#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

GenesInfo = args[4]
ProjectDir = args[5]

##### Read file with datasets information
DatasetInput=read.table(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE,row.names=1)
DatasetInput=DatasetInput[DatasetInput[,"Platform"] %in% "RNAseq",]

##### Change working directory to the project workspace
setwd(ProjectDir)

#===============================================================================
#    Load data
#===============================================================================

fileInfo = strsplit(DatasetInput[,"TargetFile"], split='/', fixed=TRUE)
targetFile <- read.table(DatasetInput[1,"TargetFile"],sep="\t",as.is=TRUE,header=TRUE)[,c(1:4)]
targetFile <- targetFile[!duplicated(targetFile[,"Name"]),]
rownames(targetFile) <- targetFile[,"Name"]
targetFile <- cbind(targetFile[,2:4],rownames(DatasetInput[1,]))
colnames(targetFile)[ncol(targetFile)] <- "Dataset"

if ( nrow(DatasetInput) > 1 ) {
    for ( i in 2:nrow(DatasetInput) ) {
    
        targetFileTmp <- read.table(DatasetInput[i,"TargetFile"],sep="\t",as.is=TRUE,header=TRUE)[,c(1:4)]
        targetFileTmp <- targetFileTmp[!duplicated(targetFileTmp[,"Name"]),]
        rownames(targetFileTmp) <- targetFileTmp[,"Name"]
        targetFileTmp <- cbind(targetFileTmp[,2:4],rownames(DatasetInput[i,]))
        colnames(targetFileTmp)[ncol(targetFileTmp)] <- "Dataset"
    
        ##### Deal with replicates
        if ( any(!is.na(targetFile[,"Replicate"])) ) {
            maxRep <- max(targetFile[!is.na(targetFile[,"Replicate"]),"Replicate"])
            targetFileTmp[,"Replicate"] <- targetFileTmp[,"Replicate"] + maxRep
        }
        targetFile <- rbind(targetFile, targetFileTmp)
    }
}


cat("Loading RNA-seq data...\n")
dat <- read.table(paste(ProjectDir,"RNAseq_count_matrix.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE,row.names=1,check.names=FALSE)

GenesInfo <- read.table(GenesInfo,header=TRUE, row.name=1, sep="\t")

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
dat <- dat[,setdiff(colnames(dat),as.vector(as.list(outlierList)))]

targetFile <- targetFile[setdiff(rownames(targetFile),as.vector(as.list(outlierList))),]

RNAseq_datasets <- targetFile$Dataset

RNAseq_datasets_No <- max(as.numeric(factor(RNAseq_datasets)))

RNAseq_datasets.colour <- getDatasetsColours(RNAseq_datasets)

pdf("RNAseq_libSize.pdf", pointsize = 10 ,width = 0.2*ncol(dat), height = 6)
par(las=3, mar=c(10,5,1,0))
barplot(colSums(dat)*1e-6, names=colnames(dat), col=RNAseq_datasets.colour[[2]], ylab="Library size (millions)")
dev.off()

dim(dat)

##### Store genes with raw counts values of '0' across all samples in separate matrix
dat0s <- dat[apply(dat, 1, function(x) sum(x)==0),]
dat <- dat[apply(dat, 1, function(x) sum(x)!=0),]
GenesInfo <- GenesInfo[ rownames(GenesInfo) %in% rownames(dat),]

##===============================================================================
##    Normalisation using CQN (conditional quantile normalisation) method
##===============================================================================

targets <- targetFile$Target


targetLevels <- unique(targets)

for ( i in 1:length(targetLevels) ) {
    
    targets <- gsub(targetLevels[i], i-1, targets)
    
}

CompComb <- combn(c(1:length(targetLevels)), 2)


###### Extract infromation about genes present in expression matrix
GenesInfo <- GenesInfo[rownames(GenesInfo) %in% rownames(dat),]
dat <- dat[rownames(dat) %in% rownames(GenesInfo),]

###### Order gene information accordingly to the order of genes present in expression matrix
GenesInfo <- GenesInfo[order(match(rownames(GenesInfo),rownames(dat))),]

###### Check if the row ordering of the count matrix is the same as the row ordering of the matrix containing length and GC-content
stopifnot(all(rownames(dat) == rownames(GenesInfo)))


RPKM.cqn <- NULL

pdf("RNAseq_cqnplot.pdf", pointsize = 8 ,width = 8, height = 4)
for (i in 1:RNAseq_datasets_No) {
	
	cat(paste("Conditional quantile normalisation of", studyIDs[i], "dataset...\n", sep=" "))
    dat_study <- dat[, RNAseq_datasets==studyIDs[i]]
	cqn.subset <- cqn(dat_study, lengths = GenesInfo$length, x = GenesInfo$GC_content, sizeFactors = NULL, verbose = TRUE)
    
    ##### Change normalised log2 values < 2 to 0 (raw counts < 4, assume that less than 4 reads are insufficient to accuretly represent gene expression)
    RPKM.cqn.log2 <- cqn.subset$y + cqn.subset$offset
    RPKM.cqn.log2[RPKM.cqn.log2 < 2] <- 0
    
	RPKM.cqn <- cbind(RPKM.cqn, RPKM.cqn.log2)

	###### Examine plots of systematic effects (gene length and GC content)
	par(mfrow=c(1,2))
	cqnplot(cqn.subset, n = 1, xlab = "GC content", lty = 1, ylim = c(1,7), main=studyIDs[i])
	cqnplot(cqn.subset, n = 2, xlab = "length", lty = 1, ylim = c(1,7), main=studyIDs[i])
    
    CompComb_study <- combn(c(1:length(unique(targets[RNAseq_datasets==studyIDs[i]]))), 2)
    
    ###### Generate MA plots
    for (j in 1:ncol(CompComb_study)) {
        
        target_study <- targets[RNAseq_datasets==studyIDs[i]]
        
        ###### Make sure that each group is represented by at least two samples
        if ( min(table(target_study)) > 1 ) {
        
            if ( table(target_study==CompComb[1,j]-1)[1] > 1 && table(target_study==CompComb[1,j]-1)[2] > 1 ) {
                CompName <- paste(sort(unique(targetFile$Target)[CompComb[,j]]),collapse=" vs ")

                ##### Plot only genes with average log2 >=2 of the CQN-normalised data
                whGenes <- which(rowMeans(RPKM.cqn.log2) >= 2 )
        
                M.std <- rowMeans( log2(dat_study[whGenes, target_study==CompComb[1,j]-1]) ) - rowMeans( log2(dat_study[whGenes, target_study==CompComb[2,j]-1]) )
                A.std <- rowMeans( log2(dat_study[whGenes,]) )
                M.cqn <- rowMeans( RPKM.cqn.log2[whGenes, target_study==CompComb[1,j]-1] ) - rowMeans( RPKM.cqn.log2[whGenes, target_study==CompComb[2,j]-1] )
                A.cqn <- rowMeans( RPKM.cqn.log2[whGenes,] )
        
                gccontent <- GenesInfo$GC_content[whGenes]
                whHigh <- which(gccontent > quantile(gccontent, 0.9))
                whLow <- which(gccontent < quantile(gccontent, 0.1))
        
                ###### MA plot Examine plots for non-normalised data with genes coloured according to whether they have high (red) or low (blue) GC-content
                plot(A.std, M.std, cex=0.3, pch=16, xlab="A", ylab="M", main=paste(CompName, "(Non-normalised)", sep=" "), ylim=c(-4,4), xlim=c(0,15), col=rgb(0, 0, 0, alpha=0.5))
                points(A.std[whHigh], M.std[whHigh], cex=0.5, pch=16, col = "red", xlab="A", ylab="M", main=paste(CompName, "(Non-normalised)", sep=" "), ylim=c(-4,4), xlim=c(0,15))
                points(A.std[whLow], M.std[whLow], cex=0.5, pch=16, col = "blue")
                abline(h=0)
        
                ###### MA plot Examine plots for CQN normalised data with genes coloured according to whether they have high (red) or low (blue) GC-content
                plot(A.cqn, M.cqn, cex=0.3, pch=16, xlab="A", ylab="M", main=paste(CompName, "(CQN normalised)", sep=" "), ylim=c(-4,4), xlim=c(0,15), col=rgb(0, 0, 0, alpha=0.5))
                points(A.cqn[whHigh], M.cqn[whHigh], cex=0.5, pch=16, col = "red", xlab="A", ylab="M", main=paste(CompName, "(CQN normalised)", sep=" "), ylim=c(-4,4), xlim=c(0,15))
                points(A.cqn[whLow], M.cqn[whLow], cex=0.5, pch=16, col = "blue")
                abline(h=0)
            }
        }
    }
}
dev.off()


##### Store genes with log2 values < 2 across all samples in separate matrix
RPKM.cqn0s <- RPKM.cqn[apply(RPKM.cqn, 1, function(x) sum(x)==0),]
dat0s <- rbind(dat0s, RPKM.cqn0s)
RPKM.cqn <- RPKM.cqn[apply(RPKM.cqn, 1, function(x) sum(x)!=0),]


pdf("RNAseq_boxplot_cqn.pdf", pointsize = 8 ,width = 0.2*ncol(dat), height = 6)
par(mar=c(13, 4, 3, 2))
#boxplot(log2(dat),col=RNAseq_datasets.colour[[2]], las = 2, ylab="") # Generates boxplot of non-normalized FPKM values.
boxplot(data.frame(RPKM.cqn),col=RNAseq_datasets.colour[[2]], main="Quantile normalised data", las = 2) # Generates boxplot of 
dev.off()


##### Convert log2 values into normalised read counts
#RPKM.cqn <- 2^RPKM.cqn
#RPKM.cqn[RPKM.cqn == 1] <- 0

cat("Writing CQN expression data to RNAseq_cqn.exp\n")
write.table(prepare2write(rbind(RPKM.cqn, dat0s)), file="RNAseq_cqn.exp",sep="\t", row.names=FALSE)


##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
