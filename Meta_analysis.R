################################################################################
#
#   File name: Meta_analysis.R
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
#   Description: Pipeline for meta-analysis of defined datasets. It performs integrative correlation analysis using functions from MergeMaid package and calcualtes between and within study variability using functions from GeneMeta package. The script generates relevant plots to help estimate between and within study variability and evaluate the impact of meta-analysis on the final results. The p-values computed for individual platforms are combined using Stouffer's method with fold-change (log2) and integrative correlation (reproducibility score for each gene) used as weighting factors. The script reads in a priori generated files with normalised and batch-effect adjusted expression data as well and differential expression analysis results obtained for each platform. It invokes Ghostscript to merge .pdf files.
#
#   Command line use example: R --file=./Meta_analysis.R --args "/scratch/jack/data/PhD/Transcriptomics_project" "/scratch/jack/data/PhD/Transcriptomics_project/PCa_assocaited_genes_DDPC-Schoenborn-COSMIC.txt"
#
#   First arg:		Project workspace. This is the directory to which cross-platform analysis results will be written
#   Second arg (OPTIONAL): List of known genes associated with studied phenotype. The first column is expected to list the Ensembl gene IDs
#
################################################################################

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

#===============================================================================
#    Functions
#===============================================================================

##### Merge .pdf files
mergePDF <- function(infiles, outfile, os = "UNIX") {
	
	version <- switch(os, UNIX = "gs", Win32 = "gswin32c", Win64 = "gswin64c")
	pre <- " -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="
	system(paste(paste(version, pre, outfile, sep = ""), infiles, collapse = " "))
}

##### Prepare object to write into a file
prepare2write <- function (x) {
	
	x2write <- cbind(rownames(x), x)
    colnames(x2write) <- c("",colnames(x))
	return(x2write)
}

##### Panel function calculating correlation coeffcient to be used for pairwise plots
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor/3, col="blue")
}


##### Fit a loess smooth curve to the data
loessSmoothing <- function (x,y) {
    
    lo <-list()
    lo[[1]] <- loess(y~x)
    lo[[2]] <- seq(min(x),max(x), (max(x) - min(x))/1000)
    
    return(lo)
}

##### Create 'not in' operator
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

##### Convert a p-value into a zscore ( formula taken from http://www.broadinstitute.org/~debakker/meta.R )
convert.pvalue <- function(pval, FC) {
    
    z <- NULL
    
    for (i in 1:length(FC)) {
        if ( FC[i] > 0 ) {
            z <- c(z, qnorm( pval[i] / 2 ))
        } else {
            z <- c(z, -(qnorm( pval[i] / 2 )))
        }
    }
    return(z)
}

##### Compute combined zscore ( formula taken from http://www.broadinstitute.org/~debakker/meta.R )
p.to.z <- function(pval, FC, weigth) {
    
    z <- NULL
    zSum <- 0
    
    for (i in 1:length(pval)) {
        
        z[i] <- convert.pvalue(pval[i], FC[i])
        zSum <- zSum + weigth[i] *z[i]
    }
    return(zSum)
}

##### Compute the weightings for meta-analysis with Stouffer's method
weighting <- function(FC, intCor) {
    
    w <- NULL
    wSum <- 0
    
    for (i in 1:length(FC)) {
        
        w[i] <- abs(FC[i]) + (abs(FC[i]) * intCor^2 )
        w[i][w[i] < 0] <- 0
        wSum <- wSum + w[i]
    }
    
    for (i in 1:length(FC)) {
        w[i] <- sqrt(abs(FC[i])) / sqrt(wSum)
    }
    return(w)
}

##### Remove rows cointaining NAs in the list
list.rm.na <- function(list) {
    
    for (i in 1:length(list)) {
    
        list[[i]] <- list[[i]][complete.cases(list[[i]]),]
    }
    return(list)
}

##### Retrieve gene annotation information
annotGenes <- function (topGenes) {
	
    
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    genesAnnot <- getGene( id=rownames(topGenes), type="ensembl_gene_id", mart = mart)
    
    genesAnnotUnique <- unique(genesAnnot$ensembl_gene_id)
    annotated <- NULL
    
    for (j in 1:length(genesAnnotUnique)) {
        
        annotated <- rbind(annotated, c(genesAnnot[genesAnnot$ensembl_gene_id %in% genesAnnotUnique[j], 1:8], topGenes[genesAnnotUnique[j],]) )
        
    }
    rownames(annotated) <- annotated[,"ensembl_gene_id"]
    
    ##### Add not annotated genes
    notAnnotated <- topGenes[rownames(topGenes) %!in% genesAnnot$ensembl_gene_id,]
    notAnnotated <- cbind(as.data.frame(matrix(NA,nrow(notAnnotated),8)), notAnnotated)
    colnames(notAnnotated) <- colnames(annotated)
    
    annotated <- as.matrix(rbind(annotated, notAnnotated))
    
	return(annotated[,2:ncol(annotated)])
}

#===============================================================================
#    Load libraries
#===============================================================================

library(gplots)
library(MergeMaid)
library(GeneMeta)
library(biomaRt)
library(org.Hs.eg.db)
library(topGO)
library(MASS)


#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

ProjectDir = args[4]
knownGenes = args[5]

##### Read file with datasets information
DatasetInput=read.table(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE,row.names=1)


##### Read file with known genes assocaited with studied phenotype
if  ( !is.na(knownGenes) ) {
    knownGenes=read.table(knownGenes,sep="\t",as.is=TRUE,header=TRUE)
}

##### Change working directory to the project workspace
setwd(ProjectDir)


#===============================================================================
#     Integrative correlation analysis using functions from MergeMaid package
#===============================================================================

fileInfo = strsplit(DatasetInput[,"TargetFile"], split='/', fixed=TRUE)
targetFile <- read.table(DatasetInput[1,"TargetFile"],sep="\t",as.is=TRUE,header=TRUE)[,c(1:4)]
targetFile <- targetFile[!duplicated(targetFile[,"Name"]),]
#rownames(targetFile) <- targetFile[,"Name"]
targetFile <- targetFile[,!colnames(targetFile) %in% "FileName"]
targetFile <- cbind(targetFile,rownames(DatasetInput[1,]), DatasetInput[1,"Platform"])
colnames(targetFile)[c(ncol(targetFile)-1,ncol(targetFile))] <- c("Dataset","Platform")

if ( nrow(DatasetInput) > 1 ) {
    for ( i in 2:nrow(DatasetInput) ) {
    
        targetFileTmp <- read.table(DatasetInput[i,"TargetFile"],sep="\t",as.is=TRUE,header=TRUE)[,c(1:4)]
        targetFileTmp <- targetFileTmp[!duplicated(targetFileTmp[,"Name"]),]
        #rownames(targetFileTmp) <- targetFileTmp[,"Name"]
        targetFileTmp <- targetFileTmp[,!colnames(targetFileTmp) %in% "FileName"]
        targetFileTmp <- cbind(targetFileTmp,rownames(DatasetInput[i,]), DatasetInput[i,"Platform"])
        colnames(targetFileTmp)[c(ncol(targetFileTmp)-1,ncol(targetFileTmp))] <- c("Dataset","Platform")
    
        ##### Deal with replicates
        if ( any(!is.na(targetFile[,"Replicate"])) ) {
            maxRep <- max(targetFile[!is.na(targetFile[,"Replicate"]),"Replicate"])
            targetFileTmp[,"Replicate"] <- targetFileTmp[,"Replicate"] + maxRep
        }
        targetFile <- rbind(targetFile, targetFileTmp)
    }
}

dataset.expFiles <- list()
dataset.exp <- list()

platforms_No <- length(unique(DatasetInput[,"Platform"]))

for ( i in 1:platforms_No ) {

    dataset.expFiles[[i]] <- paste(ProjectDir, "/Comb_", sep="")
}


for ( i in 1:platforms_No ) {
	
	if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_HuEx1ST" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "Affy_HuEx1ST",])),collapse = "_" ), sep="" )
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_U133Plus2" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "Affy_U133Plus2",])),collapse = "_" ), sep="" )
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_U133A2" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "Affy_U133A2",])),collapse = "_" ), sep="" )
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_U133A" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "Affy_U133A",])),collapse = "_" ), sep="" )
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_U133B" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "Affy_U133B",])),collapse = "_" ), sep="" )
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_U95Av2" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "Affy_U95Av2",])),collapse = "_" ), sep="" )
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_U95B" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "Affy_U95B",])),collapse = "_" ), sep="" )
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_U95C" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "Affy_U95C",])),collapse = "_" ), sep="" )
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Illum_HT_12_V3" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "Illum_HT_12_V3",])),collapse = "_" ), sep="" )
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "RNAseq" ) {
        
        dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], paste(sort(rownames(DatasetInput[DatasetInput$Platform %in% "RNAseq",])),collapse = "_" ), sep="" )
        
	} else {
		cat(paste("Platform", unique(factor(DatasetInput[,"Platform"]))[i] ,"is not yet supported\n\n", sep=" "))
		q()
	}
}

for ( i in 1:platforms_No ) {
    
    dataset.expFiles[[i]] <- paste(dataset.expFiles[[i]], "exp", sep=".")
}


#===============================================================================
#    Load normalised data from all platforms
#===============================================================================

merged <- list()
geneNames <- NULL

##### Merge gene expression datasets
for ( i in 1:platforms_No ) {
        
    cat(paste("\nUploading expression matrix for platform ", unique(factor(DatasetInput[,"Platform"]))[i] ,"...\n", sep=""))
    exp <- as.matrix(read.table(dataset.expFiles[[i]],sep="\t",as.is=TRUE,header=TRUE,row.names=1,check.names=FALSE))
    merged[[i]] <-mergeExprs(exp)
    
    geneNames <- c(geneNames, rownames(exp))
}
geneNames <- unique(geneNames)


cat("\nPerforming integrative correlation analysis...\n\n")

integrativeCors <- sapply(geneNames,function(x) NULL)

##### Compute correlation coefficients for all possible combinations of platforms
PlatformComb <- combn(c(1:platforms_No), 2)

for (i in 1:ncol(PlatformComb)) {
    
    mergeCommand <- NULL

    for ( j in 1:2 ) {

        mergeCommand <- c(mergeCommand,paste("merged[[", PlatformComb[j,i], "]]", sep=""))
    }

    mergedSub <- eval(parse(text = paste("mergeExprs(", paste(mergeCommand, collapse=","), ")")))
    
    ##### Get the genes intersection to work with datasets with the same gene sets
    geneSub <- nrow(exprs(MergeMaid::intersection(mergedSub)))

    ##### Compute integrative correlation coefficients !!!very time consuming step depending on number of genes and to some extend on sample size of each study!!!
    cat(paste("\nCalculating integrative correlations for ", geneSub, " genes common between platforms ", paste(unique(factor(DatasetInput[,3]))[PlatformComb[,i]], collapse=" and "), "...\n", sep=""))
    corcor <- intCor(mergedSub,method= "pearson",exact=TRUE)
    
    for (j in 1:length(geneNames)) {
        integrativeCors[[geneNames[j]]] <- c(integrativeCors[[geneNames[j]]], integrative.cors(corcor)[geneNames[j]])
    }
    
    ##### Save R session for meta-analysis
    save.image("Meta.R")
}

##### Get the average integrative correlation coefficients for each gene
for (i in 1:length(geneNames)) {
    integrativeCors[[geneNames[i]]] <- mean(integrativeCors[[geneNames[i]]], na.rm=TRUE)
}
integrativeCors <- unlist(integrativeCors[!is.na(integrativeCors)])


##### Compute correlation coefficients for all platforms
mergeCommand <- NULL

for ( i in 1:platforms_No ) {
    
    mergeCommand <- c(mergeCommand,paste("merged[[", i, "]]", sep=""))
}

merged <- eval(parse(text = paste("mergeExprs(", paste(mergeCommand, collapse=","), ")")))

##### Obtain the basic information about each dataset
#summary(merged)
##### Genes present in at least two of the studies
#geneStudy(merged)

##### Get platforms sizes with names and sample names
platformDim <- NULL
platformNames <- NULL
sampleNames <- NULL

for (i in 1:length(merged)) {
    
    platformDim[i] <- summary(merged)$`Number of Samples in Each Study`[2,i]
    platformNames[i] <- as.vector(unique(factor(DatasetInput[,"Platform"]))[i])
    sampleNames <- c(sampleNames, rownames(data.frame(phenoData(merged)[i])))
}

names(merged@data) <- platformNames
#summary(merged)

##### Compute integrative correlation coefficients !!!very time consuming step depending on number of genes and to some extend on sample size of each study!!!
cat("\nCalculating integrative correlations for all platforms\n\n")
corcor <- intCor(merged,method="pearson",exact=TRUE)

corcor@notes <- platformNames


##### Save R session for meta-analysis
save.image("Meta.R")

##### Draw histograms of integrative correlation
pdf("Meta_IntegrativeCor_hist.pdf", width=2*(platforms_No+1), height=2*(platforms_No+1))
hist(corcor)
dev.off()

##### Plot the distribution of the integrative correlation coefficients and the null distribution obtained by permutation
#pdf("Meta_IntegrativeCor_dens.pdf", width=2*(platforms_No+1), height=2*(platforms_No+1))
#corDens <- intcorDens(merged, method="pearson") # Error in density.default(d) : 'x' contains missing values
#dev.off()

##### Keep only samples present in normalised data files
targetFile <- targetFile[ paste(targetFile[,"Name"], targetFile[,"Platform"], sep="_") %in% paste(sampleNames, rep(platformNames,as.numeric(platformDim)), sep="_"),]

targetFile <- targetFile[order(charmatch(paste(targetFile[,"Name"], targetFile[,"Platform"], sep="_"), paste(sampleNames, rep(platformNames,as.numeric(platformDim)), sep="_"))),]

targets <- targetFile$Target

targetLevels <- unique(targets)

for ( i in 1:length(targetLevels) ) {
    
	targets[targets == targetLevels[i]] <- i-1
}

##### Annotate phenoData with targets
platformStart = 0
platformEnd = 0

##### Go through remaining platforms
for (i in 1:length(merged)) {
    
    platformStart = platformEnd + 1
    platformEnd = platformStart + as.numeric(platformDim[i]) - 1
    
    phenoData(merged)[[i]][1] <- as.factor(targets[platformStart:platformEnd])
}


#===============================================================================
#     Compute study specific regression coefficients using functions from MergeMaid package
#===============================================================================
##### Consider high-coverage platforms only including RNAseq, Affy_HuEx1ST, Affy_U133Plus2 and Illum_HT_12_V3

platforms <- NULL
dataset.expFiles.part <- list()

for ( i in 1:platforms_No ) {
	
	if ( unique(factor(DatasetInput[,"Platform"]))[i] == "RNAseq" ) {
        
        dataset.expFiles.part[[i]] <- dataset.expFiles[[i]]
        platforms <- c(platforms, "RNAseq")
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_HuEx1ST" ) {
        
        dataset.expFiles.part[[i]] <- dataset.expFiles[[i]]
        platforms <- c(platforms, "Affy_HuEx1ST")
                
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Affy_U133Plus2" ) {
        
        dataset.expFiles.part[[i]] <- dataset.expFiles[[i]]
        platforms <- c(platforms, "Affy_U133Plus2")
        
	} else if ( unique(factor(DatasetInput[,"Platform"]))[i] == "Illum_HT_12_V3" ) {
        
        dataset.expFiles.part[[i]] <- dataset.expFiles[[i]]
        platforms <- c(platforms, "Illum_HT_12_V3")
        
	} else {
		dataset.expFiles.part[[i]] <- character(0)
	}
}

dataset.expFiles.part <- dataset.expFiles.part[lapply(dataset.expFiles.part,length)>0]

platforms_No.part  <- length(dataset.expFiles.part)


#===============================================================================
#    Load normalised data from high-coverage platforms
#===============================================================================


merged.part <- list()
geneNames.part <- NULL

##### Merge gene expression datasets
for ( i in 1:platforms_No.part ) {
        
    cat(paste("\nUploading expression matrix for platform ", platforms[i] ,"...\n", sep=""))
    exp <- as.matrix(read.table(dataset.expFiles.part[[i]],sep="\t",as.is=TRUE,header=TRUE,row.names=1,check.names=FALSE))
    merged.part[[i]] <-mergeExprs(exp)
    
    geneNames.part <- c(geneNames.part, rownames(exp))
}
geneNames.part <- unique(geneNames.part)


##### Merge gene expression datasets
mergeCommand <- NULL

for ( i in 1:platforms_No.part ) {
    
    mergeCommand <- c(mergeCommand,paste("merged.part[[", i, "]]", sep=""))
}

merged.part <- eval(parse(text = paste("mergeExprs(", paste(mergeCommand, collapse=","), ")")))

##### Obtain the basic information about each dataset
#summary(merged.part)
##### Genes present in at least two of the studies
#geneStudy(merged.part)

##### Get platforms sizes with names and sample names
platformDim.part <- NULL
platformNames.part <- NULL
sampleNames.part <- NULL

for (i in 1:length(merged.part)) {
    
    platformDim.part[i] <- summary(merged.part)$`Number of Samples in Each Study`[2,i]
    platformNames.part[i] <- as.vector(platforms[i])
    sampleNames.part <- c(sampleNames.part, rownames(data.frame(phenoData(merged.part)[i])))
}

names(merged.part@data) <- platformNames.part
#summary(merged.part)

platforms.colour <- c("orange","darkcyan","firebrick","lightslategrey","dodgerblue","darkseagreen","wheat3","slateblue","skyblue","tomato","olivedrab")[1:(length(platformNames)+1)]
names(platforms.colour) <- c(platformNames,"Combined")


##### Perfrom this part only if there are > 1 high-coverage platforms
if ( platforms_No.part > 1 ) {
    
    ##### Compute integrative correlation coefficients !!!very time consuming step depending on number of genes and to some extend on sample size of each study!!!
    cat("\nCalculating integrative correlations for high-coverege platforms\n\n")
    corcor.part <- intCor(merged.part,method="pearson",exact=TRUE)

    corcor.part@notes <- platformNames.part

    ##### Save R session for meta-analysis
    save.image("Meta.R")

    ##### Draw histograms of integrative correlation
    pdf("Meta_IntegrativeCor_hist.tmp2.pdf", width=2*(platforms_No.part+1), height=2*(platforms_No.part+1))
    hist(corcor.part)
    dev.off()

    ##### Merge .pdf files
    file.rename("Meta_IntegrativeCor_hist.pdf", "Meta_IntegrativeCor_hist.tmp1.pdf")

    mergePDF("Meta_IntegrativeCor_hist.tmp1.pdf Meta_IntegrativeCor_hist.tmp2.pdf", "Meta_IntegrativeCor_hist.pdf")

    file.remove("Meta_IntegrativeCor_hist.tmp1.pdf")
    file.remove("Meta_IntegrativeCor_hist.tmp2.pdf")


    ##### Plot the distribution of the integrative correlation coefficients and the null distribution obtained by permutation
    pdf("Meta_IntegrativeCor_dens.pdf", width=2*(platforms_No.part+1), height=2*(platforms_No.part+1))
    corDens <- intcorDens(merged.part, method="pearson")
    dev.off()

    ##### Keep only samples present in normalised data files
    targetFile.part <- targetFile[ paste(targetFile[,"Name"], targetFile[,"Platform"], sep="_") %in% paste(sampleNames.part, rep(platformNames.part,as.numeric(platformDim.part)), sep="_"),]

    targetFile.part <- targetFile.part[order(charmatch(paste(targetFile.part[,"Name"], targetFile.part[,"Platform"], sep="_"), paste(sampleNames.part, rep(platformNames.part,as.numeric(platformDim.part)), sep="_"))),]

    targets.part <- targetFile.part$Target

    targetLevels.part <- unique(targets.part)

    for ( i in 1:length(targetLevels.part) ) {
    
        targets.part[targets.part == targetLevels.part[i]] <- i-1
    }

    ##### Annotate phenoData with targets
    platformStart = 0
    platformEnd = 0

    ##### Go through remaining platforms
    for (i in 1:length(merged.part)) {
    
        platformStart = platformEnd + 1
        platformEnd = platformStart + as.numeric(platformDim.part[i]) - 1
    
        phenoData(merged.part)[[i]][1] <- as.factor(targets.part[platformStart:platformEnd])
    }


    ##### Calculates study specific regression coefficients (using logistic regression) by correlating each gene expressions with outcome
    ##### Regression coefficients are used as validation of individual gene's ability to predict phenotypes of interest (see paper by Parmigiani G et al (2004), 'A cross-study comparison of gene expression studies for the molecular classification of lung cancer')
    log.coeff <- modelOutcome(merged.part, outcome=rep(1,length(merged.part)) , method="logistic")
    colnames(log.coeff@zscore) <- platformNames.part

    ##### Scatterplots to compare coefficients from different studies (use all genes)
    pdf("Meta_PhenoCor_logistic.pdf")
    pairs(zscore(log.coeff)[1], upper.panel = panel.smooth, lower.panel = panel.cor, pch = ".", main=paste("Z-scores (all ", nrow(zscore(log.coeff)[[1]]), " genes)", sep=""))

    ##### Scatterplots to compare coefficients from different studies (use 50% most reproducible genes)
    corGenesNo <- as.integer(length(integrative.cors(corcor.part))/2)
    TopCorGenes <- sort(integrative.cors(corcor.part))[corGenesNo:length(integrative.cors(corcor.part))]

    zscoreMatrix <- do.call(rbind, zscore(log.coeff)[1])
    zscoreMatrix <- zscoreMatrix[ rownames(zscoreMatrix) %in% names(TopCorGenes),]

    pairs(zscoreMatrix, upper.panel = panel.smooth, lower.panel = panel.cor, pch = ".", main=paste("Z-scores (", nrow(zscoreMatrix), " reproducible genes)", sep=""))
    dev.off()



    platforms.colour.part <- platforms.colour[c(platformNames.part,"Combined")]

} else {
    
    corGenesNo <- as.integer(length(integrative.cors(corcor))/2)
    
    TopCorGenes <- sort(integrative.cors(corcor))[corGenesNo:length(integrative.cors(corcor))]
}

##### Save R session for meta-analysis
save.image("Meta.R")

#===============================================================================
#     Calculate between and within study variability using functions from GeneMeta package
#===============================================================================

cat("\nCalculating between and within study variability...\n\n")

##### Remove genes with negative integrative correlation values
if ( platforms_No.part > 1 ) {
    for (i in 1:platforms_No.part) {
        merged.part@data[[i]]@assayData$exprs <- merged.part@data[[i]]@assayData$exprs[ setdiff( rownames(merged.part@data[[i]]@assayData$exprs), names(integrativeCors[integrativeCors < 0]) ), ]
        names(phenoData(merged.part)[[i]]) <- "Target"
    }
    
    merged.part@notes <- platformNames.part
    colnames(merged.part@geneStudy) <- platformNames.part
    
    ##### Get the genes intersection to work with datasets with the same gene sets
    mergedData <- MergeMaid::intersection(merged.part)
    
} else {
    
    for (i in 1:platforms_No) {
        merged@data[[i]]@assayData$exprs <- merged@data[[i]]@assayData$exprs[ setdiff( rownames(merged@data[[i]]@assayData$exprs), names(integrativeCors[integrativeCors < 0]) ), ]
        names(phenoData(merged)[[i]]) <- "Target"
    }
    
    merged@notes <- platformNames
    colnames(merged@geneStudy) <- platformNames
    
    ##### Get the genes intersection to work with datasets with the same gene sets
    mergedData <- MergeMaid::intersection(merged)
}


mergedExprs <- exprs(mergedData)

##### Get the reproducible genes common across all datasets
mergedExprsReprod <- exprs(mergedData)[rownames(mergedExprs) %in% names(TopCorGenes),]


##### Perform calculations for each comparison
CompComb <- combn(c(1:length(targetLevels)), 2)

for (i in 1:ncol(CompComb)) {
    
    esets <- list()
    esetsReprod <- list()
    
    CompNos <- CompComb[,i]
    CompNos <- CompNos-1
    targetsComp <- NULL
    
    cat(paste("\nRunning computations for", paste(unique(targetFile$Target)[CompComb[,i]],collapse=" vs "), "comparison...\n", sep=" "))
    
    CompName <- paste(sort(unique(targetFile$Target)[CompComb[,i]]),collapse="vs")
    platformsUsed <- NULL
    
    platformStart = 0
    platformEnd = 0
    
    ##### Clear the expression set of merged data
    if ( exists("mergedIntersect") ) {
        rm("mergedIntersect")
        rm("mergedReprod")
    }
    
    ##### Create an ExpressionSet object with genes common across all datasets and ExpressionSet object with reproducible genes only
    for (j in 1:length(merged.part)) {
    
        platformStart = platformEnd + 1
        platformEnd = platformStart + as.numeric(platformDim.part[j]) -1
        
        ##### Don't include dataset which doesn't contain samples from both groups to be compared
        if ( length(intersect(targets[platformStart:platformEnd], CompNos))==2  ) {
            
            platformsUsed <- c(platformsUsed, j)
            platformsUsedNo <- length(platformsUsed)
            
            if ( !exists("mergedIntersect") ) {
                mergedIntersect <- mergedExprs[,platformStart:platformEnd]
                mergedReprod <- mergedExprsReprod[,platformStart:platformEnd]
            
                ##### Select only samples included in current comparison
                mergedIntersect <- mergedIntersect[,unlist(phenoData(merged.part)[[j]][1]) %in% CompNos]
                mergedReprod <- mergedReprod[,unlist(phenoData(merged.part)[[j]][1]) %in% CompNos]
            
                esets[[platformsUsedNo]] <- new("ExpressionSet", exprs = mergedIntersect, phenoData = new("AnnotatedDataFrame", data=phenoData(merged.part)[[j]])[unlist(phenoData(merged.part)[[j]][1]) %in% CompNos,] )
                esetsReprod[[platformsUsedNo]] <- new("ExpressionSet", exprs = mergedReprod, phenoData = new("AnnotatedDataFrame", data=phenoData(merged.part)[[j]])[unlist(phenoData(merged.part)[[j]][1]) %in% CompNos,] )
                
                CompSampleNames <- colnames(mergedIntersect)
            
            } else {
                mergedIntersect2add <- mergedExprs[,platformStart:platformEnd]
                mergedReprod2add <- mergedExprsReprod[,platformStart:platformEnd]
            
                ##### Select only samples included in current comparison
                mergedIntersect2add <- mergedIntersect2add[,unlist(phenoData(merged.part)[[j]][1]) %in% CompNos]
                mergedReprod2add <- mergedReprod2add[,unlist(phenoData(merged.part)[[j]][1]) %in% CompNos]
            
                mergedIntersect <- mergeExprs(mergedIntersect,mergedIntersect2add)
                mergedReprod <- mergeExprs(mergedReprod,mergedReprod2add)
            
                esets[[platformsUsedNo]] <- new("ExpressionSet", exprs = mergedIntersect2add, phenoData = new("AnnotatedDataFrame", data=phenoData(merged.part)[[j]])[unlist(phenoData(merged.part)[[j]][1]) %in% CompNos,] )
                esetsReprod[[platformsUsedNo]] <- new("ExpressionSet", exprs = mergedReprod2add, phenoData = new("AnnotatedDataFrame", data=phenoData(merged.part)[[j]])[unlist(phenoData(merged.part)[[j]][1]) %in% CompNos,] )
                
                CompSampleNames <- c(CompSampleNames, colnames(mergedIntersect2add))
            }
            
            ##### Create target for datasets considered in current comparison
            targetsComp <- c(targetsComp, targets[platformStart:platformEnd])
        }
    }
    
    ##### Skip comparisons for where one gorup is present in only one type of platform
    if ( length(esets) < 2 ) {
    
        cat(paste("\nThe between and within study variability calculations are skiped for comparison", CompName, " as one of the gorups is present in only one type of high-coverage platform!\n\n",sep=" "))
        
    } else {
    
        features <- list(featureNames(esets),featureNames(esetsReprod))
        
        platformNamesUsed <- platformNames.part[platformsUsed]
        
        ##### Annotate instersected phenoData with targets
        platformDimIntersect <- summary(mergedIntersect)$`Number of Samples in Each Study`[2,]
    
        targetsComp <- targetsComp[targetsComp %in% CompNos]
    
        platformStart = 0
        platformEnd = 0
    
        ##### Go through remaining platforms
        for (j in 1:platformsUsedNo) {
    
            platformStart = platformEnd + 1
            platformEnd = platformStart + as.numeric(platformDimIntersect[j]) - 1
    
            phenoData(mergedIntersect)[[j]][1] <- as.factor(targetsComp[platformStart:platformEnd])
        }

        ##### Compute the unbiased estimates of the effect and its variance for entire set of genes and reproducible genes only
        d.merged <- list()
        d.adj.d.merged <- list()
        var.d.adj.merged <- list()
        outcome <- list()
        mns <- list(matrix(NA,length(featureNames(esets)),platformsUsedNo), matrix(NA,length(featureNames(esetsReprod)),platformsUsedNo))
        vars <- list(matrix(NA,length(featureNames(esets)),platformsUsedNo), matrix(NA,length(featureNames(esetsReprod)), platformsUsedNo))
        tau2.DL <- list()
    
        pdf(paste("Meta",CompName,"Xstudy_variance.pdf",sep="_"), height=8, width=22)
        par(mfrow = c( 2, 6 ))
    
        plotTitle <- NULL
    
        for (j in 1:2) {
            for (k in 1:platformsUsedNo) {
    
                outcome[[k]] <- as.numeric(unlist(phenoData(mergedIntersect)[[k]][1]))-1
            
                if (j==1) {
                    d.merged[[k]] <- getdF(mergedIntersect[k], outcome[[k]])
                    plotTitle[j] <- paste("(All ", length(featureNames(esets)), " genes)", sep="")
                
                } else if (j==2) {
                    d.merged[[k]] <- getdF(mergedReprod[k], outcome[[k]])
                    plotTitle[j] <- paste("(", length(featureNames(esetsReprod)), " reproducible genes)", sep="")
                }
            
                d.adj.d.merged[[k]] <- dstar(d.merged[[k]], as.numeric(summary(mergedIntersect)[[2]][2,k]))
                var.d.adj.merged[[k]] <- sigmad(d.adj.d.merged[[k]], sum(outcome[[k]]==0), sum(outcome[[k]]==1))
     
                mns[[j]][,k] <- d.adj.d.merged[[k]]
                vars[[j]][,k] <- var.d.adj.merged[[k]]
            }
        
            rownames(mns[[j]]) <- features[[j]]
            rownames(vars[[j]]) <- features[[j]]
        
            ##### Consider only genes with no missing values
            mns[[j]] <- mns[[j]][complete.cases(mns[[j]]), ]
            vars[[j]] <- vars[[j]][complete.cases(vars[[j]]), ]
        
            ##### Compute Cochran's Q statistic
            Q.values <- f.Q(mns[[j]], vars[[j]])
        
            ##### Display histogram of Q-values
            hist(Q.values,breaks=50,col="red", xlab="Q-values", main=paste("Histogram of Q-values", plotTitle[j], sep=" "))
            legend("topright",legend=paste("Q-values mean:", round(mean(Q.values), digits = 2), sep=" "), box.col="transparent")

            ##### Create and display a qq-plot for comparing the observed Q-values to a Chi square random variable
            chisqq <- qchisq(seq(0, .9999, .001), df=platformsUsedNo-1)
            tmp<-quantile(Q.values, seq(0, .9999, .001))
        
            qqplot(chisqq, tmp, ylab="Quantiles of Q-values", pch=16, xlab="Quantiles of Chi square", main=paste("QQ plot", plotTitle[j], sep=" "))
            lines(chisqq, chisqq, lty="dotted",col="red")
    
            ##### Calculate the fixed-effects model
            muFEM = mu.tau2(mns[[j]], vars[[j]])
            sdFEM = var.tau2(vars[[j]])
            ZFEM = muFEM/sqrt(sdFEM)
    
            ##### Plot the quantiles of the FEM with expected quantiles
            qqnorm(ZFEM, pch=16, ylab="Quantiles of FEM estiamtes", main=paste("Normal QQ plot", plotTitle[j], sep=" "))
            qqline(ZFEM, col="red")

            ##### Estimate the variance tau of the 'between experiments' random variable
            tau2.DL[[j]]<-tau2.DL(Q.values, length(mergedIntersect), my.weights=1/vars[[j]])
            ##### Obtain new variances s^2+tau^2
            varsDL <- vars[[j]] + tau2.DL[[j]]
            muREM <- mu.tau2(mns[[j]], varsDL)
            ##### Cumpute mu(tau)
            varREM <- var.tau2(varsDL)
            ZREM <- muREM/sqrt(varREM)
    
            ##### Plot the quantiles of the REM with expected quantiles
            qqnorm(ZREM, pch=16, ylab="Quantiles of REM estiamtes",, main=paste("Normal QQ plot", plotTitle[j], sep=" "))
            qqline(ZREM, col="red")
    
            ##### Plot the two estimates (fixed- and random effects model). If there is not much difference it's because in the REM model for most of the genes the variance tau is estimated as zero
            plot(muFEM, muREM, pch=".", main=plotTitle[j], xlab="FEM estimates", ylab="REM estimates")
            abline(0,1,col="red")

            ##### Display histogram of tau values
            hist(tau2.DL[[j]],col="red", breaks=50, xlab="Tau", main=paste("Histogram of tau", plotTitle[j], sep=" "))
            legend("topright",legend=paste("Tau mean:", round(mean(tau2.DL[[j]]), digits = 2), sep=" "), box.col="transparent")
        }
        dev.off()
    
        ##### Comapre the effect sizes between single and combined datasets, as decribed in the paper by Choi et al., 2003 ('Combining multiple microarray studies and modeling interstudy variation')
        ##### Calculate z-scores
        classes <- outcome
        classes[[platformsUsedNo+1]] <- unlist(classes)
        theScores <- list()

        esets[[platformsUsedNo+1]] <- new("ExpressionSet", exprs = mergedExprs[,CompSampleNames], phenoData = new("AnnotatedDataFrame", data=as.data.frame(do.call(rbind, phenoData(merged))[CompSampleNames,], row.names=CompSampleNames)) )
        theScores[[1]] <- zScores(esets,classes,useREM=FALSE)

        esetsReprod[[platformsUsedNo+1]] <- new("ExpressionSet", exprs = mergedExprsReprod[,CompSampleNames], phenoData = new("AnnotatedDataFrame", data=as.data.frame(do.call(rbind, phenoData(merged))[CompSampleNames,], row.names=CompSampleNames)) )
        theScores[[2]] <- zScores(esetsReprod,classes,useREM=FALSE)
    
        theScores <- list.rm.na(theScores)
    
        esetsList <- list(esets,esetsReprod)
    
        Comb <- combn(c(1:platformsUsedNo), 2)
    
        ##### Plot the z-scores the individual dataset against the z-scores of the combined dataset
        pdf(paste("Meta",CompName,"Zscores_plot.pdf",sep="_"), height=3.5, width=3*platformsUsedNo, pointsize=12)
        par(mfrow = c( 1, platformsUsedNo ))
    
        for (j in 1:2) {
            for (k in 1:platformsUsedNo) {
                plot(theScores[[j]][,k],theScores[[j]][,"zSco"],pch=".",xlab=paste("Original z-score (", platformNamesUsed[k], ")", sep=""), ylab="Combined z-score", main=plotTitle[j])
            }
        }
        dev.off()
    
        ##### Compute integration-driven discovery rates (IDRs) and generate IDR plot, as decribed in the paper by Choi et al., 2003 ('Combining multiple microarray studies and modeling interstudy variation'). For a threshold z-th this plot shows the fraction of the genes that have a higher effct size than the threshold for the combined effect z, but not for any of the experiment specific effects z-i
        pdf(paste("Meta",CompName,"IDRplot.pdf",sep="_"), height=3, width=3*(ncol(Comb)+1), pointsize=12)
        
        for (j in 1:2) {
            if (ncol(Comb) < 2) {
            
                IDRplot(theScores[[j]],Combine=1:length(mergedIntersect),colPos="red", colNeg="blue", xlab=paste("Z-threshold (", platformNamesUsed[Comb[1,1]], " + ", platformNamesUsed[Comb[2,1]], ")", sep=""), main=plotTitle[j])
            
            } else {
                par(mfrow = c( 1, ncol(Comb)+1 ))
            
                ##### Go through all possible comninations of the platforms
                for(k in 1:ncol(Comb)) {
                    theScoresTemp <- zScores(esetsList[[j]],classes,useREM=FALSE, CombineExp=c(Comb[1,k],Comb[2,k]))
                    theScoresTemp <- theScoresTemp[complete.cases(theScoresTemp),]
                    IDRplot(theScoresTemp,Combine=c(Comb[1,k],Comb[2,k]),colPos="red", colNeg="blue", xlab=paste("Z-threshold (", platformNamesUsed[Comb[1,k]], " + ", platformNamesUsed[Comb[2,k]], ")", sep=""), main=plotTitle[j])
                }
                ##### Include all platfroms
                IDRplot(theScores[[j]],Combine=1:length(mergedIntersect),colPos="red", colNeg="blue", xlab="Z-threshold (all datasets)", main=plotTitle[j])
            }
        }
        dev.off()
        
    
        ##### Estimating the false discovery rate (FDR), as decribed in the paper by Choi et al., 2003 ('Combining multiple microarray studies and modeling interstudy variation')
        pdf(paste("Meta",CompName,"FDRplot.pdf",sep="_"), height=3, width=3*(ncol(Comb)+1), pointsize=12)
    
        for (j in 1:2) {
            if  (ncol(Comb) < 2) {
        
                ScoresFDR <- zScoreFDR(esetsList[[j]], classes, useREM=FALSE, nperm=50, CombineExp=Comb)
                FDRwholeSettwo <- sort(ScoresFDR$"two.sided"[,"FDR"])
            
                experimentstwo <- list()
                for(k in 1:2){
                    experimentstwo[[k]] <- sort(ScoresFDR$"two.sided"[,paste("FDR_Ex_",Comb[k,1],sep="")])
                }
            
                ##### Two sided z-values
                plot(FDRwholeSettwo,pch="*",col=platforms.colour.part["Combined"],ylab="FDR",xlab="Number of genes", main=plotTitle[j])
                for (k in 1:2) {
                    points(experimentstwo[[k]], pch="*", col=platforms.colour.part[platformNamesUsed[k]])
                }
                legend("topleft",c(platformNamesUsed,"Combined datasets"), fill=platforms.colour.part[c(platformNamesUsed,"Combined")], box.col="transparent")
    
            } else {
                par(mfrow = c( 1, ncol(Comb)+1 ))
            
                ##### Go through all possible comninations of two platforms
                for(k in 1:ncol(Comb)) {
        
                    ScoresFDR <- zScoreFDR(esetsList[[j]], classes, useREM=FALSE, nperm=50, CombineExp=c(Comb[1,k],Comb[2,k]))
                    FDRwholeSettwo <- sort(ScoresFDR$"two.sided"[,"FDR"])
        
                    experimentstwo <- list()
                    for(l in 1:2){
                        experimentstwo[[l]] <- sort(ScoresFDR$"two.sided"[,paste("FDR_Ex_",Comb[l,k],sep="")])
                    }
    
                    ##### Two sided z-values
                    plot(FDRwholeSettwo,pch="*",col=platforms.colour.part["Combined"],ylab="FDR",xlab="Number of genes", main=plotTitle[j])
                    for (l in 1:2) {
                        points(experimentstwo[[l]], pch="*", col=platforms.colour.part[Comb[l,k]])
                    }
                    legend("topleft",c(platformNamesUsed[c(Comb[1,k],Comb[2,k])],"Combined datasets"), fill=platforms.colour.part[c(platformNamesUsed[c(Comb[1,k],Comb[2,k])],"Combined")], box.col="transparent")
                }
    
                ##### Include all platfroms
                ScoresFDR <- zScoreFDR(esetsList[[j]], classes, useREM=FALSE, nperm=50, CombineExp=1:length(mergedIntersect))
                FDRwholeSettwo <- sort(ScoresFDR$"two.sided"[,"FDR"])

                experimentstwo <- list()
                for(l in 1:length(mergedIntersect)){
                    experimentstwo[[l]] <- sort(ScoresFDR$"two.sided"[,paste("FDR_Ex_",l,sep="")])
                }
    
                ##### Two sided z-values
                plot(FDRwholeSettwo,pch="*",col=platforms.colour.part["Combined"],ylab="FDR",xlab="Number of genes", main=plotTitle[j])
                for (l in 1:length(mergedIntersect)) {
                    points(experimentstwo[[l]], pch="*", col=platforms.colour.part[l])
                }
                legend("topleft",c(platformNamesUsed, "Combined datasets"), fill=platforms.colour.part, box.col="transparent")
            }
        }
        dev.off()
    
    
        ##### Generage count plot with number of gene that are below a given threshold for the FDR. For each study (indicated by different colours) and various thresholds for the FDR (x axis) it illutrates the number of genes (y axis) that are below this threshold in the given study but above in all other studies
        pdf(paste("Meta",CompName,"FDRcount_plot.pdf",sep="_"), height=3, width=3*(ncol(Comb)+1), pointsize=12)
    
        for (j in 1:2) {
            if  (ncol(Comb) < 2) {
            
                ScoresFDR <- zScoreFDR(esetsList[[j]], classes, useREM=FALSE, nperm=50, CombineExp=Comb)
                ScoresFDR <- list.rm.na(ScoresFDR)
            
                CountPlot(ScoresFDR, Score="FDR", kindof="two.sided", cols=platforms.colour.part[c(platformNamesUsed,"Combined")], main=plotTitle[j], CombineExp=Comb)
                legend("topleft",legend=c(platformNamesUsed,"Combined datasets"), fill=platforms.colour.part[c(platformNamesUsed,"Combined")], box.col="transparent")
    
            } else {
                par(mfrow = c( 1, ncol(Comb)+1 ))
        
                ##### Go through all possible comninations of two platforms
                for(k in 1:ncol(Comb)) {
            
                    ScoresFDR <- zScoreFDR(esetsList[[j]], classes, useREM=FALSE, nperm=50, CombineExp=c(Comb[1,k],Comb[2,k]))
                    ScoresFDR <- list.rm.na(ScoresFDR)
                
                    CountPlot(ScoresFDR, Score="FDR", kindof="two.sided", cols=platforms.colour.part[c(platformNamesUsed[c(Comb[1,k],Comb[2,k])],"Combined")], main=plotTitle[j], CombineExp=c(Comb[1,k],Comb[2,k]))
                    legend("topleft",legend=c(platformNamesUsed[c(Comb[1,k],Comb[2,k])],"Combined datasets"), fill=platforms.colour.part[c(platformNamesUsed[c(Comb[1,k],Comb[2,k])],"Combined")], box.col="transparent")
                }
                ##### Include all platfroms
                ScoresFDR <- zScoreFDR(esetsList[[j]], classes, useREM=FALSE, nperm=50, CombineExp=1:length(mergedIntersect))
                ScoresFDR <- list.rm.na(ScoresFDR)
            
                CountPlot(ScoresFDR, Score="FDR", kindof="two.sided", cols=platforms.colour.part, main=plotTitle[j], CombineExp=1:length(mergedIntersect))
                legend("topleft",legend=c(platformNamesUsed,"Combined datasets"), fill=platforms.colour.part, box.col="transparent")
            }
        }
        dev.off()
    }
    
    #===============================================================================
    #     Meta-analysis using Stouffer's method ( formulas taken from http://www.broadinstitute.org/~debakker/meta.R )
    #===============================================================================
    ##### Draw plot of the combined p-values as a function of the weigthing
    ##### Calculate the combined p-values accordingly to the following formula:
    #
    # First convert p-values to z-scores
    #
    # Zcomb = sum(sqrt(w) * Z) / sqrt(sum(w))
    #
    # Z = initial z-score
    # w = weight
    #
    # Finally convert combined z-score into p-value
    #
    ####################################################################
    
    cat("\nPerforming meta-analysis...\n\n")
        
    platformStart <- 0
    platformEnd <- 0
        
    platformsUsedNo <- 0
    platformsUsed <- NULL
        
    ##### Report platforms which contain samples from both groups to be compared
    for (j in 1:platforms_No) {
            
        platformStart = platformEnd + 1
        platformEnd = platformStart + as.numeric(platformDim[j]) -1
            
        ##### Don't include dataset which doesn't contain samples from both groups to be compared
        if ( length(intersect(targets[platformStart:platformEnd], CompNos))==2  ) {
                
            platformsUsed <- c(platformsUsed, j)
                
            }
    }
        
    platformNamesUsed <- platformNames[platformsUsed]
        
    ##### Read in differential expression analysis results for each platform
    DEresults <- list()
    platformsUsedNo <- 0
        
    for ( j in seq(1:platforms_No) ) {
        
        if ( unique(factor(DatasetInput[,"Platform"]))[j] %in% platformNamesUsed ) {
                
            platformsUsedNo <- platformsUsedNo + 1
                
            cat(paste("Uploading differential expression analysis results for platform ", unique(factor(DatasetInput[,"Platform"]))[j] ,"...\n", sep=""))
            DEresults[[platformsUsedNo]] <- read.table(sub(".exp", paste("_", CompName, "_topTable.txt", sep=""), dataset.expFiles[[j]]),sep="\t",as.is=TRUE,header=TRUE,row.names=1,check.names=FALSE,quote = "")
        }
    }

    ##### Skip comparisons for where one gorup is present in only one type of platform
    if ( length(DEresults) < 2 ) {
    
        cat(paste("\nThe comparison", CompName, "is skipped as the direct groups comparison is present within only one type of platform!\n\n",sep=" "))
        
    } else {
        
        ##### Retrieve all genes
        genes <- NULL

        for ( j in 1:platformsUsedNo ) {
    
            genes <- union(genes, rownames(DEresults[[j]]) )
        }
    
        ##### Remove genes with negative integrative correlation values
        genes <- setdiff( genes, names(integrativeCors[integrativeCors < 0]) )
    
        Pval <- sapply(genes,function(x) NULL)
        FC <- sapply(genes,function(x) NULL)
        Pval2rep <- sapply(genes,function(x) NULL)
        FC2rep <- sapply(genes,function(x) NULL)
        Pval.perPlatform <- sapply(platformNamesUsed,function(x) NULL)
        FC.perPlatform <- sapply(platformNamesUsed,function(x) NULL)
        intCor.PerPlatform <- sapply(platformNamesUsed,function(x) NULL)
        tau.PerPlatform <- sapply(platformNamesUsed,function(x) NULL)
        ZvalIntersect <- NULL
    
        ##### Retrieve p-values and fold-changes for each platform
        for (j in 1:length(genes)) {
            for ( k in 1:platformsUsedNo ) {
        
                Pval2rep[[genes[j]]] <- c(Pval2rep[[genes[j]]], DEresults[[k]][genes[j],11])
                Pval[[genes[j]]] <- Pval2rep[[genes[j]]][!is.na(Pval2rep[[genes[j]]])]
                FC2rep[[genes[j]]] <- c(FC2rep[[genes[j]]], DEresults[[k]][genes[j],8] )
                FC[[genes[j]]] <- FC2rep[[genes[j]]][!is.na(FC2rep[[genes[j]]])]
                Pval.perPlatform[[k]] <- c( Pval.perPlatform[[k]], DEresults[[k]][genes[j],11] )
                FC.perPlatform[[k]] <- c( FC.perPlatform[[k]], DEresults[[k]][genes[j],8] )
                intCor.PerPlatform[[k]] <- c(intCor.PerPlatform[[k]], integrativeCors[genes[j]] )
                
                if ( platforms_No.part > 1 ) {
                    
                    tau.PerPlatform[[k]] <- c(tau.PerPlatform[[k]], tau2.DL[[1]][genes[j]] )
                }
            }
            ZvalIntersect <- c( ZvalIntersect, convert.pvalue( Pval[[genes[j]]], FC[[genes[j]]]) )
        }
    
        for ( j in 1:platformsUsedNo ) {
        
            names(Pval.perPlatform[[j]]) <- names(Pval)
            names(FC.perPlatform[[j]]) <- names(FC)
        
            Pval.perPlatform[[j]] <- Pval.perPlatform[[j]][!is.na(Pval.perPlatform[[j]])]
            FC.perPlatform[[j]] <- FC.perPlatform[[j]][!is.na(FC.perPlatform[[j]])]
            intCor.PerPlatform[[j]] <- intCor.PerPlatform[[j]][!is.na(intCor.PerPlatform[[j]])]
            
            ##### Perfrom this part only if there are > 1 high-coverage platforms
            if ( platforms_No.part > 1 ) {
                tau.PerPlatform[[j]] <- tau.PerPlatform[[j]][!is.na(tau.PerPlatform[[j]])]
            }
        }
        
        ##### Perfrom this part only if there are > 1 high-coverage platforms
        if ( platforms_No.part > 1 ) {
            tau <- unlist(tau2.DL[[1]])
        }

            
        pdf(paste("Meta",CompName,"FC_intCor_tau_hist.pdf",sep="_"), height=4, width=12)
        par(mfrow = c( 1, 3 ))
        hist(abs(unlist(FC.perPlatform)), breaks = seq(min(abs(unlist(FC.perPlatform))),max(abs(unlist(FC.perPlatform)))+0.1,0.1), main="Histogram of fold-changes", xlab="Log2 fold-change", col="bisque3")
        hist(integrativeCors[genes], breaks = seq(min(integrativeCors[genes], na.rm=TRUE),max(integrativeCors[genes], na.rm=TRUE)+0.02,0.02), main="Histogram of integrative correlation coefficients", xlab="Integrative correlation coefficient", col="olivedrab")
            
        ##### Perfrom this part only if there are > 1 high-coverage platforms
        if ( platforms_No.part > 1 ) {
                
            tauHist <- hist(tau, breaks = seq(min(tau),max(tau)+0.1,0.1), xlim=c(0,3), main="Histogram of tau values", xlab="Tau", col="lightcoral")
        }
        dev.off()
    
        ##### Plot FC (log2) obtained in each platform with p-values, intagrative correlation and tau values (add loess (local polynomial regression fitting) smoothing line)
        pdf(paste("Meta",CompName,"FC2Pvalue2intCor2tau_plot.pdf",sep="_"), height=2.7*platformsUsedNo, width=7.5)
        par(mfrow = c( platformsUsedNo, 3 ))
    
        for (j in 1:(platformsUsedNo) ) {
        
            Pval.vec <- Pval.perPlatform[[j]]
            FC.vec <- FC.perPlatform[[j]]
            intCor.vec <- intCor.PerPlatform[[j]][names(FC.vec)]
            intCor.vec <- intCor.vec[!is.na(intCor.vec)]
            
            Pval.order <- order(Pval.vec, decreasing = TRUE)
            intCor.order <- order(intCor.vec, decreasing = TRUE)
        
            ##### ...p-values
            plot(-log2(Pval.vec[Pval.order]), abs(FC.vec[Pval.order]), pch=16, cex = 0.5, col=platforms.colour[platformNamesUsed[j]],ylab="Log2 fold-change", xlab="-log2(p-value) index", main="")
            loess <- loessSmoothing(-log2(Pval.vec[Pval.order]), abs(FC.vec[Pval.order]))
            if ( !any(is.na(loess[[1]][2]$fitted)) ) {
                lines(abs(loess[[2]]), predict(loess[[1]],loess[[2]]), col="red", lwd=2)
            }
        
            ##### ...intagrative correlation coefficients
            plot(intCor.vec[intCor.order], abs(FC.vec[intCor.order]), pch=16, cex = 0.5, col=platforms.colour[platformNamesUsed[j]],ylab="Log2 fold-change", xlab="Integrative correlation index", main=platformNamesUsed[j])
            loess <- loessSmoothing(intCor.vec[intCor.order], abs(FC.vec[intCor.order]))
            lines(loess[[2]], predict(loess[[1]],loess[[2]]), col="red", lwd=2)
            
            ##### Perfrom this part only if there are > 1 high-coverage platforms
            if ( platforms_No.part > 1 ) {
                
                tau.vec <- tau.PerPlatform[[j]][names(FC.vec)]
                tau.vec <- tau.vec[!is.na(tau.vec)]
                
                ##### Set tau value with frequency > 10 (genes) as threshold for plotting
                tauThreshold <- max(tauHist$breaks[which(tauHist$counts > 50)])
                tau2plot <- tau.vec[tau.vec<tauThreshold]
                tau.order <- order(tau2plot, decreasing = TRUE)
        
                ##### ...tau values
                plot(tau2plot[tau.order], abs(FC.vec[tau.order]), pch=16, cex = 0.5, col=platforms.colour[platformNamesUsed[j]],ylab="Log2 fold-change", xlab="Tau index", main="")
                loess <- loessSmoothing(tau2plot[tau.order], abs(FC.vec[tau.order]))
                lines(loess[[2]], predict(loess[[1]],loess[[2]]), col="red", lwd=2)
                
            } else {
                plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
            }
        }
        dev.off()


        ##### Draw plot of the weighting as a function of log2 fold-change and integrative correlation coefficients. Assume there is no significant variance between studies (tau ~ 0)
        ##### Calculate the weights accordingly to the following formula:
        #
        # Weight = |FC| + (|FC| * IntCor^2)
        #
        # FC = log2 fold-change
        # IntCorr = integrative correlation
        #
        ####################################################################
    
        weigth <- sapply(seq(1,6,1),function(x) NULL)
    
        ##### Consider integrative correlation coefficients "0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9" and "1"
        for (j in seq(0,1,0.1)) {
            ##### ... and log2 fold-change values "0", "1", "2", "3", "4" and "5"
            for (k in seq(0,5,1)) {
                weigth[[k+1]] = c(weigth[[k+1]], k + (k * j^2))
            }
        }
    
        pdf(paste("Meta_weighting_plot.pdf",sep="_"), height=6, width=6)
        plot(seq(0,1,0.1), weigth[[1]], type="l", lwd=5, ylab="Weighting for Stouffer's method ", xlab="Integrative correlation coefficient", main="", ylim=c(0, max(unlist(weigth))), col=terrain.colors(13)[2])
        for (j in seq(2,6,1)) {
        
            points(seq(0,1,0.1), weigth[[j]], type="l", lwd=5, col=terrain.colors(13)[j*2])
        }
        grid(nx=NA, ny=NULL , lwd = 2)
        legend("topleft",legend=c(paste("log2 FC =", seq(0,5,1),sep=" ")), fill=terrain.colors(13)[seq(2,13,2)], box.col="transparent")
        dev.off()
    

        ##### Compute the combined z-score using weighted Stouffer's method
        Zcomb <- sapply(genes,function(x) NULL)
        Pcomb <- sapply(genes,function(x) NULL)
    
        for (j in 1:length(genes)) {
        
            ##### Report only genes present in > 1 platform
            if ( length(FC[[genes[j]]])>1 ) {
                ##### Compute weights for each study
                weigths <- weighting(FC[[genes[j]]], integrativeCors[genes[j]])
        
                ##### Compute the comnbined z-score
                Zcomb[[genes[j]]] <- p.to.z(Pval[[genes[j]]], FC[[genes[j]]],  weigths)
        
                ##### Convert the combined z-scores into p-value
                Pcomb[[genes[j]]] <- pnorm(-(abs(Zcomb[[genes[j]]]))) * 2
            }
        }
        genesComb <- names(unlist(Pcomb))

        ##### Generate histograms of combined p-values and z-scores
        pdf(paste("Meta",CompName,"Pvalue_Zscore_hist.pdf",sep="_"), height=8, width=12)
        par(mfrow = c( 2, 3 ))
    
        ##### ...p-values
        hist(unlist(Pval.perPlatform), breaks = seq(min(unlist(Pval.perPlatform)),max(unlist(Pval.perPlatform))+0.02,0.02), main="Histogram of p-values", xlab="p-value", col="slategray")
        hist(unlist(Pcomb), breaks = seq(min(unlist(Pcomb)),max(unlist(Pcomb))+0.02,0.02), main="Histogram of combined p-values", xlab="Combined p-value", col="lightcoral")

        ##### Draw density plots
        density1 <- density(unlist(Pval.perPlatform))
        density2 <- density(unlist(Pcomb))
        maxY <- max(density1$y,density2$y)

        plot(density1, xlim=c(0,1), ylim=c(0,maxY), xlab="", main="Density plot of p-values")
        par(new=TRUE)
        plot(density2, xlim=c(0,1), ylim=c(0,maxY), main="", xlab="", ylab="", xaxt='n', yaxt='n', col="lightcoral")
        legend(0,maxY, legend=c(" p-values"," combined p-values"),col=c("black","lightcoral"), lty=1, box.col="transparent")

        ##### ...z-scores
        hist(ZvalIntersect, breaks = seq(min(ZvalIntersect),max(ZvalIntersect)+1,1), main="Histogram of z-scores", xlab="Z-score", col="slategray")
        hist(unlist(Zcomb), breaks = seq(min(unlist(Zcomb)),max(unlist(Zcomb))+1,1), main="Histogram of combined z-scores", xlab="Combined z-score", col="lightcoral")

        ##### Draw density plots
        density1 <- density(ZvalIntersect)
        density2 <- density(unlist(Zcomb))
        minX <- min(density1$x,density2$x)
        maxX <- max(density1$x,density2$x)
        maxY <- max(density1$y,density2$y)

        plot(density1, xlim=c(minX,maxX), ylim=c(0,maxY), xlab="", main="Density plot of z-scores")
        par(new=TRUE)
        plot(density2, xlim=c(minX,maxX), ylim=c(0,maxY), main="", xlab="", ylab="", xaxt='n', yaxt='n', col="lightcoral")
        legend(minX,maxY, legend=c(" z-scores"," combined z-scores"),col=c("black","lightcoral"), lty=1, box.col="transparent")
        dev.off()


        ##### Draw a plots illustrating the p-value improvement after meta-analysis for each platform
        PvalRatio <- sapply(platformNames,function(x) NULL)
        PvalRatioKnownGenes <- sapply(platformNames,function(x) NULL)
    
        for (j in 1:platformsUsedNo) {
            for (k in 1:length(genesComb)) {
            
                PvalRatio[[j]] <- c(PvalRatio[[j]], log2(Pval2rep[[genesComb[k]]][j]/Pcomb[[genesComb[k]]]))
            
                ##### Report the p-value improvement for genes assocaited with studied phenotype
                if ( !is.na(knownGenes) && genesComb[k] %in% knownGenes[,1] ) {
                
                    PvalRatioKnownGenes[[j]] <- c(PvalRatioKnownGenes[[j]], log2(Pval2rep[[genesComb[k]]][j]/Pcomb[[genesComb[k]]]))
                } else {
                    PvalRatioKnownGenes[[j]] <- c(PvalRatioKnownGenes[[j]], NA)
                }
            }
        }
    
        pdf(paste("Meta",CompName,"PvalueRatio.pdf",sep="_"), height=6, width=4*platformsUsedNo, pointsize=12)
        par(mfrow = c( 1, platformsUsedNo ))
    
        for (j in 1:(platformsUsedNo) ) {
            plot(seq(0,length(genesComb)-1,1), PvalRatio[[j]][order(unlist(Pcomb), decreasing = TRUE)], pch=16, cex = 0.5, col=platforms.colour[platformNamesUsed[j]],ylab="p-value / combined p-value (log2)", xlab="Combined p-value index", main=platformNamesUsed[j], ylim=c(min(unlist(PvalRatio)[!is.na(unlist(PvalRatio))]), max(unlist(PvalRatio)[!is.na(unlist(PvalRatio))])))
        
            ##### Report the p-value improvement for genes assocaited with studied phenotype
            points(seq(0,length(genesComb)-1,1), PvalRatioKnownGenes[[j]][order(unlist(Pcomb), decreasing = TRUE)], pch=16, col="red", cex=0.8)
            points(seq(0,length(genesComb)-1,1), PvalRatioKnownGenes[[j]][order(unlist(Pcomb), decreasing = TRUE)], pch=16, cex=0.65)
            abline(h=0,col="black")
            abline(v=length(genesComb)-length(which(sort(unlist(Pcomb), decreasing = TRUE)<1e-09)),col="black", lty = 2)
            legend("topleft", legend=" Combined p-value = 1e-09", col="black", lty=2, box.col="transparent")
        }
        dev.off()
    
        pdf(paste("Meta",CompName,"Pcomb2Pvalue.pdf",sep="_"), height=6, width=12)
        ##### Use the actual p-values
        plot(seq(0,length(genesComb)-1,1), -log2(unlist(Pcomb))[order(unlist(Pcomb), decreasing = TRUE)], type="n", xlab="", ylab="", ylim=c(min(-log2(unlist(Pval))), max(-log2(unlist(Pval)))))
        for (j in 1:(platformsUsedNo) ) {
            points(seq(0,length(genesComb)-1,1), -log2(Pval.perPlatform[[j]][genesComb])[order(unlist(Pcomb), decreasing = TRUE)], type="l", col=platforms.colour[platformNamesUsed[j]], lwd=2)
        }
    
        par(new=TRUE)
        plot(seq(0,length(genesComb)-1,1), -log2(unlist(Pcomb))[order(unlist(Pcomb), decreasing = TRUE)], col=platforms.colour[length(platforms.colour)],ylab="-log2(p-value)", xlab="gene index", main="", ylim=c(min(-log2(unlist(Pval))), max(-log2(unlist(Pval)))), type="l", lwd=3)
        abline(h=-log2(1e-09),col="black", lty = 2)
        
        legend("topleft", legend=c(rep("",platformsUsedNo+2),"p-value = 1e-09"), col=c(rep("transparent",platformsUsedNo+2),"black"), lty=2, box.col="transparent")
        legend("topleft",legend=c(platformNamesUsed,"Combined datasets"), fill=platforms.colour[c(platformNamesUsed, "Combined")], box.col="transparent")
    
        ##### ... and loess smoothed values
        plot(seq(0,length(genesComb)-1,1), -log2(unlist(Pcomb))[order(unlist(Pcomb), decreasing = TRUE)], col=platforms.colour[length(platforms.colour)],ylab="-log2(p-value)", xlab="gene index (loess smoothing)", main="", ylim=c(min(-log2(unlist(Pval))), max(-log2(unlist(Pval)))), type="l", lwd=3)

        for (j in 1:(platformsUsedNo) ) {
            loess <- loessSmoothing(seq(0,length(genesComb)-1,1), -log2(Pval.perPlatform[[j]][genesComb])[order(unlist(Pcomb), decreasing = TRUE)])
            lines(loess[[2]], predict(loess[[1]],loess[[2]]), col=platforms.colour[platformNamesUsed[j]], lwd=2)
        }
        abline(h=-log2(1e-09),col="black", lty = 2)
        
        legend("topleft", legend=c(rep("",platformsUsedNo+2),"p-value = 1e-09"), col=c(rep("transparent",platformsUsedNo+2),"black"), lty=2, box.col="transparent")
        legend("topleft",legend=c(platformNamesUsed,"Combined datasets"), fill=platforms.colour[c(platformNamesUsed, "Combined")], box.col="transparent")
        dev.off()
    
    
        ##### Write combined results into a file
        cat(paste("\nWriting combined p-values to ", paste("Meta",CompName,sep="_"), ".txt\n\n",sep=""))
    
        MetaResults <- NULL
    
        for (j in 1:length(genes)) {
        
            ##### Report fold-changes, p-values and...
            MetaGene <- FC2rep[[genes[j]]]
            MetaGene <- c(MetaGene, Pval2rep[[genes[j]]])
            FC2rep[[genes[j]]]
            MetaGene <- c(MetaGene, mean(FC2rep[[genes[j]]], na.rm=TRUE))
        
            ##### ... combined p-values if available
            if (!is.null(Pcomb[[genes[j]]]) ) {
                MetaGene <- c(MetaGene, Pcomb[[genes[j]]])
            } else {
                MetaGene <- c(MetaGene, Pval[[genes[j]]])
            }
            MetaResults <- rbind(MetaResults, MetaGene)
        }
        rownames(MetaResults) <- genes
        colnames(MetaResults) <- c(paste("log2FC (", platformNamesUsed, ")", sep=""), paste("adj p-value (", platformNamesUsed, ")", sep=""), "Average log2FC", "Combined p-value")
    
    
        ##### Retrieve gene annotation information
        MetaResultsAnnot <- annotGenes(MetaResults)
       	
        write.table(prepare2write(MetaResultsAnnot), file=paste("Meta_",CompName, ".txt",sep=""),sep="\t", row.names=FALSE)
    
        ##### Perform GO enrichment analysis
        MetaResultsSorted <- sort(MetaResults[,ncol(MetaResults)])
        topPvalues <- which(MetaResultsSorted < quantile( MetaResultsSorted, 0.1))
        bottomPvalues <- which(MetaResultsSorted >= quantile( MetaResultsSorted, 0.1))
	
        SignGenes <- MetaResultsSorted[topPvalues]
        NonSignGenes <- MetaResultsSorted[bottomPvalues]
	
        allGenes <- c(rep(0, each=length(NonSignGenes)), rep(1, each=length(SignGenes)))
        names(allGenes) <- c(names(NonSignGenes), names(SignGenes))
	
        ##### Biological processes
        GOdata <- new("topGOdata", ontology = "BP", allGenes = allGenes, geneSel = function(p) p < 1e-2, description = CompName, annot = annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")

        resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        GOtable <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 500)

        write.table(prepare2write(GOtable), file=paste("Meta_",CompName, "_GO_BP.txt",sep=""),sep="\t", row.names=FALSE)
    
        ##### Molecular functions
        GOdata <- new("topGOdata", ontology = "MF", allGenes = allGenes, geneSel = function(p) p < 1e-2, description = CompName, annot = annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")

        resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        GOtable <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 500)

        write.table(prepare2write(GOtable), file=paste("Meta_",CompName, "_GO_MF.txt",sep=""),sep="\t", row.names=FALSE)
    
        ##### Cellular components
        GOdata <- new("topGOdata", ontology = "CC", allGenes = allGenes, geneSel = function(p) p < 1e-2, description = CompName, annot = annFUN.org, mapping="org.Hs.eg.db", ID="Ensembl")

        resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        GOtable <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 500)

        write.table(prepare2write(GOtable), file=paste("Meta_",CompName, "_GO_CC.txt",sep=""),sep="\t", row.names=FALSE)
    }
}

##### Remove R session for meta-analysis
#system("rm Meta.R")

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
