################################################################################
#
#   File name: GenExpressionComb_Affy_HuEx1ST.R
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
#	Description: Pipeline combining datasets produced using Affymetrix Human Exon 1.0 ST Array (gene level). The script reads in three a priori generated files: (1) list of datasets to be analysed; (2) file produced with QC_Affy_HuEx1ST.R which lists outliers identified with ArrayOutliers function and (3) '.flat' file produced with ArrayAnnot.pl which lists reliable probes annotated with Ensembl gene IDs and gene symbol. Pipeline executes 'multiGene2ProbeFilter.R' to remove multi gene to probe mappings
#   IMPORTANT: Pipeline same as 'GenExpressionComb_Affy_HuEx1ST.R' with additional detection above background (DABG) filter was applied to reduce the number of false positives expected due to high number of transcript clusters covered by this platform. (i) Transcript clusters and probesets not present in any of analysed groups, (ii) probesets that were absent in more than half of the samples of a group or transcript clusters with more than half of the probesets absent were removed prior to downstream analysis accordingly to procedure described by Rodrigo-Domingo et al. 2013, Brief in Bioinformatics
#
#	Command line use example: R --file=./GenExpressionComb_Affy_HuEx1ST.R --args "/scratch/jack/data/PhD/Transcriptomics_project"
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

library(aroma.affymetrix)
library(sva)

source("MultiGene2ProbeFilter.R")

#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

ProjectDir = args[4]
probesFlatFile = paste(paste(ProjectDir, collapse = "/"),"Affy_HuEx1ST.flat",sep="/")

##### Read file with datasets information
DatasetInput=read.table(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE,row.names=1)
DatasetInput=DatasetInput[DatasetInput[,"Platform"] %in% "Affy_HuEx1ST",]

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

cat("Loading datasets produced using Affymetrix GeneChip Human Exon 1.0 ST Array (gene level) into an AffyBatch object...\n")

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
targetFile = targetFile[setdiff(rownames(targetFile),as.vector(as.list(outlierList))),]

dirName = targetFile$FileName

##### Set up directory structure for pipeline based on input type
dir.create(paste("/scratch/jack/Aroma_R/rawData/", paste(studyIDs, collapse="_"), "/",sep=""))
dir.create(paste("/scratch/jack/Aroma_R/rawData/", paste(studyIDs, collapse="_"), "/HuEx-1_0-st-v2",sep=""))
		
cat(paste("Reading CEL files...\n"))
		
##### Copy file over from data directory to aroma directory
for (i in 1:nrow(targetFile)) {

	name = tail(unlist(strsplit(dirName[i], split='/', fixed=TRUE)), n=1)

	file.copy(dirName[i], paste("/scratch/jack/Aroma_R/rawData/", paste(studyIDs, collapse="_"),"/HuEx-1_0-st-v2/", name[i], sep=""))
	
	##### Change CEL file names to names specified in target file
	file.rename(paste("/scratch/jack/Aroma_R/rawData/", paste(studyIDs, collapse="_"),"/HuEx-1_0-st-v2/",name[i],sep=""), paste("/scratch/jack/Aroma_R/rawData/", paste(studyIDs, collapse="_"),"/HuEx-1_0-st-v2/",rownames(targetFile)[i], ".CEL", sep=""))
}

setwd("/scratch/jack/Aroma_R")

##### Get annotation files
tag = "coreR3,A20071112,EP"
annot_level = "coreR3"
chipType <- "HuEx-1_0-st-v2"
cdf <- AffymetrixCdfFile$byChipType(chipType,tags=tag)

##### Define CEL set
cs <- AffymetrixCelSet$byName(paste(studyIDs, collapse="_"), cdf=cdf)
setCdf(cs,cdf) 

setwd(paste("/scratch/jack/Aroma_R/rawData/", paste(studyIDs, collapse="_"),"/HuEx-1_0-st-v2/", sep=""))


################################################################################
##### Get Ensembl gene annotation for all probes
#probesAnnotFile = paste(paste(ProjectDir, collapse = "/"),"Affy_HuEx1ST_annot.txt",sep="/")
#Affy_HuEx1ST_probesAnnot <- read.table(probesAnnotFile,sep="\t",as.is=TRUE,header=TRUE)

##### Select an probes subset prior to normalization, based on the following post https://groups.google.com/forum/?hl=en#!topic/aroma-affymetrix/0zzGPFJPm_0

##### Indices of cells to be assigned NA
#cells2rm <-setdiff(unique(Affy_HuEx1ST_probesAnnot$cell),Affy_HuEx1ST_cells$cell) 
#signals <- rep(as.double(NA), length(cells2rm))

##### Write these signals to each CEL file 
#for (i in seq(cs)) { 
#	cf <- getFile(cs, i)
#	updateCel(getPathname(cf), indices=cells2rm, intensities=signals)
#} 
################################################################################

#===============================================================================
#    Background correction and quantile normalisation using aroma package
#===============================================================================

###computeAffinities(cdf, safe=TRUE, force=FALSE)

cat("Robust Multiarray Average normalisation with background correction.\n")

##### Background correction
bc <- RmaBackgroundCorrection(cs,tag=annot_level)
csBC <- process(bc)

##### Background correction taking into account GC content
bc <- GcRmaBackgroundCorrection(cs,tag=annot_level)
csBC <- process(bc)

##### Set up quantile normalization method
qn <- QuantileNormalization(csBC,typesToUpdate="pm")
csN <- process(qn)


#===============================================================================
#    Summarisation using "probe-level model" (PLM)
#===============================================================================

getCdf(csN)

##### Fit summary model transcript expression
##### Fit summary of the entire transcript (estimate overall expression for the transcript)
plmTr <- ExonRmaPlm(csN,mergeGroups=TRUE)

##### Run the model
fit(plmTr)

##### Generate data frame 
cesTr <- getChipEffectSet(plmTr)
trFit <- extractDataFrame(cesTr,units=NULL,addNames=TRUE)


##### Exon expression fit exon by exon
plmEx <- ExonRmaPlm(csN,mergeGroups=FALSE)

##### Run the model
fit(plmEx)

##### Generate a data frame
cesEx <- getChipEffectSet(plmEx)
ExFit <- extractDataFrame(cesEx, units=NULL,addNames=TRUE)


#===============================================================================
#    Intensity filtering
#===============================================================================

probesets.NetAffx.32 <- read.csv("/scratch/jack/Aroma_R/annotationData/chipTypes/HuEx-1_0-st-v2/HuEx-1_0-st-v2.na33.1.hg19.probeset.csv", skip=23)

##### Keep only columns of probesets.NetAffx.32 and transcript.clusters.NetAffx.32 that will be useful later
probesets.NetAffx.32 <- probesets.NetAffx.32[probesets.NetAffx.32$level == "core", c(1, 7, 12)]

targetFile$Number <- seq(1, nrow(targetFile))
groups <- split(targetFile$Number, targetFile$Target)

##### Keep only probes that do not cross-hybridise
non.crosshyb.probesets <- probesets.NetAffx.32$probeset_id[probesets.NetAffx.32$crosshyb_type == 1]

exon.intensities <- ExFit[ExFit$groupName %in% non.crosshyb.probesets, -c(3:5)]

##### Define absent/present probesets
##### Identify probesets for which half or more of probesets are present in a group
present.exons <- lapply(groups, FUN = function(group){
    if(length(group) > 1){
        apply(log2(exon.intensities[, group + 2]) < 3, 1, sum)/length(group) < 0.5
    } else {
        log2(exon.intensities[, group + 2]) > 3
    }
})

##### Convert to dataframe
present.exons <- t(do.call(rbind, present.exons))

##### Use probeset identities
rownames(present.exons) <- exon.intensities$groupName

##### Remove probesets not present in any of the groups
AllFalse <- function(vector){all(vector == FALSE)}
TwoOrMoreTrue <- function(vector){sum(vector == TRUE)>= 2}

absent.exons <- apply(present.exons, 1, AllFalse)
probesets.to.keep <- absent.exons[absent.exons == FALSE]
probesets.to.keep <- as.factor(names(probesets.to.keep))


##### Define absent/present transcripts
##### Identify transcripts for which half or more of probesets are present in sample
core.transcripts <- unique(exon.intensities$unitName)

present.genes <- lapply(core.transcripts, FUN = function(gene){
    
    apply(log2(exon.intensities[exon.intensities$unitName == gene, -c(1:2)]) < 3, 2, sum)/
    length(exon.intensities[exon.intensities$unitName == gene,]$groupName) < 0.5
}
)
names(present.genes) <- core.transcripts


##### Convert to logical matrix
present.genes <- do.call(rbind, present.genes)


##### Identify transcripts present in at least half of the samples of a group
present.genes.in.group <- lapply(groups, FUN = function(group){
    
    if(length(group) > 1){
        apply(present.genes[ , group], 1, sum)/length(group) >= 0.5
    } else {
        present.genes[ , group]
    }
}
)

##### Convert to logical matrix
present.genes.in.group <- do.call(rbind, present.genes.in.group)
present.genes.in.group <- t(present.genes.in.group)

##### Keep genes only present in at least two groups
transcripts.to.keep <- apply(present.genes.in.group, 1, TwoOrMoreTrue)
transcripts.to.keep <- names(transcripts.to.keep[transcripts.to.keep == TRUE])


probeset.intensities.dataframe.per.gene <- split(trFit, trFit$unitName)

probeset.intensities.dataframe.per.gene <- lapply(probeset.intensities.dataframe.per.gene, function(dataframe){
    n.probesets <- nrow(dataframe)
    data.frame(intensity = log2(c(as.matrix(dataframe[ , -c(1:5)]))),
    replicate = rep(targetFile$Replicate, each = n.probesets),
    treatment = rep(targetFile$Target, each = n.probesets),
    probeset = rep(dataframe$groupName)
    )
})

##### Keep only transcript clusters present in at least two groups
probeset.intensities.filtered <-  probeset.intensities.dataframe.per.gene[names(probeset.intensities.dataframe.per.gene) %in% transcripts.to.keep]


##### Remove absent probesets and transcripts with a single probeset to study
probeset.intensities.filtered <- lapply(probeset.intensities.filtered, function(dataframe){
    dataframe <- dataframe[dataframe$probeset %in% probesets.to.keep, ]
})

##### Remove transcript clusters that are empty or have a single probeset
probeset.intensities.filtered[lapply(probeset.intensities.filtered, function(dataframe){
    length(unique(dataframe$probeset))}) < 1] <- NULL


##### Keep only intensity info
probeset.intensities.only <- lapply(probeset.intensities.filtered, function(dataframe){
    
    dataframe <- dataframe$intensity
})

##### Convert to logical matrix
Affy_HuEx1ST_data <- matrix(unlist(probeset.intensities.only), ncol = length(probeset.intensities.only[[1]]), byrow = TRUE)


rownames(Affy_HuEx1ST_data) <- trFit$groupName[trFit$unitName %in% names(probeset.intensities.filtered)]
colnames(Affy_HuEx1ST_data) <- colnames(trFit)[6:ncol(trFit)]

##### Re-order samples accordingly to target file
Affy_HuEx1ST_data <- Affy_HuEx1ST_data[,rownames(targetFile)]


#===============================================================================
#    Boxplot on normalised data
#===============================================================================

##### Change working directory to the project workspace
setwd(ProjectDir)

Affy_HuEx1ST_datasets <- targetFile$Dataset

Affy_HuEx1ST_datasets_No <- max(as.numeric(factor(Affy_HuEx1ST_datasets)))

Affy_HuEx1ST_datasets.colour <- getDatasetsColours(Affy_HuEx1ST_datasets)


pdf("Affy_HuEx1ST_boxplot.pdf", pointsize = 8 ,width = 0.2*length(targetFile$FileName), height = 6)
par(mar=c(13, 4, 3, 2))
boxplot(data.frame(Affy_HuEx1ST_data),col=Affy_HuEx1ST_datasets.colour[[2]], main="Normalised data", las = 2) # Generates boxplot of normalized log intensity values.
dev.off()

#===============================================================================
#    Remove multi gene to probeset mappings
#===============================================================================

##### Keep only reliable probes based on cigar strings derived from Ensembl microarray probe mapping
Affy_HuEx1ST_probes <- read.table(probesFlatFile,sep="\t",as.is=TRUE,header=TRUE,row.names=1)

Affy_HuEx1ST_data <- Affy_HuEx1ST_data[intersect(rownames(Affy_HuEx1ST_data), rownames(Affy_HuEx1ST_probes)),]


##### Dataset is presumed to contain non-unique gene names, so write the gene names into the first column rather then into the rownames
Affy_HuEx1ST_data <- cbind(Affy_HuEx1ST_probes[intersect(rownames(Affy_HuEx1ST_data), rownames(Affy_HuEx1ST_probes)),1], Affy_HuEx1ST_data)

Affy_HuEx1ST_data  <- multiGene2ProbeFilter(Affy_HuEx1ST_data)

dim(Affy_HuEx1ST_data)

cat("Writing normalised expression data to Affy_HuEx1ST_qn.exp\n")
write.table(prepare2write(Affy_HuEx1ST_data), file="Affy_HuEx1ST_qn.exp",sep="\t", row.names=FALSE)


##### Remove redundant folders
unlink(paste("/scratch/jack/Aroma_R/rawData/", paste(studyIDs, collapse="_"), sep=""), recursive = TRUE)

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
