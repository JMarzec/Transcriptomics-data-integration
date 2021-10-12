################################################################################
#
#   File name: QC_Affy_HuEx1ST.R
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
#	Description: Quality check pipeline for Affymetrix Human Exon 1.0 ST Array expression array. The script writes in the project workspace a file listing outliers identified with ArrayOutliers function
#
#	Command line use example: R --file=./QC_Affy_HuEx1ST.R --args "/scratch/jack/data/PhD/raw/Affy_HuEx1ST/GSE30521/target_GSE30521.txt" "/scratch/jack/data/PhD/Transcriptomics_project"
#
#	First arg:      Full path with name of the target file listing data to be analysed. The file is expected to include the following columns: (1) sample name; (2) file name with full path, (3) target and (4) technical replicates indicated by same number
#	Second arg:		Project workspace. This is the directory to which cross-platform analysis results will be written
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

#===============================================================================
#    Load libraries
#===============================================================================

library(aroma.affymetrix)

#===============================================================================
#    Main
#===============================================================================

cat("Affymetrix GeneChip Human Exon 1.0 ST Array data quality control.\n")

args <- commandArgs()

fileInfo = unlist(strsplit(args[4], split='/', fixed=TRUE))
targetFile = fileInfo[length(fileInfo)]
dataDir = paste(fileInfo[-length(fileInfo)], collapse = '/')
studyID = fileInfo[length(fileInfo)-1]
QCdir = paste(dataDir, paste("QC", studyID, sep = '_'), sep = '/')
ProjectDir = args[5]

##### Change working directory to the directory with the data
setwd(dataDir)

cat( paste("Processing dataset", studyID, "in", dataDir, ".\n", sep = " "))

##### Set/create a directory for files to be generate
if (file.exists(QCdir)){
	cat( "Folder with QC files already exists.\n" )
} else {
	dir.create(QCdir);
}

##### Set/create a directory for the project
if (file.exists(ProjectDir)){
	cat( paste("The project folder already exists. Writing files to the", ProjectDir, "\n", sep = " ") )
} else {
	dir.create(ProjectDir);
    cat( paste("Writing files to the", ProjectDir, "\n", sep = " ") )
}

##### Add dataset to the lists of datasets to be integrated or create such file if it does not exist yet
if (file.exists(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"))){
    
    Dataset2add <- read.table(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE)[,c(1:4)]
    
    if ( studyID %in% Dataset2add[,"DatasetName"] ) {
        
        cat( paste("Dataset", studyID, "is already in the datasets list.\n", sep = " ") )
    } else {
        Dataset2add <- rbind(Dataset2add, matrix(c(studyID,paste(fileInfo,collapse="/"),"raw","Affy_HuEx1ST"),1,4,dimnames=list(studyID,colnames(Dataset2add))) )
        write.table( Dataset2add, paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"), row.names = FALSE, col.names = colnames(Dataset2add),  quote = FALSE, sep="\t" )
    }
} else {
    Dataset2add <- matrix(c(studyID,paste(fileInfo,collapse="/"),"raw","Affy_HuEx1ST"),1,4,dimnames=list(studyID,c("DatasetName","TargetFile","Type","Platform")))
	write.table( Dataset2add, paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"), row.names = FALSE, col.names = colnames(Dataset2add),  quote = FALSE, sep="\t" )
}

#===============================================================================
#    Load data from affymetrix '.CEL' files
#===============================================================================

pd<-read.table(targetFile, header=T, row.name="Name", sep="\t")
name=pd$FileName

##### Set up directory structure for pipeline based on input type
dir.create(paste("/scratch/jack/Aroma_R/rawData/", studyID,"/",sep=""))
dir.create(paste("/scratch/jack/Aroma_R/rawData/", studyID,"/HuEx-1_0-st-v2",sep=""))
		
cat(paste("Reading CEL files from :", dataDir, "...\n"))
		
##### Copy file over from data directory to aroma directory
for (i in 1:length(name)) {

	file.copy(paste(dataDir,"/",name[i],sep=""), paste("/scratch/jack/Aroma_R/rawData/", studyID,"/HuEx-1_0-st-v2/", name[i], sep=""))
	
	##### Change CEL file names to names specified in target file
	file.rename(paste("/scratch/jack/Aroma_R/rawData/", studyID,"/HuEx-1_0-st-v2/",name[i],sep=""), paste("/scratch/jack/Aroma_R/rawData/", studyID,"/HuEx-1_0-st-v2/",rownames(pd)[i], ".CEL", sep=""))
}

setwd("/scratch/jack/Aroma_R")

##### Specify full messages to the R console
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

##### Get annotation files
tag = "coreR3,A20071112,EP"
annot_level = "coreR3"
chipType <- "HuEx-1_0-st-v2"
cdf <- AffymetrixCdfFile$byChipType(chipType,tags=tag)

##### Define CEL set
cs <- AffymetrixCelSet$byName(studyID, cdf=cdf)
setCdf(cs,cdf) 

##### Change working directory to the folder with QC results
setwd(QCdir)


#===============================================================================
#    Background correction and quantile normalisation using aroma package
#===============================================================================

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

write.table(trFit, "trfit.txt", sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)

################################################################################
##### Exon expression fit exon by exon
#plmEx <- ExonRmaPlm(csN,mergeGroups=FALSE)

##### Run the model
#fit(plmEx)

##### Generate a data frame
#cesEx <- getChipEffectSet(plmEx)
#ExFit <- extractDataFrame(cesEx, units=NULL,addNames=TRUE)

#write.table(ExFit, "exfit.txt", sep = "\t", col.names=TRUE,row.names = FALSE, quote = FALSE)
################################################################################

#===============================================================================
#    Quality assessment
#===============================================================================

##### Box plot of overall signal intensity per sample
pdf("box.gcrma.pdf", pointsize = 8 ,width = 0.2*length(pd$FileName), height = 4)
par(mar=c(8,5,2,1)+0.1)
boxplot(log2(trFit[,6:ncol(trFit)]),col="blue", las = 2, main='Normalised GC-RMA data') # Generates boxplot of normalized intensity values.
dev.off()

##### Density plots
pdf("density.gcrma.pdf")
aroma.affymetrix::plotDensity(cs, ylim=c(0, 1), main='Raw data')
aroma.affymetrix::plotDensity(csBC, ylim=c(0, 1), main='Background corrected data')
aroma.affymetrix::plotDensity(csN, ylim=c(0, 1), main='Normalised GC-RMA data')
dev.off()

##### Quality assessment of PLM fit (NUSE and RLE plots)
qamTr <- QualityAssessmentModel(plmTr)

pdf("NUSE.pdf", pointsize = 8 ,width = 0.2*length(pd$FileName), height = 4)
par(mar=c(8,5,2,1)+0.1)
plotNuse(qamTr)
dev.off()

pdf("RLE.pdf", pointsize = 8 ,width = 0.2*length(pd$FileName), height = 4)
par(mar=c(8,5,2,1)+0.1)
plotRle(qamTr)
dev.off()

pdf("RLE.pdf", pointsize = 8 ,width = 0.2*length(pd$FileName), height = 4)
par(mar=c(8,5,2,1)+0.1)
plotRle(qamTr)
dev.off()

##### Create file for reporting outliers
samples2exclude = data.frame(c(studyID, dataDir, ""))

write.table( t(samples2exclude), paste(ProjectDir,"/outliers_", studyID, ".txt", sep = ""), row.names = FALSE, col.names = c("DatasetName", "DataDir", "Samples2exclude"),  quote = FALSE, sep="\t" )

##### Remove redundant folders
unlink(paste(QCdir, "probeData", sep="/"), recursive = TRUE)
unlink(paste(QCdir, "plmData", sep="/"), recursive = TRUE)
unlink(paste("/scratch/jack/Aroma_R/rawData/", studyID, sep=""), recursive = TRUE)

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
