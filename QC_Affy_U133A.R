################################################################################
#
#   File name: QC_Affy_U133A.R
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
#	Description: Quality check pipeline for Affymetrix Human Genome U133A Array expression array. The script writes in the project workspace a file listing outliers identified with ArrayOutliers function
#
#	Command line use example: R --file=./QC_Affy_U133A.R --args "/scratch/jack/data/PhD/raw/Affy_U133A/GSE32269/target_GSE32269.txt" "/scratch/jack/data/PhD/Transcriptomics_project"
#
#	First arg:      Full path with name of the target file listing data to be analysed. The file is expected to include the following columns: (1) sample name; (2) file name with full path, (3) target and (4) technical replicates indicated by same number
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

#===============================================================================
#    Load libraries
#===============================================================================

library(affy)
library(simpleaffy)
library(arrayQualityMetrics)
library(affyPLM)
library(arrayMvout)
library(hgu133acdf)
library(limma)
library(hgu133a.db)

#===============================================================================
#    Main
#===============================================================================

cat("Affymetrix GeneChip Human Genome U133A Array data quality control.\n")

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
        Dataset2add <- rbind(Dataset2add, matrix(c(studyID,paste(fileInfo,collapse="/"),"raw","Affy_U133A"),1,4,dimnames=list(studyID,colnames(Dataset2add))) )
        write.table( Dataset2add, paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"), row.names = FALSE, col.names = colnames(Dataset2add),  quote = FALSE, sep="\t" )
    }
} else {
    Dataset2add <- matrix(c(studyID,paste(fileInfo,collapse="/"),"raw","Affy_U133A"),1,4,dimnames=list(studyID,c("DatasetName","TargetFile","Type","Platform")))
	write.table( Dataset2add, paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"), row.names = FALSE, col.names = colnames(Dataset2add),  quote = FALSE, sep="\t" )
}

#===============================================================================
#    Load data from affymetrix '.CEL' files
#===============================================================================

pd<-read.AnnotatedDataFrame(targetFile ,header=T, row.name="Name", sep="\t")
dat<-ReadAffy(filenames = paste(pd$FileName,sep=""), sampleNames = sampleNames(pd), phenoData=pd)

##### Change working directory to the folder with QC results
setwd(QCdir)


#===============================================================================
#    Data quality control using arrayQualityMetrics package
#===============================================================================

##### Customise arrayQualityMetrics reports
preparedData = prepdata(expressionset = dat, intgroup = "Target", do.logtransform = TRUE)

QCboxplot <- aqm.boxplot(preparedData, subsample=20000, outlierMethod = "KS")
QCdensity <- aqm.density(preparedData)
QCheatmap <- aqm.heatmap(preparedData)
QCpca <- aqm.pca(preparedData)
QCmaplot <- aqm.maplot(preparedData, subsample=20000, Dthresh=0.15, maxNumArrays=8, nrColumns=4)
QCmeansd <- aqm.meansd(preparedData)

preparedAffy = prepaffy(expressionset = dat, x = preparedData)
##### Affymetrix specific sections
QCrle <- aqm.rle(preparedAffy, outlierMethod = "KS")
QCnuse <- aqm.nuse(preparedAffy, outlierMethod = "upperquartile")

qm = list("Boxplot"=QCboxplot, "Density"=QCdensity, "Heatmap"=QCheatmap, "PCA"=QCpca, "MAplot"=QCmaplot, "MeanSD"=QCmeansd, "RLE"=QCrle, "NUSE"=QCnuse)

aqm.writereport(modules = qm, reporttitle = paste("arrayQualityMetrics report for", studyID, sep=" "), outdir = QCdir, arrayTable = pData(dat))

#arrayQuality <- arrayQualityMetrics(expressionset = dat, outdir = QCdir, force = TRUE, do.logtransform = TRUE, intgroup = "Target")

##### RNA degradation
deg <- AffyRNAdeg(dat, log.it=FALSE)
pdf("deg.pdf")
plotAffyRNAdeg(deg)
dev.off()

##### Data QC with simpleaffy
qc.data <- qc(dat)

pdf("qc.pdf", pointsize = 8 ,width = 6, height = 0.25*length(pd$FileName))
plot(qc.data)
plot(qc.data, usemid=T)
dev.off()

##### multivariate outlier detection
ii = ArrayOutliers(dat, alpha = 0.001, qcOut = qc.data, pc2use=1:3)

write.table(prepare2write(ii[[1]]),"outliers.txt",sep="\t", row.names=FALSE)

##### Report detected outliers in the project workspace
samples2exclude = data.frame(c(studyID, dataDir, paste(ii[[1]]$samp, collapse = ',')))

write.table( t(samples2exclude), paste(ProjectDir,"/outliers_", studyID, ".txt", sep = ""), row.names = FALSE, col.names = c("DatasetName", "DataDir", "Samples2exclude"),  quote = FALSE, sep="\t" )

pdf("outliers.pdf", pointsize = 6)
plot(ii, choices = c(1, 3))
dev.off()


#===============================================================================
#    Normalisation with GC-RMA
#===============================================================================

##### Remove samples identified as outliers
#pd=pd[setdiff(sampleNames(pd),as.vector(ii[[1]]$samp))]
#dat<-ReadAffy(filenames = paste(pd$FileName,sep=""), sampleNames = sampleNames(pd), phenoData=pd)

gcrma=NULL
gcrma = gcrma(dat)

##### signal intensity before and after GC-RMA normalisation
pdf("density.gcrma.pdf")
hist(dat,lty=1)
hist(gcrma,lty=1)
dev.off()

##### Write expression file
write.exprs(gcrma, file="gcrma.exp",sep="\t")

##### Box-and-whisker plot
pdf("box.gcrma.pdf", pointsize = 8 ,width = 0.2*length(pd$FileName), height = 6)
par(mar=c(13, 4, 3, 2))
boxplot(dat,col="red", las = 2) # Generates boxplot of non-normalized intensity values.
boxplot(data.frame(exprs(gcrma)),col="blue", main="Normalised GC-RMA data", las = 2) # Generates boxplot of normalized log intensity values.
dev.off()

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
