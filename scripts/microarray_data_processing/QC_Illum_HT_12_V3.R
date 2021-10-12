################################################################################
#
#   File name: QC_Illum_HT_12_V3.R
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
#	Description: Quality check pipeline for Illumina Human HT-12 v3 expression beadchip. The input expression file is expected to be in Illumina BeadStudio output format and is required to include the columns with the following headers per sample: (1) AVG_Signal; (2) BEAD_STDEV (if not available the values for all probes may be set to 'NaN' and (3) Detection (refers to detection p-value), as reported in 'lumi' package reference manual ( http://bioconductor.org/packages/2.12/bioc/vignettes/lumi/inst/doc/lumi.pdf ). The script writes in the project workspace a file listing outliers identified with ArrayOutliers function.
#
#	Command line use example: R --file=./QC_Illum_HT_12_V3.R --args "/scratch/jack/data/PhD/raw/Illum_HT_12_V3/GSE32571/target_GSE32571.txt" "/scratch/jack/data/PhD/Transcriptomics_project"
#
#	First arg:      Full path with name of the target file listing data to be analysed. The file is expected to include the following columns: (1) sample name; (2) file name, (3) target and (4) technical replicates indicated by same number
#	Second arg:      Project workspace. This is the directory to which cross-platform analysis results will be written
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

##### Replace signal intensities of '0' with minimal value detected on the array
replace0s <- function (lumiBatch) {
	
	for (i in 1:ncol(exprs(lumiBatch))) {
		
		exprs.ordered <- sort(exprs(lumiBatch)[,i])
		exprs.ordered <- exprs.ordered[ exprs.ordered != 0 ]
				
		exprs(lumiBatch)[exprs(lumiBatch)[,i]==0] <- exprs.ordered[1]
	}	
	return(lumiBatch)
}

#===============================================================================
#    Load libraries
#===============================================================================

library(lumi)
library(arrayQualityMetrics)

#===============================================================================
#    Main
#===============================================================================

cat("Illumina HumanHT-12 v3.0 Expression Beadchip data quality control.\n")

args <- commandArgs()

fileInfo = unlist(strsplit(args[4], split='/', fixed=TRUE))
targetFile = fileInfo[length(fileInfo)]
dataDir = paste(fileInfo[-length(fileInfo)], collapse = '/')
studyID = fileInfo[length(fileInfo)-1]
QCdir = paste(dataDir, paste("QC", studyID, sep = '_'), sep = '/')
ProjectDir = args[5]

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
        Dataset2add <- rbind(Dataset2add, matrix(c(studyID,paste(fileInfo,collapse="/"),"raw","Illum_HT_12_V3"),1,4,dimnames=list(studyID,colnames(Dataset2add))) )
        write.table( Dataset2add, paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"), row.names = FALSE, col.names = colnames(Dataset2add),  quote = FALSE, sep="\t" )
    }
} else {
    Dataset2add <- matrix(c(studyID,paste(fileInfo,collapse="/"),"raw","Illum_HT_12_V3"),1,4,dimnames=list(studyID,c("DatasetName","TargetFile","Type","Platform")))
	write.table( Dataset2add, paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"), row.names = FALSE, col.names = colnames(Dataset2add),  quote = FALSE, sep="\t" )
}

##### Change working directory to the folder with QC results
setwd(QCdir)

#===============================================================================
#    Load expression data
#===============================================================================

target <- as.data.frame(read.table(paste(dataDir,targetFile,sep="/"),sep="\t",as.is=TRUE,header=TRUE))
fileName = paste(dataDir, target[1,"FileName"], sep="/")

target <- as.data.frame(read.table(paste(dataDir,targetFile,sep="/"),sep="\t",as.is=TRUE,header=TRUE)[,"Target"])

lumiData = lumiR( fileName, sep = "\t", detectionTh = 0.01, na.rm = TRUE, lib = NULL)

rownames(target) <- read.table(paste(dataDir,targetFile,sep="/"),sep="\t",as.is=TRUE,header=TRUE)[,"Name"]
colnames(target) <- "Target"
sampleNames(lumiData) <- rownames(target)

data <- new("ExpressionSet", exprs = exprs(lumiData), phenoData = new("AnnotatedDataFrame", data=target) )

#===============================================================================
#    Data quality control using arrayQualityMetrics package
#===============================================================================

##### Basic data quality using lumi package
QCsummary <- lumiData@QC$sampleSummary

##### Write QC summary
write.table(prepare2write(QCsummary), file="QC_summary.txt",sep="\t", row.names=FALSE)

##### Customise arrayQualityMetrics reports
preparedData = prepdata(expressionset = data, intgroup = "Target", do.logtransform = TRUE)

QCboxplot <- aqm.boxplot(preparedData, subsample=20000, outlierMethod = "KS")
QCdensity <- aqm.density(preparedData)
QCheatmap <- aqm.heatmap(preparedData)
QCpca <- aqm.pca(preparedData)
QCmaplot <- aqm.maplot(preparedData, subsample=20000, Dthresh=0.15, maxNumArrays=8, nrColumns=4)
QCmeansd <- aqm.meansd(preparedData)

qm = list("Heatmap"=QCheatmap, "PCA"=QCpca, "Boxplot"=QCboxplot, "Density"=QCdensity, "MeanSD"=QCmeansd, "MAplot"=QCmaplot)

aqm.writereport(modules = qm, reporttitle = paste("arrayQualityMetrics report for", studyID, sep=" "), outdir = QCdir, arrayTable = pData(data))


##### multivariate outlier detection
library(arrayMvout)
ii = ArrayOutliers(lumiData, alpha = 0.001, pc2use=1:3)

write.table(prepare2write(ii[[1]]),"outliers.txt",sep="\t", row.names=FALSE)

##### Report detected outliers in the project workspace
samples2exclude = data.frame(c(studyID, dataDir, paste(rownames(ii[[1]]), collapse = ',')))

write.table( t(samples2exclude), paste(ProjectDir,"/outliers_", studyID, ".txt", sep = ""), row.names = FALSE, col.names = c("DatasetName", "DataDir", "Samples2exclude"),  quote = FALSE, sep="\t" )

pdf("outliers.pdf", pointsize = 6)
plot(ii, choices = c(1, 3))
dev.off()
 
#===============================================================================
#    Normalisation with Robust Spline Normalization (RSN)
#===============================================================================

##### Remove samples identified as outliers
#data.exprs=data.exprs[setdiff(sampleNames(lumiData),as.vector(rownames(ii[[1]])))]

##### RSN method combines features from loess and quantile normalization methods. It was recommended by Schmid et al ( http://www.biomedcentral.com/1471-2164/11/349 ) as it  combines the positive effects of quantile normalization, i.e. preservation of the rank order, and spline interpolation, i.e. continuous mapping of the values, and at the same time circumvents their drawbacks, i.e. discontinuous mapping of intensity values and no rank preservation, respectively 

##### Replace signal intensities of '0', which introduce errors in RSN, with minimal value detected on the array
data  <- replace0s(lumiData)
data.N <- lumiN(data, method='rsn')

#===============================================================================
#    Quality control with lumi package on normalised data
#===============================================================================

data.N.Q <- lumiQ(data.N)

QCsummary <- data.N.Q@QC$sampleSummary

##### Write QC summary
write.table(prepare2write(QCsummary), file="QC_summary.rsn.txt",sep="\t", row.names=FALSE)

##### signal intensity before and after RSN normalisation
pdf("density.rsn.pdf")
plot(lumiData, what='density')
plot(data.N.Q, what='density')
dev.off()

##### Sample relations using hierarchical clustering before and after RSN normalisation
pdf("Sample_relations_cluster.rsn.pdf", pointsize = 8 ,width = 0.2*length(sampleNames(data)), height = 6)
plot(lumiData, what='sampleRelation', hang=-1)
plot(data.N.Q, what='sampleRelation', hang=-1)
dev.off()

##### Sample relations using multidimensional scaling before and after RSN normalisation
pdf("Sample_relations_MDS.rsn.pdf")
plot(lumiData, what='sampleRelation', method='mds')
plot(data.N.Q, what='sampleRelation', method='mds')
dev.off()

##### Box-and-whisker plot
pdf("box.rsn.pdf", pointsize = 8 ,width = 0.2*length(sampleNames(data)), height = 6)
par(mar=c(13, 4, 3, 2))
boxplot(log2(exprs(lumiData)),col="red", las = 2) # Generates boxplot of non-normalized intensity values.
boxplot(data.frame(log2(exprs(data.N.Q))),col="blue", main="Normalised RSN data", las = 2) # Generates boxplot of normalized log intensity values.
dev.off()

##### Write expression file
write.table(prepare2write(log2(exprs(data.N.Q))), file="rsn.exp",sep="\t", row.names=FALSE)

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
