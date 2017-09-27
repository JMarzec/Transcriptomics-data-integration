################################################################################
#
#   File name: Study_effect.R
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
#   Description: Pipeline for performing study effect analysis for defined datasets. It investigates the distribution of all pairwise correlations (1) between the samples in all datasets; (2) between the sample subgroups from different studies within the same biological group and (3) between the different biological groups within the same study. It also performs differential expression analysis on selected datasets. The script reads in a priori generated files with normalised expression data
#
#   Command line use example: R --file=./Study_effect.R --args "/scratch/jack/data/PhD/Transcriptomics_project" "Affy_U133Plus2" "GSE17951,GSE37199" "nonparam"
#
#   First arg:      Project workspace. This is the directory to which cross-platform analysis results will be written
#   Second arg:     Platforms to be considered for the analysis. Currently suported platforms are [RNAseq], [Affy_HuEx1ST], [Affy_U133Plus2], [Affy_U133A2], [Affy_U133A], [Affy_U133B], [Affy_U95Av2], [Affy_U95B], [Affy_U95C], [Illum_HT_12_V3]. Type 'all' for all supported platforms to be considered
#   Third arg:      Datasets, processed with seleceted platforms in third arg, to be considered for the analysis. Each datasets ID is expected to be separated by comma. Type 'all' for all listed datasets to be analysed
#   Forth arg (OPTIONAL): If "nonparam", use nonparametric adjustments in ComBat. AS default, parametric adjustments are used

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
    dataBound <- read.table(datasetList[[1]],sep="\t",as.is=TRUE,header=TRUE,row.names=1,check.names=FALSE)
    
	for (i in 2:length(datasetList)) {
		
        cat(paste("Loading expression data", datasetList[[i]], "...\n"), sep="")
        data2bind <- read.table(datasetList[[i]],sep="\t",as.is=TRUE,header=TRUE,row.names=1,check.names=FALSE)
        
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

##### Assign colours to analysed groups
getTargetsColours <- function(targets) {
    
    ##### Predefined selection of colours for groups
    targets.colours <- c("red","blue","green","darkgoldenrod","darkred","deepskyblue")
    
    f.targets <- factor(targets)
    vec.targets <- targets.colours[1:length(levels(f.targets))]
    targets.colour <- rep(0,length(f.targets))
    for(i in 1:length(f.targets))
    targets.colour[i] <- vec.targets[ f.targets[i]==levels(f.targets)]
    
    return( list(vec.targets, targets.colour) )
}

#===============================================================================
#    Load libraries
#===============================================================================

library(sva)
library(limma)
library(gplots)
library(geneplotter)
source("a2R_code.R")

#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

ProjectDir = args[4]
platform = args[5]
datasets2corr = args[6]
datasets2corr = gsub("\\s","", datasets2corr)
datasets2corr =  sort(unlist(strsplit(datasets2corr, split=',', fixed=TRUE)))
parPrior = args[7]

if (!is.na(parPrior) && parPrior == "nonparam") {
	parPrior = FALSE
} else {
	parPrior = TRUE
}

##### Read file with datasets information
DatasetInput=read.table(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE,row.names=1)

if ( platform != "all" && platform != "ALL" ) {
    
    DatasetInput=DatasetInput[DatasetInput[,"Platform"] %in% platform,]
} else {
    platform=unique(DatasetInput[,"Platform"])
}

if ( datasets2corr != "all" && datasets2corr != "ALL" ) {
    
    if ( all(datasets2corr %in% rownames(DatasetInput)) ) {
        DatasetInput=DatasetInput[rownames(DatasetInput) %in% datasets2corr,]
        
    } else {
        cat(paste("\nSome of selected datasets (", paste(datasets2corr, collapse = ", "), ") were not generated by selected platforms (", paste(platform, collapse = ", "), ")!\n\n", sep=" "))
        q()
    }
} else {
     datasets2corr <- sort(rownames(DatasetInput))
}

##### Change working directory to the project workspace
setwd(ProjectDir)

##### Report used parameters to a file
write(args, file = paste("Study_effect", paste(datasets2corr, collapse="_"), "parameters.txt", sep = "_"), append = FALSE, sep="\t")


#===============================================================================
#    Load normalised data
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

datasetList <- NULL

for ( i in 1:length(levels(factor(DatasetInput[,3]))) ) {
	
	if ( unique(DatasetInput[,"Platform"])[i] == "Affy_HuEx1ST" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"Affy_HuEx1ST_qn.exp",sep="/")
		filterThreshold <- 10000
		
	} else if ( unique(DatasetInput[,"Platform"])[i] == "Affy_U133Plus2" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"Affy_U133Plus2_gcrma.exp",sep="/")
		filterThreshold <- 10000
		
	} else if ( unique(DatasetInput[,"Platform"])[i] == "Affy_U133A2" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"Affy_U133A2_gcrma.exp",sep="/")
		filterThreshold <- 5000
		
	} else if ( unique(DatasetInput[,"Platform"])[i] == "Affy_U133A" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"Affy_U133A_gcrma.exp",sep="/")
		filterThreshold <- 5000
		
	} else if ( unique(DatasetInput[,"Platform"])[i] == "Affy_U133B" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"Affy_U133B_gcrma.exp",sep="/")
		filterThreshold <- 5000
		
	} else if ( unique(DatasetInput[,"Platform"])[i] == "Affy_U95Av2" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"Affy_U95Av2_gcrma.exp",sep="/")
		filterThreshold <- 3000
		
	} else if ( unique(DatasetInput[,"Platform"])[i] == "Affy_U95B" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"Affy_U95B_gcrma.exp",sep="/")
		filterThreshold <- 3000
		
	} else if ( unique(DatasetInput[,"Platform"])[i] == "Affy_U95C" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"Affy_U95C_gcrma.exp",sep="/")
		filterThreshold <- 3000
		
	} else if ( unique(DatasetInput[,"Platform"])[i] == "Illum_HT_12_V3" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"Illum_HT_12_V3_rsn.exp",sep="/")
        filterThreshold <- 10000
        
	} else if ( unique(DatasetInput[,"Platform"])[i] == "RNAseq" ) {
        
		datasetList[[i]] <- paste(ProjectDir,"RNAseq_cqn.exp",sep="/")
        filterThreshold <- 30000
        
	} else {
		cat(paste("Platform", unique(DatasetInput[,"Platform"]) ,"is not yet supported!\n\n", sep=" "))
        q()
	}
}


if ( length(datasetList) > 1 ) {
    
    cat("Binding datasets...\n")
    data <- bindDatasets(datasetList)
    
} else if ( length(datasetList) == 1 ) {
    
     cat(paste("Loading expression data", datasetList[[1]], "...\n"), sep="")
     data <- read.table(datasetList[[1]],sep="\t",as.is=TRUE,header=TRUE,row.names=1,check.names=FALSE)
     
} else {
	cat("No datasets provided!\n\n")
    q()
}

dim(data)

##### Keep only samples present in normalised data files
targetFile <- targetFile[intersect(colnames(data), rownames(targetFile)),]


##### Select data to be analysed
if ( datasets2corr != "all" && datasets2corr != "ALL" ) {
    data <- data[,rownames(targetFile)]
    targetFile <- targetFile[setdiff(colnames(data),as.vector(list(targetFile$Name))),]
    dim(data)
}


datasets <- targetFile$Dataset

datasets_No <- max(as.numeric(factor(datasets)))

datasets.colour <- getDatasetsColours(datasets)

targets <- targetFile$Target

targets.colour <- getTargetsColours(targets)


#===============================================================================
#     View the trend of the correlation between the genes average expression level among studies
#===============================================================================

if (datasets_No > 1 ) {
    rowMeans <- matrix(nrow=nrow(data), ncol=datasets_No)

    for (i in 1:datasets_No) {
    
        datasetSamples <- rownames( targetFile[targetFile[,4]==unique(targetFile[,4])[i],] )

        for (j in 1:nrow(data)) {
            rowMeans[j,i] <- apply(data[j,intersect(colnames(data),as.vector(datasetSamples))], 1, mean)
    
        }
    }

    rownames(rowMeans) <- rownames(data)
    colnames(rowMeans) <- unique(targetFile[,4])

    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_pairwise_plot.pdf", sep=""),width = 3*datasets_No, height = 3*datasets_No)
    pairs(rowMeans, gap = 0, pch = ".")
    dev.off()


    #===============================================================================
    #     Pairwise correlations between samples
    #===============================================================================

    ##### Get datasets sizes
    datasetDim <- NULL

    for (i in 1:datasets_No) {
    
        datasetDim[i] <- nrow( targetFile[targetFile[,4]==unique(targetFile[,4])[i],] )
    }
    datasetDim <- t(data.frame(datasetDim))
    colnames(datasetDim) <- unique(targetFile[,4])


    ##### Calculate Pearson correlation coefficients between sample subgroups from different studies within the same biological group
    group_corr <- NULL
    group_corr.list <- NULL

    datasetStart = 0
    datasetEnd = 0
    ##### Go through selected datasets
    for (i in seq(1,datasets_No-1,by=1)) {
    
        datasetStart = datasetEnd + 1
        datasetEnd = datasetStart + datasetDim[colnames(datasetDim)==unique(targetFile[,4])[i]] - 1
    
        cat(paste("dataset:",targetFile[datasetStart,4],"start:",datasetStart,"end",datasetEnd,"\n",sep=" "))
    
        ##### Go through each sample in a dataset and...
        for (j in seq(datasetStart,datasetEnd,by=1)) { #cat(paste("\nsample:",colnames(data)[j],"\n",sep=" "))

            ##### ...correlate with remaining samples from other datasets...
            for(k in seq(datasetEnd+1,ncol(data),by=1)) { #cat(paste("agaist:",colnames(data)[k],"\n",sep=" "))
  
                ##### ... consider samples of the same subgroup only
                if(targets[j] == targets[k]) { #cat(paste(colnames(data)[k],"group match:",targets[j],targets[k],"\n",sep=" "))
                
                    group_corr <- c(group_corr, cor(x=data[,j],y=data[,k], use="all.obs", method="pearson"))
                
                }
            }
        }
        group_corr.list[[i]] <- group_corr
    }


    ##### Calculate Pearson correlation coefficients between the different biological groups within the same study
    study_corr <- NULL
    study_corr.list <- NULL

    datasetStart = 0
    datasetEnd = 0
    ##### Go through selected datasets
    for (i in seq(1,datasets_No,by=1)) {
    
        datasetStart = datasetEnd + 1
        datasetEnd = datasetStart + datasetDim[colnames(datasetDim)==unique(targetFile[,4])[i]] - 1
    
        cat(paste("dataset:",targetFile[datasetStart,4],"start:",datasetStart,"end",datasetEnd,"\n",sep=" "))
    
        ##### Go through each sample in a dataset and...
        for (j in seq(datasetStart,datasetEnd-1,by=1)) {  #cat(paste("\nsample:",colnames(data)[j],"\n",sep=" "))

            ##### ...correlate with remaining samples from same dataset...
            for(k in seq(j+1,datasetEnd,by=1)) {
            
                #cat(paste("agaist:",colnames(data)[k],"\n",sep=" "))
            
                ##### ... consider samples of the different subgroup only
                if(targets[j] != targets[k]) { #cat(paste(colnames(data)[k],"group not match:",targets[j],targets[k],"\n",sep=" "))

                    study_corr <- c(study_corr, cor(x=data[,j],y=data[,k], use="all.obs", method="pearson"))
                }
            }
        }
        study_corr.list[[i]] <- study_corr
    }


    ##### Calculate Pearson correlation coefficients between the samples in all datasets
    all_corr <- NULL

    ##### Go through each sample in all datasets and...
    for (i in 1:ncol(data)) { #cat(paste("\nsample:",colnames(data)[i],"\n",sep=" "))
    
        j=i
        ##### ...correlate with remaining samples
        while(j < ncol(data)) { #cat(paste("agaist:",colnames(data)[j+1],"\n",sep=" "))
        
            j = j + 1
            all_corr <- c(all_corr, cor(x=data[,i],y=data[,j], use="all.obs", method="pearson"))
        }
    }

    ##### Draw density plots of Pearson correlation coefficients
    if ( is.null(group_corr.list) ) {
    
        cat("The bilogical effect cannot not be estimated!\n\n")
    
    } else if ( is.null(study_corr.list) ) {
    
        cat("The study effect cannot be estimated!\n\n")
    } else {

        density1 <- density(unlist(group_corr.list))

        density2 <- density(unlist(study_corr.list))

        density3 <- density(unlist(all_corr))

        minX <- min(unlist(group_corr.list),unlist(study_corr.list),all_corr)-0.2
        maxY <- max(density1$y,density2$y,density3$y)

        pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_density.pdf", sep="") ,width = 6, height = 6)
        plot(density1, xlim=c(minX,1), ylim=c(0,maxY), main="", xlab="Pearson correlation coefficient", col="darkolivegreen4")
        par(new=TRUE)
        plot(density2, xlim=c(minX,1), ylim=c(0,maxY), main="", xlab="", ylab="", xaxt='n', yaxt='n', col="red")
        par(new=TRUE)
        plot(density3, xlim=c(minX,1), ylim=c(0,maxY), main="", xlab="", ylab="", xaxt='n', yaxt='n',, col="black")
        legend(minX,maxY, legend=c(" ALL"," Study"," Group"),col=c("black","red","darkolivegreen4"), lty=1, box.col="black")
        dev.off()
    }
} else {
    
    #===============================================================================
    #     Non-specific filtering
    #===============================================================================
    
    ##### Depending on the platform, keep 10,000, 5,000 or 3,000 of all genes with the highest expression variance across samples
    ##### ... unless the initial data contains less genes (possible for RNA-seq datasets with different number of genes measured), then use 60% of these genes
    if ( nrow(data) < filterThreshold ) {
        
        filterThreshold <- 0.6 * nrow(data)
    }
    
    rsd<-apply(data,1,sd)
    sel<-order(rsd, decreasing=TRUE)[1:filterThreshold]
    data <-data[sel,]
}

#===============================================================================
#     Unsupervised clustering
#===============================================================================

d.usa <- dist(t(data), "euc")
h.usa <- hclust(d.usa, method="ward")

#####  Adopt cluster plots width
if ( ncol(data) < 20 ) {
    
	clust.width = 0.3*ncol(data)
    
} else {
	clust.width = 0.1*ncol(data)
}


#####  Set.seed returns NULL, invisibly
set.seed(1)

if (datasets_No > 1 ) {
    
    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_cluster_datasets.pdf", sep=""),width = clust.width, height = 6)
    h.usa$labels = colnames(data)
    par(mar=c(2,2,2,6))
    A2Rplot(h.usa,
    k=length(levels(factor(datasets))), # k = changes the detail/number of subgroups shown.
    fact.sup=datasets, box=FALSE,
    show.labels=TRUE, col.up = "black", main=" ")
    dev.off()

    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_cluster_targets.pdf", sep=""),width = clust.width, height = 6)
    h.usa$labels = colnames(data)
    par(mar=c(2,2,2,6))
    A2Rplot(h.usa,
    k=length(levels(factor(datasets))), # k = changes the detail/number of subgroups shown.
    fact.sup=targets, box=FALSE,
    show.labels=TRUE, col.up = "black", main=" ")
    dev.off()
    
} else {
    
    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_cluster_targets.pdf", sep=""),width = clust.width, height = 6)
    h.usa$labels = colnames(data)
    par(mar=c(2,2,2,6))
    A2Rplot(h.usa,
    k=length(levels(factor(targets))), # k = changes the detail/number of subgroups shown.
    fact.sup=targets, box=FALSE,
    show.labels=TRUE, col.up = "black", main=" ")
    dev.off()
}

pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"),"_dendrogram.pdf", sep = ""), width = 0.2*ncol(data)+2, height = 6, pointsize = 12)
par(mar=c(2,5,2,0))
plot(h.usa, xlab="", labels=paste(colnames(data), targets, sep="       "), hang = -1, main="")
dev.off()


#===============================================================================
#     K-means clustering
#===============================================================================

fit_2D <- cmdscale(d.usa,eig=TRUE, k=2) # k is the number of dim

#####  Plot solution
x <- fit_2D$points[,1]
y <- fit_2D$points[,2]

fit_2D.points = cbind(x, y)
colnames(fit_2D.points) <- c("Coordinate 1", "Coordinate 2")
cl <- kmeans(fit_2D.points, 1)

#####  Set pch symbols so that each represents single dataset on the plots
pchs <- rep(list(c(16,1),c(15,0),c(17,2),c(16,10),c(15,7),c(16,13),c(16,1)),10)

pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_K_means.pdf", sep=""), width = 7.5, height = 7.5)
plot(fit_2D.points, col = cl$cluster, type="n", ylim=c( min(y),max(y) + (abs(min(y))+abs(max(y)))/4 ))
#####  Use different shape for each dataset
for (i in 1:datasets_No) {
    points(x[as.numeric(factor(datasets))==i], y[as.numeric(factor(datasets))==i], cex=1, pch = pchs[[i]][1], col=targets.colour[[2]][as.numeric(factor(datasets))==i])
    points(x[as.numeric(factor(datasets))==i], y[as.numeric(factor(datasets))==i], cex=1, pch = pchs[[i]][2], col="black")
}
#####  Add the legend
legend("topleft", legend=levels(factor(targets)),pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()


#===============================================================================
#     Principal components analysis
#===============================================================================

#####  Keep only probes with variance > 0 across all samples
rsd<-apply(data,1,sd)
data <- data[rsd>0,]

#####  Perform principal components analysis after scaling the data
data_pca <- prcomp(t(data), scale=TRUE)


pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_PCA_screeplot.pdf", sep=""), width = 14, height = 6)
screeplot(data_pca, npcs = length(data_pca$sdev), main="The variances captured by principal components")
dev.off()

#####  Get variance importance for all principal components
importance_pca <- summary(data_pca)$importance[2,]
importance_pca <- paste(round(100*importance_pca, 2), "%", sep="")


pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_PCA.pdf", sep=""))
plot(data_pca$x[,1], data_pca$x[,2], type="n", xlab = paste("PC1 (",importance_pca[1],")",sep=""), ylab = paste("PC2 (",importance_pca[2],")",sep=""), ylim=c( min(data_pca$x[,2]),max(data_pca$x[,2]) + (abs(min(data_pca$x[,2]))+abs(max(data_pca$x[,2])))/4 ))

#####  Use different shape for each dataset
for (i in 1:datasets_No) {
    points(data_pca$x[,1][as.numeric(factor(datasets))==i],data_pca$x[,2][as.numeric(factor(datasets))==i], pch=pchs[[i]][1], col=targets.colour[[2]][as.numeric(factor(datasets))==i])
    points(data_pca$x[,1][as.numeric(factor(datasets))==i],data_pca$x[,2][as.numeric(factor(datasets))==i], pch = pchs[[i]][2], col="black")
}
#####  Adding the legend
legend("topleft", legend=levels(factor(targets)),pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()


pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_PCA_labeled.pdf", sep=""))
plot(data_pca$x[,1], data_pca$x[,2], type="n", xlab = paste("PC1 (",importance_pca[1],")",sep=""), ylab = paste("PC2 (",importance_pca[2],")",sep=""), ylim=c( min(data_pca$x[,2]),max(data_pca$x[,2]) + (abs(min(data_pca$x[,2]))+abs(max(data_pca$x[,2])))/4 ))

#####  Use different shape for each dataset
for (i in 1:datasets_No) {
    points(data_pca$x[,1][as.numeric(factor(datasets))==i],data_pca$x[,2][as.numeric(factor(datasets))==i], pch=pchs[[i]][1], col=targets.colour[[2]][as.numeric(factor(datasets))==i])
    points(data_pca$x[,1][as.numeric(factor(datasets))==i],data_pca$x[,2][as.numeric(factor(datasets))==i], pch = pchs[[i]][2], col="black")
}
#####  Adding sample names above symbols
text(data_pca$x[,1], data_pca$x[,2]+(abs(min(data_pca$x[,2]))+abs(max(data_pca$x[,2])))/50, labels= colnames(data), cex=0.5)
#####  Adding the legend
legend("topleft", legend=levels(factor(targets)),pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()


pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_PCA_pairs.pdf", sep=""))
pairs(data_pca$x[,1:4], pch=20, col=targets.colour[[2]]) # Plots all combinations between the first 4 pc
dev.off()


#####  Smotther scatter plot - shows the density of the data points
pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_smoothScatter.pdf", sep=""))
smoothScatter(data_pca$x) #  A smotther scatter plot - shows the density of the data points.
dev.off()


#===============================================================================
#     Remove batch effects with ComBat
#===============================================================================

if (datasets_No > 1 ) {
    
    #####  Obtain known batch effects
    batch = as.vector(datasets)

    #####  Create model matrix for outcome of interest and other covariates besides batch
    f.model <- factor(targets, levels = levels(factor(targets)))

    if ( length(levels(f.model))>1 ) {
        mod <- model.matrix(~f.model)
    } else {
        mod = NULL
    }

    #####  Create an environment and assign target categories (covariates) to samples from each batch
    batchenv<-new.env()

    for(i in 1:length(batch)) {
        batchenv[[ batch[i] ]] <- c(batchenv[[ batch[i] ]],f.model[i])
    }

    #####  Check whether any of known batches are confounded with any of specified covariates ( problem affecting design matrix discussed at glokbase forum: http://grokbase.com/t/r/bioconductor/138km86v3q/bioc-fwd-combat-error-in-solve-default-t-design-design-lapack-routine-dgesv-system-is-exactly-singular )
    for(i in 1:length(ls(batchenv))) {
    
        #####  If any of the covariates are confounded with any batch then stop the analysis
        if ( length(unique(batchenv[[ls(batchenv)[i]]])) == 1 ) {
        
            cat("Data combinationed stopped due to presence of covariates which seem to be confounded with one of the provided batched" )
        
            #################### THIS PART NEEDS MORE TESTING ####################
            #####  If any of the covariates are confounded with any batch then consider only "Normal" covariate in the design matrix
            #if ( length(intersect(colnames(mod),c("f.modelNormal")))>0 ) {
            #
            #    mod=mod[,intersect(colnames(mod),c("f.modelNormal"))]
            #
            #####  If any of the covariates are confounded with any batch and there in no "Normal" covariate then igonre the design matrix
            #} else if ( length(intersect(colnames(mod),c("f.modelNormal")))==0 ) {
            #
            #    mod=NULL
            #}
            ######################################################################
        }
    }

    cat("Adjusting data for batch effects using ComBat...\n")
    if ( isTRUE(parPrior) ) {
        ##### Generate prior plots to assess if parametric or nonparametric adjustments should be used in ComBat. If the red (parametric estimate of the empirical batch effect density) and black (kernel estimate of the empirical batch effect density) lines don't match up well, use the nonparametric ajustments
        pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_PriorPlot.pdf", sep=""), pointsize = 12, bg = "transparent")
        data_combat = ComBat(dat=data, batch=batch, mod=mod, par.prior=parPrior, prior.plots=TRUE)
        dev.off()
    } else {
        data_combat = ComBat(dat=data, batch=batch, mod=mod, par.prior=parPrior, prior.plots=FALSE)
    }

    #===============================================================================
    #     Boxplot of combined datasets
    #===============================================================================

    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_boxplot.pdf", sep=""), pointsize = 8 ,width = 0.2*ncol(data), height = 6)
    par(mar=c(13, 4, 3, 2))
    boxplot(data.frame(data),col=datasets.colour[[2]], las = 2) # Generates boxplot of all datasets.
    boxplot(data.frame(data_combat),col=datasets.colour[[2]], main="Batch effect adjusted", las = 2) # Generates boxplot of batch effect corrected data.
    dev.off()

    #===============================================================================
    #     View the trend of the correlation between the genes average expression level among combined studies
    #===============================================================================

    rowMeans <- matrix(nrow=nrow(data_combat), ncol=datasets_No)

    for (i in 1:datasets_No) {
    
        datasetSamples <- rownames( targetFile[targetFile[,4]==unique(targetFile[,4])[i],] )
    
        for (j in 1:nrow(data_combat)) {
            rowMeans[j,i] <- apply(data_combat[j,intersect(colnames(data_combat),as.vector(datasetSamples))], 1, mean)
        
        }
    }

    rownames(rowMeans) <- rownames(data_combat)
    colnames(rowMeans) <- unique(targetFile[,4])

    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_pairwise_plot_ComBat.pdf", sep=""),width = 3*datasets_No, height = 3*datasets_No)
    pairs(rowMeans, gap = 0, pch = ".")
    dev.off()


    #===============================================================================
    #     Non-specific filtering
    #===============================================================================

    ##### Depending on the platform, keep 10,000, 5,000 or 3,000 of all genes with the highest expression variance across samples
    ##### ... unless the data contains less genes (possible for RNA-seq datasets with different number of genes measured), then use 60% of these genes
    if ( nrow(data_combat) < filterThreshold ) {
        
        filterThreshold <- 0.6 * nrow(data_combat)
    }
    
    rsd<-apply(data_combat,1,sd)
    sel<-order(rsd, decreasing=TRUE)[1:filterThreshold]
    data_combat <-data_combat[sel,]


    #===============================================================================
    #     Pairwise correlations between samples
    #===============================================================================

    ##### Calculate Pearson correlation coefficients between sample subgroups from different studies within the same biological group
    group_corr <- NULL
    group_corr.list <- NULL

    datasetStart = 0
    datasetEnd = 0
    ##### Go through selected datasets
    for (i in seq(1,datasets_No-1,by=1)) {
    
        datasetStart = datasetEnd + 1
        datasetEnd = datasetStart + datasetDim[colnames(datasetDim)==unique(targetFile[,4])[i]] - 1
    
        cat(paste("dataset:",targetFile[datasetStart,4],"start:",datasetStart,"end",datasetEnd,"\n",sep=" "))
    
        ##### Go through each sample in a dataset and...
        for (j in seq(datasetStart,datasetEnd,by=1)) { #cat(paste("\nsample:",colnames(data_combat)[j],"\n",sep=" "))
        
            ##### ...correlate with remaining samples from other datasets...
            for(k in seq(datasetEnd+1,ncol(data_combat),by=1)) { #cat(paste("agaist:",colnames(data_combat)[k],"\n",sep=" "))
            
                ##### ... consider samples of the same subgroup only
                if(targets[j] == targets[k]) { #cat(paste(colnames(data_combat)[k],"group match:",targets[j],targets[k],"\n",sep=" "))
                
                    group_corr <- c(group_corr, cor(x=data_combat[,j],y=data_combat[,k], use="all.obs", method="pearson"))
                }
            }
        }
        group_corr.list[[i]] <- group_corr
    }


    ##### Calculate Pearson correlation coefficients between the different biological groups within the same study
    study_corr <- NULL
    study_corr.list <- NULL

    datasetStart = 0
    datasetEnd = 0
    ##### Go through selected datasets
    for (i in seq(1,datasets_No,by=1)) {
    
        datasetStart = datasetEnd + 1
        datasetEnd = datasetStart + datasetDim[colnames(datasetDim)==unique(targetFile[,4])[i]] - 1
    
        cat(paste("dataset:",targetFile[datasetStart,4],"start:",datasetStart,"end",datasetEnd,"\n",sep=" "))
    
        ##### Go through each sample in a dataset and...
        for (j in seq(datasetStart,datasetEnd-1,by=1)) { #cat(paste("\nsample:",colnames(data_combat)[j],"\n",sep=" "))

            ##### ...correlate with remaining samples from same dataset...
            for(k in seq(j+1,datasetEnd,by=1)) { #cat(paste("agaist:",colnames(data_combat)[k],"\n",sep=" "))
  
                ##### ... consider samples of the different subgroup only
                if(targets[j] != targets[k]) { #cat(paste(colnames(data_combat)[k],"group not match:",targets[j],targets[k],"\n",sep=" "))
                
                    study_corr <- c(study_corr, cor(x=data_combat[,j],y=data_combat[,k], use="all.obs", method="pearson"))
                }
            }
        }
        study_corr.list[[i]] <- study_corr
    }


    ##### Calculate Pearson correlation coefficients between the samples in all datasets
    all_corr <- NULL

    ##### Go through each sample in all datasets and...
    for (i in 1:ncol(data_combat)) { #cat(paste("\nsample:",colnames(data_combat)[i],"\n",sep=" "))
    
        j=i
        ##### ...correlate with remaining samples
        while(j < ncol(data_combat)) { #cat(paste("agaist:",colnames(data_combat)[j+1],"\n",sep=" "))
        
            j = j + 1
            all_corr <- c(all_corr, cor(x=data_combat[,i],y=data_combat[,j], use="all.obs", method="pearson"))
        }
    }

    ##### Draw density plots of Pearson correlation coefficients
    if ( is.null(group_corr.list) ) {
    
        cat("The bilogical effect cannot not be estimated!\n\n")
    
    } else if ( is.null(study_corr.list) ) {
    
        cat("The study effect cannot be estimated!\n\n")
    } else {

        density1 <- density(unlist(group_corr.list))

        density2 <- density(unlist(study_corr.list))

        density3 <- density(unlist(all_corr))

        minX <- min(unlist(group_corr.list),unlist(study_corr.list),all_corr)
        maxY <- max(density1$y,density2$y,density3$y)

        ##### Perform a two-sample Kolmogorov-Smirnov test
        KS_test <- ks.test(unlist(group_corr.list), unlist(study_corr.list), alternative = "two.sided", exact = NULL)
	
        if ( KS_test$p.value == 0 ) {
            KS_P <- " < 2.2e-16"
        } else {
            KS_P <- paste(" = ", KS_test$p.value, sep="")
        }

        pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_density_ComBat.pdf", sep="") ,width = 6, height = 6)
        plot(density1, xlim=c(minX,1), ylim=c(0,maxY), main="", xlab="Pearson correlation coefficient", col="darkolivegreen4")
        par(new=TRUE)
        plot(density2, xlim=c(minX,1), ylim=c(0,maxY), main="", xlab="", ylab="", xaxt='n', yaxt='n', col="red")
        par(new=TRUE)
        plot(density3, xlim=c(minX,1), ylim=c(0,maxY), main="", xlab="", ylab="", xaxt='n', yaxt='n',, col="black")
        legend(minX,maxY, legend=c(" ALL"," Study"," Group"),col=c("black","red","darkolivegreen4"), lty=1, box.col="black")
        text(minX+0.2, maxY, paste("p-value", KS_P, sep=""))
        dev.off()
    }

    #===============================================================================
    #     Unsupervised clustering
    #===============================================================================

    d.usa <- dist(t(data_combat), "euc")
    h.usa <- hclust(d.usa, method="ward")

    #####  Set.seed returns NULL, invisibly
    set.seed(1)

    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_cluster_datasets_ComBat.pdf", sep=""),width = clust.width, height = 6)
    h.usa$labels = colnames(data_combat)
    par(mar=c(2,2,2,6))
    A2Rplot(h.usa,
    k=length(levels(factor(datasets))), # k = changes the detail/number of subgroups shown.
    fact.sup=datasets, box=FALSE,
    show.labels=TRUE, col.up = "black", main=" ")
    dev.off()

    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_cluster_targets_ComBat.pdf", sep=""),width = clust.width, height = 6)
    h.usa$labels = colnames(data_combat)
    par(mar=c(2,2,2,6))
    A2Rplot(h.usa,
    k=length(levels(factor(datasets))), # k = changes the detail/number of subgroups shown.
    fact.sup=targets, box=FALSE,
    show.labels=TRUE, col.up = "black", main=" ")
    dev.off()
    
    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"),"_dendrogram_ComBat.pdf", sep = ""), width = 0.2*ncol(data)+2, height = 6, pointsize = 12)
    par(mar=c(2,5,2,0))
    plot(h.usa, xlab="", labels=paste(colnames(data_combat), targets, sep="       "), hang = -1, main="")
    dev.off()
    
    #===============================================================================
    #     K-means clustering
    #===============================================================================

    fit_2D <- cmdscale(d.usa,eig=TRUE, k=2) # k is the number of dim

    #####  Plot solution
    x <- fit_2D$points[,1]
    y <- fit_2D$points[,2]

    fit_2D.points = cbind(x, y)
    colnames(fit_2D.points) <- c("Coordinate 1", "Coordinate 2")
    cl <- kmeans(fit_2D.points, 1)


    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_K_means_ComBat.pdf", sep=""), width = 7.5, height = 7.5)
    plot(fit_2D.points, col = cl$cluster, type="n", ylim=c( min(y),max(y) + (abs(min(y))+abs(max(y)))/4 ))
    #####  Use different shape for each dataset
    for (i in 1:datasets_No) {
        points(x[as.numeric(factor(datasets))==i], y[as.numeric(factor(datasets))==i], cex=1, pch = pchs[[i]][1], col=targets.colour[[2]][as.numeric(factor(datasets))==i])
        points(x[as.numeric(factor(datasets))==i], y[as.numeric(factor(datasets))==i], cex=1, pch = pchs[[i]][2], col="black")
    }
    #####  Add the legend
    legend("topleft", legend=levels(factor(targets)),pch = 16, col = targets.colour[[1]], box.col = "transparent")
    legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
    dev.off()


    #===============================================================================
    #     Principal components analysis
    #===============================================================================

    #####  Keep only probes with variance > 0 across all samples
    rsd<-apply(data_combat,1,sd)
    data_combat <- data_combat[rsd>0,]

    #####  Perform principal components analysis after scaling the data
    data_combat_pca <- prcomp(t(data_combat), scale=TRUE)

    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_PCA_screeplot_ComBat.pdf", sep=""), width = 14, height = 6)
    screeplot(data_combat_pca, npcs = length(data_combat_pca$sdev), main="The variances captured by principal components")
    dev.off()

    #####  Get variance importance for all principal components
    importance_pca <- summary(data_combat_pca)$importance[2,]
    importance_pca <- paste(round(100*importance_pca, 2), "%", sep="")


    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_PCA_ComBat.pdf", sep=""))
    plot(data_combat_pca$x[,1], data_combat_pca$x[,2], type="n", xlab = paste("PC1 (",importance_pca[1],")",sep=""), ylab = paste("PC2 (",importance_pca[2],")",sep=""), ylim=c( min(data_combat_pca$x[,2]),max(data_combat_pca$x[,2]) + (abs(min(data_combat_pca$x[,2]))+abs(max(data_combat_pca$x[,2])))/4 ))

    #####  Use different shape for each dataset
    for (i in 1:datasets_No) {
        points(data_combat_pca$x[,1][as.numeric(factor(datasets))==i],data_combat_pca$x[,2][as.numeric(factor(datasets))==i], pch=pchs[[i]][1], col=targets.colour[[2]][as.numeric(factor(datasets))==i])
        points(data_combat_pca$x[,1][as.numeric(factor(datasets))==i],data_combat_pca$x[,2][as.numeric(factor(datasets))==i], pch = pchs[[i]][2], col="black")
    }
    #####  Adding the legend
    legend("topleft", legend=levels(factor(targets)),pch = 16, col =targets.colour[[1]], box.col = "transparent")
    legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
    dev.off()


    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_PCA_ComBat_labeled.pdf", sep=""))
    plot(data_combat_pca$x[,1], data_combat_pca$x[,2], type="n", xlab = paste("PC1 (",importance_pca[1],")",sep=""), ylab = paste("PC2 (",importance_pca[2],")",sep=""), ylim=c( min(data_combat_pca$x[,2]),max(data_combat_pca$x[,2]) + (abs(min(data_combat_pca$x[,2]))+abs(max(data_combat_pca$x[,2])))/4 ))

    #####  Use different shape for each dataset
    for (i in 1:datasets_No) {
        points(data_combat_pca$x[,1][as.numeric(factor(datasets))==i],data_combat_pca$x[,2][as.numeric(factor(datasets))==i], pch=pchs[[i]][1], col=targets.colour[[2]][as.numeric(factor(datasets))==i])
        points(data_combat_pca$x[,1][as.numeric(factor(datasets))==i],data_combat_pca$x[,2][as.numeric(factor(datasets))==i], pch = pchs[[i]][2], col="black")
    }
    #####  Adding sample names above symbols
    text(data_combat_pca$x[,1], data_combat_pca$x[,2]+(abs(min(data_combat_pca$x[,2]))+abs(max(data_combat_pca$x[,2])))/50, labels= colnames(data_combat), cex=0.5)
    #####  Adding the legend
    legend("topleft", legend=levels(factor(targets)),pch = 16, col = targets.colour[[1]], box.col = "transparent")
    legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
    dev.off()


    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_PCA_pairs_ComBat.pdf", sep=""))
    pairs(data_combat_pca$x[,1:4], pch=20, col=targets.colour[[2]]) # Plots all combinations between the first 4 pc
    dev.off()


    #####  Smotther scatter plot - shows the density of the data points.
    pdf(paste("Study_effect_", paste(datasets2corr, collapse="_"), "_smoothScatter_ComBat.pdf", sep=""))
    smoothScatter(data_combat_pca$x) #  A smotther scatter plot - shows the density of the data points.
    dev.off()
}

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
