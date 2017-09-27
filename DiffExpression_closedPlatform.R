################################################################################
#
#   File name: DiffExpression_closedPlatform.R
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
#   Description: Pipeline correcting for batch effect with ComBat ( http://www.bu.edu/jlab/wp-assets/ComBat/Abstract.html ) and performing differential expression analysis on microarray data. The script reads in a priori generated files with normalised expression data
#
#   Command line use example: R --file=./DiffExpression_closedPlatform.R --args "/scratch/jack/data/PhD/Transcriptomics_project" "Affy_U133Plus2" "GSE17951,GSE37199" "0.000000001" "1" "BH" "nonparam" "/scratch/jack/data/PhD/Transcriptomics_project/stromal_genes_all.txt"
#
#   First arg:      Project workspace. This is the directory to which cross-platform analysis results will be written
#   Second arg:     Platforms to be considered for the analysis. Currently suported microarray platforms are [Affy_HuEx1ST], [Affy_U133Plus2], [Affy_U133A2], [Affy_U133A], [Affy_U133B], [Affy_U95Av2], [Affy_U95B], [Affy_U95C] and [Illum_HT_12_V3]
#   Third arg:      Datasets, processed with seleceted platform in third arg to be considered for the analysis. Each datasets ID is expected to be separated by comma. Type 'all' for all listed datasets to be analysed
#   Fourth arg:     Adjustep p-value to be used as a threshold for detecting differentially expressed genes
#   Fifth arg:      Absolute value of logarithmic fold change to be used as a threshold for detecting differentially expressed genes
#   Sixth arg:      P-value adjustment method. Possible values are "none", "BH", "fdr" (equivalent to "BH"), "BY" and "holm"
#   Seventh arg (OPTIONAL): If "nonparam", use nonparametric adjustments in ComBat. AS default, parametric adjustments are used
#   Eight arg (OPTIONAL): List of genes to be excluded from the differential expression analysis. The first column is expected to list the Ensembl gene IDs
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

##### Create 'not in' operator
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0

##### Retrieve gene annotation information
annotGenes <- function (topGenes) {
	
    
    mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
    #genesAnnot <- getGene( id=rownames(topGenes), type="ensembl_gene_id", mart = mart)  ##### Currently does not work
    genesAnnot <- getBM(c("ensembl_gene_id","hgnc_symbol","description","chromosome_name","band","strand","start_position","end_position"), "ensembl_gene_id", rownames(topGenes), mart = mart)
    
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

##### Modification of the Monte Carlo integration function to find the nonparametric adjustments for a particular batch used in ComBat (see https://stat.ethz.ch/pipermail/bioconductor/2012-November/049411.html ). When using large batch sizes and the non-parametric parameter estimation some probe intensity values after ComBat correction would be NAs. To correct this, we use the log-likelihood, which is less susceptible to numerical precision errors.

int.eprior <- function(sdat,g.hat,d.hat){
    #Get the number of probes
    r <- nrow(sdat);
    #The number of probes removing the probe of interest
    G = r - 1;
    #g.star and d.star are the output matrices
    g.star <- d.star <- matrix(NA,nrow=G,ncol=1);
    #Loop over each probe
    for(i in 1:r){
        #Remove the beta value of the probe of interest
        g <- g.hat[-i];
        #Remove the variance of the probe of interest
        d <- d.hat[-i];
        #Isolate the value of this probe across all samples
        x <- sdat[i,!is.na(sdat[i,])];
        #Number of samples in this batch
        n <- length(x)
        #A matrix of dimension gxn which contains the probe values repeated in each row
        dat <- matrix(as.numeric(x),length(g),n,byrow=T)
        #The residual between the probe values and the regression coefficient
        resid2 <- (dat-g)^2
        #Calculate Log-likelihood using Gaussian
        LL = matrix(0,nrow=length(G),ncol=1);
        for (j in 1:n) {
            LL = LL + -0.5*log(2*pi*d)-(resid2[,j]/(2*d));
        }
        #Order by the log-likelihood values to avoid Inf in parameter estimates
        orderind = order(LL,decreasing=TRUE);
        LL = LL[orderind];
        g = g[orderind];
        d = d[orderind];
        #Determine the parameter estimates
        numerator.g = 1;
        denominator = 1;
        numerator.d = 1;
        for (k in 2:G) {
            numerator.g = numerator.g + (g[k]/g[1]) * exp(LL[k]-LL[1]);
            denominator = denominator + exp(LL[k]-LL[1]);
            numerator.d = numerator.d + (d[k]/d[1]) * exp(LL[k]-LL[1]);
        }
        g.star[i] = g[1]*(numerator.g/denominator);
        d.star[i] = d[1]*(numerator.d/denominator);
    }
    adjust <- rbind(g.star,d.star)
    rownames(adjust) <- c("g.star","d.star")
    adjust
}

#===============================================================================
#    Load libraries
#===============================================================================

library(sva)
library(limma)
library(gplots)
library(biomaRt)
source("a2R_code.R")

##### Override the original int.eprior function in the svn package namespace with our functio based on log-likelihood (http://stackoverflow.com/questions/8743390/how-do-i-override-a-non-visible-function-in-the-package-namespace )
assignInNamespace("int.eprior", int.eprior, ns="sva")

#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

ProjectDir = args[4]
platform = args[5]
datasets2corr = args[6]
datasets2corr = gsub("\\s","", datasets2corr)
datasets2corr =  sort(unlist(strsplit(datasets2corr, split=',', fixed=TRUE)))
pThreshold = as.numeric(args[7])
lfcThreshold = as.numeric(args[8])
adjMethod = args[9]

if (!is.na(args[10]) && args[10] == "nonparam") {
    
	parPrior = FALSE
    
    if (!is.na(args[11])) {
        
        genesList = args[11]
        
    } else {
        
        genesList = NA
    }
    
} else if (!is.na(args[10]) && args[10] != "nonparam") {
    
    parPrior = TRUE
	genesList = args[10]
    
} else {
    
    parPrior = TRUE
    genesList = NA
}


##### Read file with datasets information
DatasetInput=read.table(paste(ProjectDir,"GenExpression_InputFiles.txt",sep="/"),sep="\t",as.is=TRUE,header=TRUE,row.names=1)
DatasetInput=DatasetInput[DatasetInput[,"Platform"] %in% platform,]


##### Read file with genes to exclude if provided
if (!is.na(genesList)) {
    
    genes2excl=read.table(genesList,sep="\t",as.is=TRUE,header=TRUE)
    genes2excl <- unique(genes2excl[,1])
    
} else {
	genes2excl = NULL
}


if ( datasets2corr != "all" && datasets2corr != "ALL" ) {
    
    if ( all(datasets2corr %in% rownames(DatasetInput)) ) {
        DatasetInput=DatasetInput[rownames(DatasetInput) %in% datasets2corr,]
        
    } else{
        cat(paste("\nSome of selected datasets (", paste(datasets2corr, collapse = ", "), ") were not generated by selected platforms (", paste(platform, collapse = ", "), ")!\n\n", sep=" "))
        q()
    }
} else {
    datasets2corr <- sort(rownames(DatasetInput))
}

##### Change working directory to the project workspace
setwd(ProjectDir)

##### Report used parameters to a file
write(args, file = paste("Comb", paste(datasets2corr, collapse="_"), "parameters.txt", sep = "_"), append = FALSE, sep="\t")


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

if ( unique(DatasetInput[,"Platform"]) == "Affy_HuEx1ST" ) {
    
    datasetList[[1]] <- paste(ProjectDir,"Affy_HuEx1ST_qn.exp",sep="/")
    #filterThreshold <- 10000
    
} else if ( unique(DatasetInput[,"Platform"]) == "Affy_U133Plus2" ) {
    
    datasetList[[1]] <- paste(ProjectDir,"Affy_U133Plus2_gcrma.exp",sep="/")
    #filterThreshold <- 10000
    
} else if ( unique(DatasetInput[,"Platform"]) == "Affy_U133A2" ) {
    
    datasetList[[1]] <- paste(ProjectDir,"Affy_U133A2_gcrma.exp",sep="/")
    #filterThreshold <- 5000
    
} else if ( unique(DatasetInput[,"Platform"]) == "Affy_U133A" ) {
    
    datasetList[[1]] <- paste(ProjectDir,"Affy_U133A_gcrma.exp",sep="/")
    #filterThreshold <- 5000
    
} else if ( unique(DatasetInput[,"Platform"]) == "Affy_U133B" ) {
    
    datasetList[[1]] <- paste(ProjectDir,"Affy_U133B_gcrma.exp",sep="/")
    #filterThreshold <- 5000
    
} else if ( unique(DatasetInput[,"Platform"]) == "Affy_U95Av2" ) {
    
    datasetList[[1]] <- paste(ProjectDir,"Affy_U95Av2_gcrma.exp",sep="/")
    #filterThreshold <- 3000
    
} else if ( unique(DatasetInput[,"Platform"]) == "Affy_U95B" ) {
    
    datasetList[[1]] <- paste(ProjectDir,"Affy_U95B_gcrma.exp",sep="/")
    #filterThreshold <- 3000
    
} else if ( unique(DatasetInput[,"Platform"]) == "Affy_U95C" ) {
    
    datasetList[[1]] <- paste(ProjectDir,"Affy_U95C_gcrma.exp",sep="/")
    #filterThreshold <- 3000
    
} else if ( unique(DatasetInput[,"Platform"]) == "Illum_HT_12_V3" ) {
    
    datasetList[[1]] <- paste(ProjectDir,"Illum_HT_12_V3_rsn.exp",sep="/")
    #filterThreshold <- 10000
    
} else {
    cat(paste("Platform", unique(DatasetInput[,"Platform"]) ,"is not yet supported!\n", sep=" "))
    q()
}


if ( length(datasetList) == 1 ) {
     cat(paste("Loading expression data", datasetList[[1]], "...\n"), sep="")
     data <- read.table(datasetList[[1]],sep="\t",as.is=TRUE,header=TRUE,row.names=1,check.names=FALSE)
} else {
	cat("No datasets provided!\n")
    q()
}
dim(data)


##### Exclude user-defined genes
if (!is.na(genesList)) {
    
    data <- data[rownames(data) %!in% genes2excl,]
}


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
#     Remove batch effects with ComBat (if more than one dataset selected)
#===============================================================================

if ( datasets_No > 1 ) {
    
    #####  Keep only probes with variance > 0 across all samples
    rsd<-apply(data,1,sd)
    data <- data[rsd>0,]

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
		##### Generate prior plots to assess if parametric or nonparametric adjustments should be used in ComBat. If the red lines (parametric estimate of the empirical batch effect density) and black lines (kernel estimate of the empirical batch effect density) don't match up well, use the nonparametric adjustments
		pdf(file = paste("Comb_", paste(datasets2corr, collapse="_"),"_PriorPlot.pdf", sep = ""))
    	data_combat = ComBat(dat=data, batch=batch, mod=mod, par.prior=parPrior, prior.plots=TRUE)
		dev.off()
	} else {
		data_combat = ComBat(dat=data, batch=batch, mod=mod, par.prior=parPrior, prior.plots=FALSE)
	}
	
    cat(paste("Writing combined expression data to ", paste("Comb_", paste(datasets2corr, collapse="_"), ".exp", sep=""), "\n",sep=""))
    write.table(prepare2write(data_combat), file=paste("Comb_", paste(datasets2corr, collapse="_"), ".exp", sep=""),sep="\t", row.names=FALSE)

} else {
    
    data_combat <- data
    write.table(prepare2write(data_combat), file=paste("Comb_", paste(datasets2corr, collapse="_"), ".exp", sep=""),sep="\t", row.names=FALSE)
}

#===============================================================================
#     Non-specific filtering
#===============================================================================

##### Use 60% of all genes with the highest expression variance across samples in the non-specific filtering step
filterThreshold <- round(0.6*nrow(data_combat), digits = 0)


##### Depending on the platform, keep 60% of all genes with the highest expression variance across samples
rsd<-apply(data_combat,1,sd)
sel<-order(rsd, decreasing=TRUE)[1:filterThreshold]
data_combat <-data_combat[sel,]

#===============================================================================
#     Differenatial expression analysis
#===============================================================================

#####  Estimate relative quality weights for each array
cat("Estimating relative quality weights for each array...\n")
arrayw <- arrayWeights(data_combat)

pdf(paste("Comb_", paste(datasets2corr, collapse="_"),"_arrayWeights.pdf", sep = ""), width = 0.2*ncol(data)+2, height = 6, pointsize = 12)
par(las=3, mar=c(11,5,1,0))
barplot(arrayw, xlab="", ylab="Weight", col=datasets.colour[[2]], names=colnames(data_combat), las=2)
abline(h=1, lwd=1, lty=2)
dev.off()

if (length(unique(targets)) > 1) {
	
	cat(paste("There are", length(unique(targets)), "groups specified in this study. Perfroming differenatial expression analysis using limma...\n"), sep=" ")
} else {
	cat("There is only one group specified in this study. No differential expression analysis will be performed!\n")
	q()
}

#####  Create model matrix for outcome of interest
f.model <- factor(targets, levels = levels(factor(targets)))
mod <- model.matrix(~0 + f.model)
colnames(mod) <- levels(factor(targets))

#####  Add to the model relevant principle components to adjust the data for surrogate variables detected in PCA analysis -> divide this code to two parts, data combination and differential expression analysis to enable users to select PCs to be used to adjust the data for surrogate variables
#data_combat_pca$rotation
#data_combat_pca$sdev
#mod = cbind(mod,data_combat_pca$sdev)
#mod <- model.matrix(~0 + f.model + data_combat_pca$sdev)
#colnames(mod) <- c(levels(factor(targets)),"PC1")

#####  Fit linear model to combined data, given the above design
if ( all(is.na(targetFile$Replicate)) ) {
    
    fit <- lmFit(data_combat, mod, weights=arrayw)
    
} else {
    cat("Dealing with technical replicates...\n")
    
    biolrep <- targetFile$Replicate
    biolrep[is.na(biolrep)] <- 0
    maxRep <- max(biolrep)
    
    cumSum = maxRep
    for (i in 1:length(biolrep)) {
        if ( biolrep[i]==0 ) {
            cumSum = cumSum+1
            biolrep[i] <- cumSum
        }
    }

    corfit <- duplicateCorrelation(data_combat, design=mod, ndups = 1, block = biolrep)
    fit <- lmFit(data_combat, mod, block = biolrep, cor = corfit$consensus, weights=arrayw)
}

#####  Create matrix of possible comparisons
comb <- combn(levels(factor(targets)), 2)

#####  Get number of possible comparisons using the following formula:
#
# n!/((n-r)!(r!))
#
# n = the number of classes to compare
# r = the number of elements for single comparison
#
################################################################################

targetsNo <- length(levels(factor(targets)))
combNo <- factorial(targetsNo)/(factorial(targetsNo-2)*(factorial(2))) # n!/((n-r)!(r!))

contrasts <- NULL
contrastNames <- NULL

#####  Create string with possible contrasts
for (i in 1:combNo) {
    
    contrasts <- c(contrasts, paste(paste(comb[1,i], comb[2,i], sep="vs"),paste(comb[1,i], comb[2,i], sep="-"), sep="="))
    contrastNames[i] <- paste(comb[1,i], comb[2,i], sep=" vs ")
}
contrasts <- paste(contrasts, collapse=", ")

#####  Create contrasts of interest
func = "makeContrasts"
arguments = paste(contrasts, "levels=mod",sep=", ")

contrast.matrix <- eval(parse(text = paste(func, "(", arguments, ")",sep="")))

#####  Fit contrasts to linear model
fitContrasts <- contrasts.fit(fit, contrast.matrix)

#####  Apply empirical Bayes statistics
eb <- eBayes(fitContrasts)


#===============================================================================
#     Differenatial expression analysis results visualisation
#===============================================================================

#####  P-value histograms
for (i in 1:ncol(eb$p.value) ) {
    pdf(file = paste("Comb_", paste(datasets2corr, collapse="_"), "_", colnames(eb$p.value)[i], "_P_hist.pdf", sep=""))
    histogram <- hist(eb$p.value[,i], breaks=seq(0,1,by= 0.01), main=contrastNames[i], xlab="p-value")
    exprected_p.value <- mean(histogram$counts)+(sd(histogram$counts)*1)
    abline(v=0.05,col="red")
    abline(h=exprected_p.value,col="blue")
    dev.off()
}

#####  Volcano plots of log2 fold-changes versus significance (adjusted p-values, +label top 10)
for (i in 1:ncol(eb$p.value) ) {
    
    topGenes <- topTable(eb, coef=colnames(eb)[i], adjust=adjMethod, sort.by="none", number=nrow(data_combat))
    
    pdf(file = paste("Comb_", paste(datasets2corr, collapse="_"), "_", colnames(eb$p.value)[i], "_volcano_plot", ".pdf", sep = ""))
    plot(topGenes[,"logFC"],-log2(topGenes[,"adj.P.Val"]),pch=16,cex=0.5,xlab="Log2 fold-change",ylab="-log2(adjusted p-value)",main=contrastNames[i],col="grey")
    #####  Highlight genes with logFC above specified threshold
    points(topGenes[abs(topGenes[,"logFC"])>lfcThreshold,"logFC"],-log2(topGenes[abs(topGenes[,"logFC"])>lfcThreshold,"adj.P.Val"]),cex=0.5,pch=16)
    abline(h=-log2(pThreshold),col="red", lty = 2)
    #####  Label top 10 most significant genes
    ord <- order(-log2(topGenes[,"adj.P.Val"]),decreasing=TRUE)
    top10 <- ord[1:10]
    text(topGenes[top10,"logFC"],-log2(topGenes[top10,"adj.P.Val"]), labels=rownames(data_combat[top10,]),cex=0.6,col="blue")
    dev.off()
}

eb.decideTests <- decideTests(eb, method="separate", adjust.method=adjMethod, p.value=pThreshold, lfc=lfcThreshold)
summary(eb.decideTests)

write.table(prepare2write(eb.decideTests), file=paste("Comb_", paste(datasets2corr, collapse="_"),"_decideTests.txt", sep = ""), sep="\t", row.names=FALSE)


#####  Generate Venn diagrams
venn.counts <- vennCounts(eb.decideTests, include = "both")
colnames(venn.counts)[1:ncol(venn.counts)-1] <- contrastNames

if ( length(contrastNames) < 4) {
    pdf(file = paste("Comb_", paste(datasets2corr, collapse="_"),"_DE_venn.pdf", sep = ""), width = 8, height = 8, pointsize = 7, bg = "transparent")
    vennDiagram(venn.counts)
    dev.off()
}


#####  Write the analysis results into a file
for (i in 1:ncol(eb$p.value) ) {
    
	topGenes <- topTable(eb, coef=colnames(eb)[i], adjust=adjMethod,  sort.by="p", number=nrow(data_combat))
    topGenes <- topGenes[,colnames(topGenes) %in% c("logFC","t","P.Value","adj.P.Val")]
    colnames(topGenes) <- c("log2FC", "t-statistic", "p-value", paste("p-value", adjMethod,sep="_"))

    cat(paste("Writing differential expression analysis results for", contrastNames[i], "comparison to", paste("Comb_", paste(datasets2corr, collapse="_"), "_", colnames(eb$p.value)[i], "_topTable.txt", sep=""), "\n", sep=" "))
    ##### Retrieve gene annotation information
    topGenes <- annotGenes(topGenes)
    
    write.table(prepare2write(topGenes), file=paste("Comb_", paste(datasets2corr, collapse="_"), "_", colnames(eb$p.value)[i], "_topTable.txt", sep=""),sep="\t", row.names=FALSE)
}

#####  ... write into a file only differnetially expressed genes
topGenes.index <- NULL

for (i in 1:ncol(eb$p.value) ) {
    
	topGenes <- topTable(eb, coef=colnames(eb)[i], adjust=adjMethod,  sort.by="p", p.value=pThreshold, lfc=lfcThreshold, number=nrow(data_combat))
    
    if ( length(topGenes) != 0 ) {
        
        topGenes <- topGenes[,colnames(topGenes) %in% c("logFC","t","P.Value","adj.P.Val")]
        colnames(topGenes) <- c("log2FC", "t-statistic", "p-value", paste("p-value", adjMethod,sep="_"))
        
        cat(paste("Writing differential expression analysis results for", contrastNames[i], "comparison to", paste("Comb_", paste(datasets2corr, collapse="_"), "_", colnames(eb$p.value)[i], "_DE_genes.txt", sep=""), "\n", sep=" "))
        ##### Retrieve gene annotation information
        topGenes <- annotGenes(topGenes)
        
        write.table(prepare2write(topGenes), file=paste("Comb_", paste(datasets2corr, collapse="_"), "_", colnames(eb$p.value)[i], "_DE_genes.txt", sep=""),sep="\t", row.names=FALSE)
        
        topGenes.index <- c(topGenes.index,rownames(topGenes))
    } else {
        cat(paste("No differentially expressed genes detected in ", contrastNames[i], "\n", sep=""))
    }
}


#===============================================================================
#     Hierarchical clustering
#===============================================================================

#####  Select top discriminating genes
topGenes <- data_combat[unique(topGenes.index),]
dim(topGenes)

if (nrow(topGenes) > 0) {
    hr <- hclust(as.dist(1-cor(t(topGenes), method="pearson")), method="ward")
    hc <- hclust(as.dist(dist(t(topGenes), method="euclidean")), method="ward")

    pdf(paste("Comb_", paste(datasets2corr, collapse="_"),"_DE_dendrogram.pdf", sep = ""), width = 0.2*ncol(data)+2, height = 6, pointsize = 12)
    par(mar=c(2,5,2,0))
    plot(hc, xlab="", labels=paste(colnames(data_combat), targets, sep="       "), hang = -1, main="")
    dev.off()

    pdf(paste("Comb_", paste(datasets2corr, collapse="_"),"_DE_heatmap.pdf", sep = ""), width = 6, height = 10, pointsize = 12)
    heatmap.2(as.matrix(topGenes), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), scale="row", ColSideColors=targets.colour[[2]] ,margins = c(2, 6), labRow="", labCol="", trace="none", key = FALSE)
    #####  Add the legend
    legend("topright", legend=levels(factor(targets)),fill=targets.colour[[1]], box.col = "transparent", cex=1.2)
    dev.off()

    #####  Set.seed returns NULL, invisibly
    set.seed(1)

    pdf(paste("Comb_", paste(datasets2corr, collapse="_"),"_DE_cluster.pdf", sep = ""),width = 0.1*ncol(data)+2, height = 6)
    hc$labels = colnames(data)
    par(mar=c(2,2,2,6))
    A2Rplot(hc,
    k=targetsNo, # k = changes the detail/number of subgroups shown.
    fact.sup=targets, box=FALSE,
    show.labels=TRUE, col.up = "black", main=" ")
    dev.off()

    #===============================================================================
    #     Pair-wise comparison
    #===============================================================================

    distance <- as.dist(dist(t(topGenes), method="euclidean"))
    plot.top.matrix <- as.matrix(distance)
    hc <- hclust(distance, method="ward")

    pdf(paste("Comb_", paste(datasets2corr, collapse="_"),"_DE_pairwise.pdf", sep = ""), width = 9, height = 8, pointsize = 12)
    heatmap.2(plot.top.matrix, symm = TRUE, distfun = as.dist, ColSideColors=targets.colour[[2]] , RowSideColors=targets.colour[[2]] , Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc),margins = c(2, 6),labCol="",labRow="", trace="none", key = FALSE)
    #####  Add the legend
    legend("topright", legend=levels(factor(targets)),fill=targets.colour[[1]], box.col = "transparent", cex=1.2)
    dev.off()
}

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()
