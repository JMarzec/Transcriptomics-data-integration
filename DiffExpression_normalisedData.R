################################################################################
#
#   File name: DiffExpression_normalisedData.R
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
#   Description: Pipeline performing differential expression analysis on normalised expression data. NOTE: the script allowes to process gene matrix with duplicated gene IDs.
#
#   Command line use example: R --file=./DiffExpression_normalisedData.R --args "/Users/marzec01/data/PhD/Transcriptomics_project_HGPIN_sign/GSE6099_expr_matrix.txt" "/Users/marzec01/data/PhD/Transcriptomics_project_HGPIN_sign/GSE6099_annot.txt" "Normal_Epithelium,BPH_Epithelium,PIN,Localised_PCa,Metastatic_PCa" "0.05" "1" "BH" "/Users/marzec01/data/PhD/Transcriptomics_project_HGPIN_sign/GSE6099_HGPIN_sign_DE" "/Users/marzec01/data/PhD/Transcriptomics_project_HGPIN_sign/Meta_top500_tracking_4_genes_symbols.txt"
#
#   First arg:      Full path with name of the normalised expression matrix of the dataset of interest
#   Second arg:     Full path with name of the text file with samples annotation. The file is expected to include the following columns: (1) sample name, (2) file names and (3) sample type
#   Third arg:      String of sample types to be used for the analysis (e.g. Normal, Tumour, Metastatic)
#   Fourth arg:     Adjustep p-value to be used as a threshold for detecting differentially expressed genes
#   Fifth arg:      Absolute value of logarithmic fold change to be used as a threshold for detecting differentially expressed genes
#   Sixth arg:      P-value adjustment method. Possible values are "none", "BH", "fdr" (equivalent to "BH"), "BY" and "holm"
#   Seventh arg:    Full path with name of the output folder
#   Eight arg (OPTIONAL): Full path with name of the text file listing genes to be considered for the analysis. Individual genes are expected to be listed per row
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


##### Assign colours to analysed groups
getTargetsColours <- function(targets) {
    
    ##### Predefined selection of colours for groups
    #targets.colours <- c("blue","red","darkgoldenrod","green","darkred","deepskyblue")
    targets.colours <- c("darkred","red","darkgoldenrod","darkgreen","darkolivegreen3","burlywood","blue","lightskyblue")
    
    f.targets <- factor(targets)
    vec.targets <- targets.colours[1:length(levels(f.targets))]
    targets.colour <- rep(0,length(f.targets))
    for(i in 1:length(f.targets))
    targets.colour[i] <- vec.targets[ f.targets[i]==unique(f.targets)]
    
    return( list(vec.targets, targets.colour) )
}

##### Create 'not in' operator
"%!in%" <- function(x,table) match(x,table, nomatch = 0) == 0


#===============================================================================
#    Load libraries
#===============================================================================

library(limma)
library(Amelia)
library(gplots)
source("a2R_code.R")


#===============================================================================
#    Main
#===============================================================================

args <- commandArgs()

expFile <- args[4]
annFile <- args[5]
groups <- args[6]
groups <- gsub("\\s","", groups)
groups <-  unlist(strsplit(groups, split=',', fixed=TRUE))
pThreshold = as.numeric(args[7])
lfcThreshold = as.numeric(args[8])
adjMethod = args[9]
outFolder <- args[10]


##### Read file with normalised expression data
expData <- read.table(expFile,sep="\t",as.is=TRUE,header=TRUE,row.names=NULL)

##### Retieve the expression data file name
coreName <- strsplit(expFile, "/")
coreName <- coreName[[1]][length(coreName[[1]])]


##### Read sample annotation file
annData <- read.table(annFile,sep="\t",as.is=TRUE,header=TRUE,row.names=1)
rownames(annData) <- gsub("-", ".", rownames(annData))


##### Keep only annotated samples
expData <- expData[,colnames(expData) %in% c(colnames(expData)[1],rownames(annData))]
annData <- subset(annData, rownames(annData) %in% colnames(expData))


##### Match samples in expression matrix and sample annotation file
annData <- annData[ colnames(expData)[2:ncol(expData)], ]


##### Keep only samples assigned to the groups of interest
expData <- expData[, c(TRUE, annData[,2] %in% groups) ]
annData <- subset(annData, rownames(annData) %in% colnames(expData))


##### Keep only genes of interest if such list is provided
if ( !is.na(args[11]) ) {
    
    ##### Read file with the gene list
    genes <- unlist(read.table(args[11],sep="\t",as.is=TRUE,header=FALSE))
    
    ##### Identify genes of interest not present in the expression matrix
    absentGenes <- genes[genes %!in% expData[,1]]
    
    ##### Keep only genes of interest
    expData <- expData[ expData[,1] %in% genes, ]
}
    

#===============================================================================
#    Imputation of missing data using Ameila package
#===============================================================================

##### Remove genes with missing data in >50% of samples
genes2remove <- NULL
allowedNAs <- round(0.5*ncol(expData), digits = 0)

for ( i in 1:nrow(expData) ) {
    
    if ( sum(is.na(expData[i,])) > allowedNAs ) {
        
        genes2remove <- c(genes2remove, i)
    }
}

if ( !is.null(genes2remove) ) {
    expData <- expData[-genes2remove,]
}


##### Extract expression matrix gene list
expData_annot <- expData[,1]
expData <- expData[,-1]


##### Use amelia function to impute missing date. NOTE: the parameter 'empri' is set to the number of rows of the imput data. If error is produces by 'amelia' funciton try to increase this number
##### First check if there is any missing data
if ( any(is.na(expData)) ) {
    
    expDataImput <- amelia(x = expData, empri = nrow(expData))
    expData <- expDataImput$imputations$imp1
}

#===============================================================================


##### Reverse the groups so that the later groups are compared to the prior groups in the DE analysis step
groups <- rev(groups)

##### Order samples in the expression matrix accrodingly to provided sample type oreder
samplesOrdered <- NULL

for ( i in 1:length(groups) ) {
    
    for ( j in 1:nrow(annData) ) {
        
        if ( annData[j,2] == groups[i] ) {
            
            samplesOrdered <- c(samplesOrdered,j)
            
        }
    }
}

expData <- expData[,samplesOrdered]
annData <- annData[samplesOrdered,]

##### Set/create a directory for the output files
if (file.exists(outFolder)){
    cat( paste("The output folder already exists. Writing files to the", outFolder, "\n", sep = " ") )
} else {
    dir.create(outFolder);
    cat( paste("Writing files to the", outFolder, "\n", sep = " ") )
}

##### Change working directory to the project workspace
setwd(outFolder)

##### Report used parameters to a file
write(args, file = paste(coreName, "parameters.txt", sep = "_"), append = FALSE, sep="\t")


#===============================================================================
#     Non-specific filtering
#===============================================================================

##### If not gene list is provided then use 60% of all genes with the highest expression variance across samples for the non-specific filtering step
if ( is.na(args[11]) ) {
    ##### Use 60% of all genes with the highest expression variance across samples in the non-specific filtering step
    filterThreshold <- round(0.6*nrow(expData), digits = 0)

    ##### Depending on the platform, keep 60% of all genes with the highest expression variance across samples
    rsd<-apply(expData,1,sd)
    sel<-order(rsd, decreasing=TRUE)[1:filterThreshold]
    expData <-expData[sel,]

##### ...otherwise use all use-defined gene list
} else {
    
    ##### Report genes of interest not present in the expression matrix
    write(absentGenes, file = paste(coreName, "absent_genes.txt", sep = "_"), append = FALSE, sep="\t")
}


#===============================================================================
#     Unsupervised clustering
#===============================================================================

d.usa <- dist(t(expData), "euc")
h.usa <- hclust(d.usa, method="ward")

#####  Adopt cluster plots width
if ( ncol(expData) < 20 ) {
    
    clust.width = 0.3*ncol(expData)
    
} else {
    clust.width = 0.1*ncol(expData)
}


#####  Set.seed returns NULL, invisibly
set.seed(1)

pdf(paste(coreName, "_cluster.pdf", sep=""),width = clust.width, height = 6)
h.usa$labels = colnames(expData)
par(mar=c(2,2,2,6))
A2Rplot(h.usa,
k=length(groups), # k = changes the detail/number of subgroups shown.
fact.sup=annData[,2], box=FALSE,
show.labels=TRUE, col.up = "black", main=" ")
dev.off()


pdf(paste(coreName, "_dendrogram.pdf", sep = ""), width = 0.2*ncol(expData)+2, height = 6, pointsize = 12)
par(mar=c(2,5,2,0))
plot(h.usa, xlab="", labels=paste(colnames(expData), annData[,2], sep="       "), hang = -1, main="")
dev.off()


#===============================================================================
#     K-means clustering
#===============================================================================

##### Prepare colours for sample groups
targets.colour <- getTargetsColours(annData[,2])


fit_2D <- cmdscale(d.usa,eig=TRUE, k=2) # k is the number of dim

#####  Plot solution
x <- fit_2D$points[,1]
y <- fit_2D$points[,2]

fit_2D.points = cbind(x, y)
colnames(fit_2D.points) <- c("Coordinate 1", "Coordinate 2")
cl <- kmeans(fit_2D.points, 1)


pdf(paste(coreName, "_K_means.pdf", sep=""), width = 7.5, height = 7.5)
plot(fit_2D.points, col = cl$cluster, type="n", ylim=c( min(y),max(y) + (abs(min(y))+abs(max(y)))/4 ))

points(x, y, cex=1, pch = 16, col=targets.colour[[2]])
points(x, y, cex=1, pch = 1, col="black")

#####  Add the legend
legend("topleft", legend=groups,pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()


#===============================================================================
#     Principal components analysis
#===============================================================================

#####  Keep only probes with variance > 0 across all samples
rsd<-apply(expData,1,sd)
expData <- expData[rsd>0,]

#####  Perform principal components analysis after scaling the data
data_pca <- prcomp(t(expData), scale=TRUE)


pdf(paste(coreName, "_PCA_screeplot.pdf", sep=""), width = 14, height = 6)
screeplot(data_pca, npcs = length(data_pca$sdev), main="The variances captured by principal components")
dev.off()

#####  Get variance importance for all principal components
importance_pca <- summary(data_pca)$importance[2,]
importance_pca <- paste(round(100*importance_pca, 2), "%", sep="")


pdf(paste(coreName, "_PCA.pdf", sep=""))
plot(data_pca$x[,1], data_pca$x[,2], type="n", xlab = paste("PC1 (",importance_pca[1],")",sep=""), ylab = paste("PC2 (",importance_pca[2],")",sep=""), ylim=c( min(data_pca$x[,2]),max(data_pca$x[,2]) + (abs(min(data_pca$x[,2]))+abs(max(data_pca$x[,2])))/4 ))

points(data_pca$x[,1],data_pca$x[,2], pch=16, col=targets.colour[[2]])
points(data_pca$x[,1],data_pca$x[,2], pch = 1, col="black")

#####  Adding the legend
legend("topleft", legend=groups,pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()


pdf(paste(coreName, "_PCA_labeled.pdf", sep=""))
plot(data_pca$x[,1], data_pca$x[,2], type="n", xlab = paste("PC1 (",importance_pca[1],")",sep=""), ylab = paste("PC2 (",importance_pca[2],")",sep=""), ylim=c( min(data_pca$x[,2]),max(data_pca$x[,2]) + (abs(min(data_pca$x[,2]))+abs(max(data_pca$x[,2])))/4 ))

points(data_pca$x[,1],data_pca$x[,2], pch=16, col=targets.colour[[2]])
 points(data_pca$x[,1],data_pca$x[,2], pch = 1, col="black")

#####  Adding sample names above symbols
text(data_pca$x[,1], data_pca$x[,2]+(abs(min(data_pca$x[,2]))+abs(max(data_pca$x[,2])))/50, labels= colnames(expData), cex=0.5)
#####  Adding the legend
legend("topleft", legend=groups,pch = 16, col = targets.colour[[1]], box.col = "transparent")
legend("topleft", legend=c(rep("",length(targets.colour[[1]]))), pch = 1, box.col = "transparent")
dev.off()


pdf(paste(coreName, "_PCA_pairs.pdf", sep=""))
pairs(data_pca$x[,1:4], pch=20, col=targets.colour[[2]]) # Plots all combinations between the first 4 pc
dev.off()


#####  Smotther scatter plot - shows the density of the data points
pdf(paste(coreName, "_smoothScatter.pdf", sep=""))
smoothScatter(data_pca$x) #  A smotther scatter plot - shows the density of the data points.
dev.off()


#===============================================================================
#     Differenatial expression analysis
#===============================================================================

##### Deal with the duplicated genes
genesList <- NULL
genesRepl <- NULL

for ( j in 1:nrow(expData) ) {
    
    geneName <- expData_annot[j]
    
    ##### Distingish duplicated genes by adding duplicate number
    if ( geneName %in% genesList ) {
        
        ##### Report genes with more than one duplicates
        if ( geneName %in% names(genesRepl) ) {
            
            genesRepl[[ geneName ]] = genesRepl[[ geneName ]]+1
            
            geneName <- paste(geneName, ".", genesRepl[[ geneName ]], sep="")
            
        } else {
            genesRepl[[ geneName ]] <- 2
            
            geneName <- paste(geneName, ".2", sep="")
        }
    }
    genesList <- c(genesList,geneName)
}

rownames(expData) <- genesList


#####  Estimate relative quality weights for each array
cat("Estimating relative quality weights for each array...\n")
arrayw <- arrayWeights(expData)

pdf(paste(coreName, "_arrayWeights.pdf", sep = ""), width = 0.2*ncol(expData)+2, height = 6, pointsize = 12)
par(las=3, mar=c(11,5,1,0))
barplot(arrayw, xlab="", ylab="Weight", col=targets.colour[[2]], names=colnames(expData), las=2, ylim=c(0, max(arrayw)+0.5))
abline(h=1, lwd=1, lty=2)
legend("topleft", legend=groups,fill=targets.colour[[1]], box.col = "transparent", bty = "n", horiz=TRUE, cex=1)
dev.off()


if (length(groups) > 1) {
	
	cat(paste("There are", length(groups), "groups specified in this study. Perfroming differenatial expression analysis using limma...\n"), sep=" ")
} else {
	cat("There is only one group specified in this study. No differential expression analysis will be performed!\n")
	q()
}

#####  Create model matrix for outcome of interest
f.model <- factor(annData[,2], levels = levels(factor(annData[,2])))
mod <- model.matrix(~0 + f.model)
colnames(mod) <- levels(factor(groups))

#####  Add to the model relevant principle components to adjust the data for surrogate variables detected in PCA analysis
#data_combat_pca$rotation
#data_combat_pca$sdev
#mod = cbind(mod,data_combat_pca$sdev)
#mod <- model.matrix(~0 + f.model + data_combat_pca$sdev)
#colnames(mod) <- c(levels(factor(targets)),"PC1")

#####  Fit linear model to combined data, given the above design
fit <- lmFit(expData, mod, weights=arrayw)
    

#####  Create matrix of possible comparisons
comb <- combn(groups, 2)

#####  Get number of possible comparisons using the following formula:
#
# n!/((n-r)!(r!))
#
# n = the number of classes to compare
# r = the number of elements for single comparison
#
################################################################################

groupsNo <- length(groups)
combNo <- factorial(groupsNo)/(factorial(groupsNo-2)*(factorial(2))) # n!/((n-r)!(r!))

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
    pdf(file = paste(coreName, "_", colnames(eb$p.value)[i], "_P_hist.pdf", sep=""))
    histogram <- hist(eb$p.value[,i], breaks=seq(0,1,by= 0.01), main=contrastNames[i], xlab="p-value")
    exprected_p.value <- mean(histogram$counts)+(sd(histogram$counts)*1)
    abline(v=0.05,col="red")
    abline(h=exprected_p.value,col="blue")
    dev.off()
}

#####  Volcano plots of log2 fold-changes versus significance (adjusted p-values, +label top 10)
for (i in 1:ncol(eb$p.value) ) {
    
    topGenes <- topTable(eb, coef=colnames(eb)[i], adjust=adjMethod, sort.by="none", number=nrow(expData))
    
    pdf(file = paste(coreName, "_", colnames(eb$p.value)[i], "_volcano_plot", ".pdf", sep = ""))
    plot(topGenes[,"logFC"],-log2(topGenes[,"adj.P.Val"]),pch=16,cex=0.5,xlab="Log2 fold-change",ylab="-log2(adjusted p-value)",main=contrastNames[i],col="grey")
    #####  Highlight genes with logFC above specified threshold
    points(topGenes[abs(topGenes[,"logFC"])>lfcThreshold,"logFC"],-log2(topGenes[abs(topGenes[,"logFC"])>lfcThreshold,"adj.P.Val"]),cex=0.5,pch=16)
    abline(h=-log2(pThreshold),col="red", lty = 2)
    #####  Label top 10 most significant genes
    ord <- order(-log2(topGenes[,"adj.P.Val"]),decreasing=TRUE)
    top10 <- ord[1:10]
    text(topGenes[top10,"logFC"],-log2(topGenes[top10,"adj.P.Val"]), labels=rownames(expData)[top10],cex=0.6,col="blue")
    dev.off()
}

eb.decideTests <- decideTests(eb, method="separate", adjust.method=adjMethod, p.value=pThreshold, lfc=lfcThreshold)
summary(eb.decideTests)

write.table(prepare2write(eb.decideTests), file=paste(coreName,"_decideTests.txt", sep = ""), sep="\t", row.names=FALSE)


#####  Generate Venn diagrams
if ( combNo < 4) {
    
    venn.counts <- vennCounts(eb.decideTests, include = "both")
    colnames(venn.counts)[1:ncol(venn.counts)-1] <- contrastNames
    
    pdf(file = paste(coreName,"_DE_venn.pdf", sep = ""), width = 8, height = 8, pointsize = 7, bg = "transparent")
    vennDiagram(venn.counts)
    dev.off()
}


#####  Write the analysis results into a file
for (i in 1:ncol(eb$p.value) ) {
    
	topGenes <- topTable(eb, coef=colnames(eb)[i], adjust=adjMethod,  sort.by="p", number=nrow(expData))
    topGenes <- topGenes[,colnames(topGenes) %in% c("logFC","t","P.Value","adj.P.Val")]
    colnames(topGenes) <- c("log2FC", "t-statistic", "p-value", paste("p-value", adjMethod,sep="_"))

    cat(paste("Writing differential expression analysis results for", contrastNames[i], "comparison to", paste(coreName, "_", colnames(eb$p.value)[i], "_topTable.txt", sep=""), "\n", sep=" "))
    
    write.table(prepare2write(topGenes), file=paste(coreName, "_", colnames(eb$p.value)[i], "_topTable.txt", sep=""),sep="\t", row.names=FALSE)
}

#####  ... write into a file only differnetially expressed genes
topGenes.index <- NULL

for (i in 1:ncol(eb$p.value) ) {
    
	topGenes <- topTable(eb, coef=colnames(eb)[i], adjust=adjMethod,  sort.by="p", p.value=pThreshold, lfc=lfcThreshold, number=nrow(expData))
    
    if ( length(topGenes) != 0 ) {
        
        topGenes <- topGenes[,colnames(topGenes) %in% c("logFC","t","P.Value","adj.P.Val")]
        colnames(topGenes) <- c("log2FC", "t-statistic", "p-value", paste("p-value", adjMethod,sep="_"))
        
        cat(paste("Writing differential expression analysis results for", contrastNames[i], "comparison to", paste(coreName, "_", colnames(eb$p.value)[i], "_DE_genes.txt", sep=""), "\n", sep=" "))
        
        write.table(prepare2write(topGenes), file=paste(coreName, "_", colnames(eb$p.value)[i], "_DE_genes.txt", sep=""),sep="\t", row.names=FALSE)
        
        topGenes.index <- c(topGenes.index,rownames(topGenes))
    } else {
        cat(paste("No differentially expressed genes detected in ", contrastNames[i], "\n", sep=""))
    }
}


#===============================================================================
#     Hierarchical clustering
#===============================================================================

#####  Select top discriminating genes
topGenes <- expData[unique(topGenes.index),]
dim(topGenes)

if (nrow(topGenes) > 0) {
    hr <- hclust(as.dist(1-cor(t(topGenes), method="pearson")), method="ward")
    hc <- hclust(as.dist(dist(t(topGenes), method="euclidean")), method="ward")

    pdf(paste(coreName,"_DE_dendrogram.pdf", sep = ""), width = 0.2*ncol(expData)+2, height = 6, pointsize = 12)
    par(mar=c(2,5,2,0))
    plot(hc, xlab="", labels=paste(colnames(expData), groups, sep="       "), hang = -1, main="")
    dev.off()

    pdf(paste(coreName,"_DE_heatmap.pdf", sep = ""), width = 6, height = 10, pointsize = 12)
    heatmap.2(as.matrix(topGenes), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=colorRampPalette(c("darkblue","darkblue","darkslateblue","darkslateblue","white","firebrick3","firebrick3","firebrick4","firebrick4"))(100), scale="row", ColSideColors=targets.colour[[2]] ,margins = c(2, 6), labRow=rownames(topGenes), labCol="", trace="none", key = FALSE)
    #####  Add the legend
    legend("topright", legend=groups,fill=targets.colour[[1]], box.col = "transparent", cex=0.6)
    dev.off()

    #####  Set.seed returns NULL, invisibly
    set.seed(1)

    pdf(paste(coreName,"_DE_cluster.pdf", sep = ""),width = 0.1*ncol(expData)+2, height = 6)
    hc$labels = colnames(expData)
    par(mar=c(2,2,2,6))
    A2Rplot(hc,
    k=groupsNo, # k = changes the detail/number of subgroups shown.
    fact.sup=annData[,2], box=FALSE,
    show.labels=TRUE, col.up = "black", main=" ")
    dev.off()

    #===============================================================================
    #     Pair-wise comparison
    #===============================================================================

    distance <- as.dist(dist(t(topGenes), method="euclidean"))
    plot.top.matrix <- as.matrix(distance)
    hc <- hclust(distance, method="ward")

    pdf(paste(coreName,"_DE_pairwise.pdf", sep = ""), width = 9, height = 8, pointsize = 12)
    heatmap.2(plot.top.matrix, symm = TRUE, distfun = as.dist, ColSideColors=targets.colour[[2]] , RowSideColors=targets.colour[[2]] , Rowv=as.dendrogram(hc), Colv=as.dendrogram(hc),margins = c(2, 6),labCol="",labRow="", trace="none", key = FALSE)
    #####  Add the legend
    legend("topright", legend=groups,fill=targets.colour[[1]], box.col = "transparent", cex=1.2)
    dev.off()
}

##### Clear workspace
rm(list=ls())
##### Close any open graphics devices
graphics.off()

