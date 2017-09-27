################################################################################
#
#   File name: MultiGene2ProbeFilter.R
#
#   Authors: Jacek Marzec ( j.marzec@qmul.ac.uk ),  code adopted from 'screenExpr.R' by Syed Heider
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

################################################################################
#
#	Description: Removes multi gene to probe/probeset mappings by selecting probe/probeset with highest standard deviation across samples. The first column of the input expression matrix is expected to contain names of genes mapping to corresponding probes/probesets.
#
################################################################################

multiGene2ProbeFilter <- function(exprsMatrix) {
    
	cat(paste("Processing expression data for ",dim(exprsMatrix)[1]," genes.\n",sep=""))
    
	#####  Remove probes that do not map to any gene
	exprs.mapped <- exprsMatrix[toupper(exprsMatrix[,1]) != "UNMAPPED" & !is.na(exprsMatrix[,1]),]
	cat(paste("Probes mapped to ",dim(exprs.mapped)[1]," genes names.\n",sep=""))
    
	#####  Compute sample variance
	exprs <- apply(exprs.mapped[,-1], 2,as.numeric)
	rownames(exprs) <- rownames(exprs.mapped)
	nLastCol <- ncol(exprs)
	sampleVar <- apply(exprs,1,var)
    
	nGenes <- nrow(exprs)
    
	genesUnique <- as.vector(unique(exprs.mapped[,1]))
    
	#####  [Syed: to avoid NA issues]
	genesUnique <- as.vector(na.omit(genesUnique))
	nGenesUnique <- length(genesUnique)
    
	nSamples <- ncol(exprs)
    
	#####  Create matrix for filtered data
	exprs.filtered <- array(dim=c(nGenesUnique,nLastCol))
	rownames(exprs.filtered) <- genesUnique
	colnames(exprs.filtered) <- colnames(exprs)
    
	cat("Removing duplicate genes (selecting for highest standard deviation)...\n")
    
	for (gene in genesUnique) {
        
            #####  Index/indices of all genes matching gene (detect duplicates)
            idxGenes <- seq(along=1:nGenes)[exprs.mapped[,1]==gene];
        
            exprs.slice <- exprs[idxGenes,];
        
            #####  Find duplicates with highest variance
            if (length(idxGenes)>1) {
                idxMaxVar <- which.max(sampleVar[idxGenes]);
            
                exprs.slice <- exprs.slice[idxMaxVar,];
            }
            exprs.filtered[rownames(exprs.filtered)==gene,] <- as.matrix(exprs.slice);
	}
    
	cat(paste(nGenesUnique," unique genes names remain.\n",sep=""))
    
	rm(exprsMatrix, exprs, exprs.mapped, exprs.slice)
    
	return(exprs.filtered)
}
