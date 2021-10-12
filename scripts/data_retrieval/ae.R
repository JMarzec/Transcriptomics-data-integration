################################################################################
#
#   File name: ae.R
#
#   Authors: Jacek Marzec ( j.marzec@qmul.ac.uk ), code modified from Rossalind Cutts
#
#   Barts Cancer Institute,
#   Queen Mary, University of London
#   Charterhouse Square, London EC1M 6BQ
#
################################################################################

################################################################################
#
#	Description: Getting data from ArrayExpress using ArrayExpress package
#
#	Command line use: R --file=./ae.R --args "E-MEXP-993" "raw|processed|full" "/scratch/jack/data/PhD/raw/Affy_U133Plus2" - this will get the raw data files if they exist
#	
#	First arg: ArrayExpress id
#	Second arg: whether data is raw or matrix format
#	Third arg: the directory for the data
#	
################################################################################

#===============================================================================
#    Main
#===============================================================================

library("ArrayExpress")

Args <- commandArgs();

ArrayExpressID=Args[4]
FileType=Args[5]
DataDir=paste(Args[6], ArrayExpressID, sep="/")

print(ArrayExpressID)
print(FileType)
print(DataDir)

##### Set/create a directory for data download
if (file.exists(DataDir)){
	setwd(DataDir)
} else {
	dir.create(DataDir);
	setwd(DataDir)
}

##### Get experiment info and data
rawset = ArrayExpress(ArrayExpressID, path = DataDir, save = FALSE)

##### Get raw|processed or both (full) files from ArrayExpress
files=getAE(ArrayExpressID,type=FileType,path=DataDir)

phenodata = NULL

##### Get experiment phenodata
if (is.list(rawset)) {
    
    samples = NULL
    
    for (i in 1:length(rawset)) {
        
        phenodata = rbind(phenodata,pData(phenoData(rawset[[i]])))
        samples = c(samples, sampleNames(rawset[[i]]))
    }
} else {
    phenodata=pData(phenoData(rawset))
    samples = sampleNames(rawset)
}

targetFile=NULL

##### Generate target file
for (i in 1:length(names(phenodata))) {
	
	if ( length(grep("Source.Name",names(phenodata)[i])) != 0 ) {
		
		targetFile=data.frame(phenodata[i])
						  
	} else if ( length(grep("Characteristics",names(phenodata)[i])) != 0 ) {
		
		targetFile=data.frame(targetFile, phenodata[i])
		
	} else if ( length(grep("Factor",names(phenodata)[i])) != 0 ) {
		
		targetFile=data.frame(targetFile, phenodata[i])
	}
}

targetFile=data.frame(samples,files$rawFiles,targetFile)

names(targetFile)=gsub("[.]{2}","_",sapply(names(targetFile), function(elt) elt[length(elt)]))
names(targetFile)=gsub("[.]","",sapply(names(targetFile), function(elt) elt[length(elt)]))
names(targetFile)[1:2]=c("Name","FileName")

write.table(targetFile,file="target.txt",sep="\t",row.names=FALSE,quote=FALSE)

##### Remove the compressed raw|processed data files
system("rm *.zip")	

