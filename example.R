#Ran this on R3.4.1 and that was the only way to get the code to work to cluster. I then will reload the input to get the final result. 

library("ggplot2")
library('CAGEr')
library('BSgenome.Hsapiens.UCSC.hg38')


# Example script that will use the CTSS files output from the python script to do the Capfilter. 
# Importing the CTSS files. These include annotated and unannotated G reads. 

ifiles=list.files(pattern="m\\.CTSS$")

cs = new("CAGEset", genomeName = "BSgenome.Hsapiens.UCSC.hg38", inputFiles = ifiles, inputFilesType = "ctss", sampleLabels = ifiles)

getCTSS(cs)
normalizeTagCount(cs, method="simpleTpm")

clusterCTSS(object = cs, method = "distclu", threshold = 0.1, nrPassThreshold = 2, thresholdIsTpm = TRUE, removeSingletons = TRUE, keepSingletonsAbove = 0.2, minStability = 2, maxLength = 100, reduceToNonoverlapping = TRUE, useMulticore = TRUE, nrCores = 12)

aggregateTagClusters(cs, tpmThreshold = 0.3, qLow = NULL, qUp = NULL, maxDist = 100, excludeSignalBelowThreshold = FALSE)


cc.df <- consensusClusters(cs, sample = NULL, returnInterquantileWidth = FALSE, qLow = NULL, qUp = NULL)
cc.df.tpm <- consensusClustersTpm(cs)

ctss.tags <- CTSStagCount(cs)

library("snow")

clus <- makeCluster(12)

clusterExport(clus,"ctss.tags")

resultapply <- parApply(clus,cc.df[,c("chr", "start", "end", "strand")],1,function(x) colSums(ctss.tags[ctss.tags$chr == x[1] & ctss.tags$pos >= as.numeric(x[2]) & ctss.tags$pos <= as.numeric(x[3]) & ctss.tags$strand == x[4], 4:ncol(ctss.tags)]))
stopCluster(clus)

# Importing the CTSS files for the unannotated G files

ifiles2=list.files(pattern="G\\.CTSS$")

cs2 = new("CAGEset", genomeName = "BSgenome.Hsapiens.UCSC.hg38", inputFiles = ifiles2, inputFilesType = "ctss", sampleLabels = ifiles2)

getCTSS(cs2)

ctss.tags.2 <- CTSStagCount(cs2)

library("snow")

clus <- makeCluster(12)

clusterExport(clus,"ctss.tags.2")

resultapply2 <- parApply(clus,cc.df[,c("chr", "start", "end", "strand")],1,function(x) colSums(ctss.tags.2[ctss.tags.2$chr == x[1] & ctss.tags.2$pos >= as.numeric(x[2]) & ctss.tags.2$pos <= as.numeric(x[3]) & ctss.tags.2$strand == x[4], 4:ncol(ctss.tags.2)]))
stopCluster(clus)i

allCTSSpeakcount <- as.data.frame(t(resultapply))
uGCTSSpeakcount <- as.data.frame(t(resultapply2))

uGPercent <- uGCTSSpeakcount/allCTSSpeakcount

uGPercent[is.na(uGPercent)] <- 0

cc.df$maxUGPercentage <- apply(uGPercent, 1, function(x) max(x))

cc.df.tpm <- as.data.frame(cc.df.tpm)

cc.df <- cbind(cc.df, cc.df.tpm)

# This will filter peaks that have an unannotated g percentage o greater than .15

threshold <- .15

cc.df.fil <- cc.df[cc.df$maxUGPercentage >= threshold,]
