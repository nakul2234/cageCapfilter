# cageCapfilter

The python script will take a CAGE BAM file and convert it into a CTSS file. It makes two files:

1. *CTSS <- This is a traditional CTSS file that can be used by CAGEr
2. *_unannotatedG.CTSS <- This is a CTSS with only the reads with unannotated G's with more confidence of actually being the start of a read. 

There is an R script that shows how these can be imported into R with the CAGEr package to be able to filter the peaks to get the ones that only have a certain thresshold of unannotated G's.
