####################
#MEDIPS - create set
####################
library(tidyverse)
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("path/to/bam/files/")
sample_bam <- list.files(path = ".", pattern = "dedup.bam$", recursive = TRUE) # list BAM files from MEDIPIPE.

chr <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
bsgenome <- "BSgenome.Hsapiens.UCSC.hg38"

# create MEDIPs sets from all BAM files in a directory
for (i in sample_bam){
  temp <- MEDIPS.createSet(file = i,
                   BSgenome = bsgenome,
                   paired = TRUE, # for paired end files only
                   # paired = FALSE, # for single end files only
                   window_size = 10000, # change, depending on desired window size.
                   uniq = 0,
                   chr.select = chr)
  setwd("path/to/counts/")
  MEDIPS.exportWIG(Set = temp, file = paste0(substr(i, 1, nchar(i)-4),"_counts_10kb.txt"), format = "count", descr = i)
  setwd("path/to/bam/files/")
  }
