####################
#MEDIPS - create set
####################
library(tidyverse)
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("~/Library/CloudStorage/OneDrive-UniversityofToronto/1 MBP/1 Bratman Lab/1 Experiments/11 Sequencing/19-08-2022_k4_k27me3_troubleshooting/31-10-2023 letter to the editor/bam")
sample_bam <- list.files(path = ".", pattern = "dedup.bam$", recursive = TRUE) #list BAM files from MEDIPIPE.

chr <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
bsgenome <- "BSgenome.Hsapiens.UCSC.hg38"

## MEDIPs set
for (i in sample_bam){
  temp <- MEDIPS.createSet(file = i,
                   BSgenome = bsgenome,
                   paired = TRUE, #for paired end files only
                   #paired = FALSE,
                   window_size = 300,# change, depending on desired window size.
                   uniq = 0,
                   chr.select = chr)
  setwd("~/Library/CloudStorage/OneDrive-UniversityofToronto/1 MBP/1 Bratman Lab/1 Experiments/11 Sequencing/19-08-2022_k4_k27me3_troubleshooting/31-10-2023 letter to the editor/counts")
  MEDIPS.exportWIG(Set = temp, file = paste0(substr(i, 1, nchar(i)-4),"_counts_300bp.txt"), format = "count", descr = i)
  setwd("~/Library/CloudStorage/OneDrive-UniversityofToronto/1 MBP/1 Bratman Lab/1 Experiments/11 Sequencing/19-08-2022_k4_k27me3_troubleshooting/31-10-2023 letter to the editor/bam")
  }

#################Continue working on this; export all non-normalized counts extracted from short BAM files. Change working directories back and forth or
###set variables for "import" and "export".

##Once complete, create new DESeq2 script to recognize short fragment counts files, and create a new sample list (using chunk 7 of the .Rmd)
