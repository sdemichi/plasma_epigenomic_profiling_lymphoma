# Goal; to generate a 10 kilobase (10kb) bin reference, of genome wide, adjacent 10kb bins.
# Reference will be extracted from a MEDIPS object, generated from a sample BAM file.

# libraries
library(tidyverse)
library(GenomicRanges)
library(IRanges)

# set working directory for listing BAM files
setwd("path/to/bam/files/")

# list one of the bam files
sample_bam <- list.files(path = ".", pattern = "sample_name.bam$") # list BAM file

# save chromosome attributes and genome name as variables
chr <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
bsgenome <- "BSgenome.Hsapiens.UCSC.hg38"

# create medips set
temp <- MEDIPS.createSet(file = sample_bam,
                   BSgenome = bsgenome,
                   paired = TRUE, #for paired end files only
                   window_size = 10000, # change, depending on desired window size.
                   uniq = 0,
                   chr.select = chr)

# extract GRanges object from MEDIPS object.
chr.select = temp@chr_names
window_size = window_size(temp)
chr_lengths = unname( seqlengths(BSgenome.Hsapiens.UCSC.hg38)[ seqnames(BSgenome.Hsapiens.UCSC.hg38@seqinfo)%in%chr.select ] )
no_chr_windows = ceiling(chr_lengths/window_size)
supersize_chr = cumsum(no_chr_windows)
chromosomes = chr.select

all.Granges.genomeVec = MEDIPS.GenomicCoordinates(supersize_chr, no_chr_windows, chromosomes, chr_lengths, window_size)

# convert to data frame
final_granges_10kb <- as.data.frame(all.Granges.genomeVec)

# mutate row names as bin_id
final_granges_10kb <- final_granges_10kb %>% 
  mutate(bin_id = rownames(final_granges_10kb))

# export as txt file
setwd("/cluster/projects/scottgroup/people/steven/reference_sets/other_references")

# write to text file. This file will be used as input for regulatory element annotation, using annotatR.
write.table(final_granges_10kb,
            file = "10kb_bin_genome_wide_genomic_ranges_reference_raw.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
