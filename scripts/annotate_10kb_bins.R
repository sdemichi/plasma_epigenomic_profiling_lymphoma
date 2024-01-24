# Goal of this script is to annotate coding (including gene names) and non-coding regulatory features in all 10kb, genome wide adjacent bins across the human genome.
# This reference file will be used for project analyses, related to annotating features from machine learning models, and other computational tasks.
# This script can be used to annotate different bin sizes.

# libraries
library(annotatr)

# import 10kb reference. Reference was created prior.
setwd("path/to/references/")
ref_10kb <- read.delim("10kb_bin_genome_wide_genomic_ranges_reference_raw.txt", header = TRUE, sep = "\t")
# ref_100kb <- read.delim("100kb_bin_genome_wide_genomic_ranges_reference_raw.txt", header = TRUE, sep = "\t")

# annotate using annotatR. NOTE; annotatR requires an internet connection to access TxDb and other annotation databases. 
# subset to only chrom, start, end
bin_annotations_regions <- ref_10kb[,1:3]
# bin_annotations_regions <- ref_100kb[,1:3]

# rename columns
bin_annotations_regions <- setNames(bin_annotations_regions, c("seqnames","start","end"))

# create GRanges object from the 10kb reference
bin_annotations_regions_gr <- GRanges(bin_annotations_regions)

##
# Create three separate annotation sets as character vectors; CpG annotations, noncoding regulatory, coding regulatory.
annotations_hg38_basicgenes <- c("hg38_basicgenes")
annotations_hg38_genic <- c("hg38_genes_1to5kb",
                            "hg38_genes_promoters", 
                            "hg38_genes_5UTRs", 
                            "hg38_genes_exons",
                            "hg38_genes_firstexons",
                            "hg38_genes_introns",
                            "hg38_genes_intronexonboundaries",
                            "hg38_genes_exonintronboundaries",
                            "hg38_genes_3UTRs",
                            "hg38_genes_intergenic")
annotations_hg38_enhancers <- c("hg38_enhancers_fantom")

# build annotations for each category
annotation_build_hg38_cpg <- build_annotations(genome = "hg38", annotations = "hg38_cpgs") # basic annotation preset; CpG features
annotation_build_hg38_genic <- build_annotations(genome = "hg38", annotations = annotations_hg38_genic) # genic regulatory features
annotation_build_hg38_basicgenes <- build_annotations(genome = "hg38", annotations = annotations_hg38_basicgenes) # gene names
annotation_build_hg38_enhancers <- build_annotations(genome = "hg38", annotations = annotations_hg38_enhancers) # enhancers (FANTOM5 w liftover)


##
# annotate regions, for each of the desired annotation sets

# CpG annotations
annotate_300bp_bins_hg38_cpg <- annotate_regions(regions = bin_annotations_regions_gr,
                                                 annotations = annotation_build_hg38_cpg,
                                                 quiet = FALSE)
unique_annotate_300bp_bins_hg38_cpg <- unique(annotate_300bp_bins_hg38_cpg) # remove duplicate annotations

# gene regulatory element annotations
annotate_300bp_bins_hg38_genic <- annotate_regions(regions = bin_annotations_regions_gr,
                                                 annotations = annotation_build_hg38_genic,
                                                 quiet = FALSE)
unique_annotate_300bp_bins_hg38_genic <- unique(annotate_300bp_bins_hg38_genic)

# enhancer annotations
annotate_300bp_bins_hg38_enhancers <- annotate_regions(regions = bin_annotations_regions_gr,
                                                 annotations = annotation_build_hg38_enhancers,
                                                 quiet = FALSE)
unique_annotate_300bp_bins_hg38_enhancers <- unique(annotate_300bp_bins_hg38_enhancers)


##
# gene annotations (requires additional code, due to duplicate genes in bins)
annotate_300bp_bins_hg38_basicgenes <- annotate_regions(regions = bin_annotations_regions_gr,
                                                 annotations = annotation_build_hg38_basicgenes,
                                                 quiet = FALSE)

# convert to data frame
annotate_300bp_bins_hg38_basicgenes <- as.data.frame(annotate_300bp_bins_hg38_basicgenes) # the 10kb reference now has bin IDs from start

# remove rows where annot.symbol == "NA"
annotate_300bp_bins_hg38_basicgenes <- annotate_300bp_bins_hg38_basicgenes[!(is.na(annotate_300bp_bins_hg38_basicgenes$annot.symbol)),]

# use dplyr for this, group by bin_ID, and remove duplicates of annot.symbol
gene_names_per_bin <- annotate_300bp_bins_hg38_basicgenes %>% 
  group_by(bin_id) %>% 
  summarize(gene_names = paste(unique(annot.symbol), collapse = ", "))

# convert to data frame
gene_names_per_bin <- as.data.frame(gene_names_per_bin)


##
# Convert the GRanges objects back to data frames
df_unique_annotate_300bp_bins_hg38_cpg <- as.data.frame(unique_annotate_300bp_bins_hg38_cpg)
df_unique_annotate_300bp_bins_hg38_genic <- as.data.frame(unique_annotate_300bp_bins_hg38_genic)
df_unique_annotate_300bp_bins_hg38_enhancers <- as.data.frame(unique_annotate_300bp_bins_hg38_enhancers)

# Match the annotated data frame to the 10kb reference, for all references
df_merged_cpg <- merge(
  x = bin_annotations_regions, # the left data frame to be merged
  y = df_unique_annotate_300bp_bins_hg38_cpg, # The right data frame to be merged
  by = c("seqnames", "start", "end"), # columns for merging
  all = TRUE # include all rows from both data frames
)

df_merged_genic <- merge(
  x = bin_annotations_regions, 
  y = df_unique_annotate_300bp_bins_hg38_genic, 
  by = c("seqnames", "start", "end"), 
  all = TRUE
)

df_merged_enhancers <- merge(
  x = bin_annotations_regions, 
  y = df_unique_annotate_300bp_bins_hg38_enhancers, 
  by = c("seqnames", "start", "end"), 
  all = TRUE
)

# Bind new annotated features to original data frame of annotated bins.
bin_annotations_final <- cbind(ref_10kb,df_merged_cpg$annot.type,df_merged_genic$annot.type,df_merged_enhancers$annot.type)
# bin_annotations_final <- cbind(ref_100kb,df_merged_cpg$annot.type,df_merged_genic$annot.type,df_merged_enhancers$annot.type) # for 100kb reference

# merge bin_annotations_final and gene_names_per_bin, as there are many bins missing from gene_names_per_bin, since not all bins have genes.
bin_annotations_final <- merge(bin_annotations_final, gene_names_per_bin, by = "bin_id", all.x = TRUE)

# make all NAs in the data frame "other"
bin_annotations_final <- data.frame(lapply(bin_annotations_final, function(x) ifelse(is.na(x), "other", x)))

# set column names of annotated reference                                           
bin_annotations_final <- setNames(bin_annotations_final,
                                  c("chrom","start","end","width","strand",
                                    "bin_id","cpg_annotations","genic_annotations","enhancer_annotations", "gene_name"))

##
#export annotation set at text file
setwd("path/to/references/")

# write annotated reference to text file                                           
write.table(bin_annotations_final,
            file = "10kb_bin_genome_wide_annotated_reference_with_regulatory_features_and_gene_names.txt",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
