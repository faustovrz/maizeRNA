# Kallisto RNA-seq Statistics and Analysis
# This script processes Kallisto quantification results

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# Read sample metadata 
sample_info <- read.csv("../data/sample_info.csv", header = TRUE)
plots <- read.csv("../data/plots_info.csv")

# Join metadata
metadata <- sample_info %>%
  dplyr::inner_join(plots, by= c(row ="PHO22")) 

# Process time information
metadata <- within(metadata, {
  TIME = sub("1:","13:",metadata$TIME)
  TIME = hm(TIME)
  decimal_time <- hour(TIME) + minute(TIME)/60 + second(TIME) / 3600
})

# Save processed metadata
write.csv(
  metadata, 
  file="../data/inv4mRNAseq_metadata.csv", 
  row.names = FALSE
)

# Read kallisto QC data
kalqc <- read.table("../results/kallisto_qc.tab", sep = "\t", 
                   header = FALSE, skip = 1)
colnames(kalqc) <- c("gsl_sample","pseudoaligned_pct",
                    "pseudoaligned", "processed")
kalqc$tube <- substr(kalqc$gsl_sample, 1, 3)

# Plot QC statistics
kalqc %>%
  tidyr::pivot_longer(pseudoaligned_pct:processed, 
                     names_to = "var", values_to = "value") %>%
  ggplot2::ggplot(aes(x=value)) +
  ggtitle("Kallisto pseudoalignment") +
  ylab("# tissue libraries") +
  xlab("read stat") +
  geom_histogram() +
  facet_wrap(~ var, scales = "free") +
  theme_classic()

# Save the plot
ggsave("../results/kallisto_alignment_stats.png", width = 8, height = 6)

# Get sample directories
samples <- dir("../results/quant_out")

# Process gene expression data
all_exp <- lapply(samples, function(x) {
  print(x)
  sample_exp <- read.table(
    file.path("../results/quant_out", x, "abundance.tsv"), 
    sep = "\t", header = TRUE
  )
  sample_exp$est_counts <- as.integer(sample_exp$est_counts)
  sample_exp$gene = factor(sub("_T.+$", "", sample_exp$target_id, perl = TRUE))

  out <- sample_exp %>%
    dplyr::group_by(gene) %>%
    dplyr::summarize(counts = sum(est_counts))
  out$sample <- x
  out
}) %>% dplyr::bind_rows()

# Create gene-sample expression matrix
gene_sample_exp <- all_exp %>% 
  pivot_wider(names_from = "sample", values_from = "counts")

# Save expression matrix
write.csv(gene_sample_exp, "../results/gene_sample_exp.csv")
