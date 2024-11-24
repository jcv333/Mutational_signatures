#################################################################################################.
#################################################################################################.
###
###   Author: Jose Camacho Valenzuela.
###   Date: XX.
###   R version: R version 4.1.3 (2022-03-10) -- "One Push-Up"
###   Script: Analyze mutational signatures with MutationalPatterns of a cancer patient.
###           - SBS and INDELs signatures - COSMIC catalogue.
###           - Input: 3 vcfs of 3 tumors of a breast cancer patient.
###
#################################################################################################.
#################################################################################################.

# 1. Load the required libraries.
library(MutationalPatterns)
library(BSgenome)
library(gridExtra)
library(data.table)
library(dplyr)
library(openxlsx)

# 2. Select the hg19 and import the input vcfs.
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
dir <- "C:/Users/path/to/dir/Input"
vcf_files <- list.files(path = dir, pattern = "*.vcf", full.names = TRUE)
list.files(dir)

# 3. Give the names to your input files.
list.files(dir)
name.files <- as.data.table(list.files(dir))
head(name.files$V1)
name.files$V1 <- gsub(".filtered.vcf", "", name.files$V1)
head(name.files$V1)
sample_names <- name.files$V1

# 4. GRangesList.
grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
summary(grl)
tissue <- c(rep("Breast", 3))

# 5. Get the matrix.
muts <- mutations_from_vcf(grl[[1]])
types <- mut_type(grl[[1]])
context <- mut_context(grl[[1]], ref_genome)
type_context <- type_context(grl[[1]], ref_genome)
lapply(type_context, head, 12)
type_occurrences <- mut_type_occurrences(grl, ref_genome)
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)

# 6. Plot the mutational profile.
setwd("C:/Users/path/to/dir/Output")
plot_96_profile(mut_mat[, c(1:3)], condensed = TRUE, ymax = 0.2)


# --------------------------------------------------------------------------------------------- #
# -------------------             Single Base Substitution (SBS)            ------------------- #
# --------------------------------------------------------------------------------------------- #

# 7. Use the most updated COSMIC catalog.
signatures = get_known_signatures(muttype = c("snv"), source = c("COSMIC_v3.2"))
setwd("C:/Users/path/to/dir/Output")

# 8. Strict delta value.
{strict_refit_sbs <- fit_to_signatures_strict(mut_mat, signatures, max_delta = 0.04)
  head(strict_refit_sbs)
  
  fit_res_strict_sbs <- strict_refit_sbs$fit_res
}

# 9. Plotting. Manual save of the plots.
# ----- BARPLOT. Absolute contribution.
plot_contribution(fit_res_strict_sbs$contribution,
                  coord_flip = F,
                  mode = "absolute")

# ----- BARPLOT. Relative contribution.
plot_contribution(fit_res_strict_sbs$contribution,
                  coord_flip = F,
                  mode = "relative")

# ----- HEATMAP. Relative contribution.
plot_contribution_heatmap(fit_res_strict_sbs$contribution,
                          sample_order = NA,
                          cluster_samples = FALSE,
                          method = "complete",
                          plot_values = FALSE)


# --------------------------------------------------------------------------------------------- #
# -------------------              INDELs signatures (SBS)                  ------------------- #
# --------------------------------------------------------------------------------------------- #

# 10. Rename the sample names for ID signatures.
sample_names_indels <- sample_names

# 11. Get the matrix.
indel_grl <- read_vcfs_as_granges(vcf_files, sample_names_indels, ref_genome, type = "indel", predefined_dbs_mbs = TRUE)
indel_grl <- get_indel_context(indel_grl, ref_genome)
indel_counts <- count_indel_contexts(indel_grl)
head(indel_counts)

# 12. Plot the INDEL spectra.
plot_main_indel_contexts(indel_counts)
plot_indel_contexts(indel_counts, condensed = TRUE)

# 13. Use the most updated COSMIC catalog.
signatures_indels = get_known_signatures(muttype = c("indel"), source = c("COSMIC_v3.2"))
setwd("C:/Users/path/to/dir/Output")

# 14. Strict delta value.
{strict_refit_IDs <- fit_to_signatures_strict(indel_counts, signatures_indels, max_delta = 0.004)
  head(strict_refit_IDs)

  fit_res_strict_IDs <- strict_refit_IDs$fit_res
}

# 15. Plotting. Manual save of the plots.
# ----- BARPLOT. Absolute contribution.
plot_contribution(fit_res_strict_IDs$contribution,
                  coord_flip = FALSE,
                  mode = "absolute")

# ----- BARPLOT. Relative contribution.
plot_contribution(fit_res_strict_IDs$contribution,
                  coord_flip = FALSE,
                  mode = "relative")

# ----- HEATMAP. Relative contribution.
plot_contribution_heatmap(fit_res_strict_IDs$contribution,
                          sample_order = NA,
                          cluster_samples = FALSE,
                          method = "complete",
                          plot_values = FALSE)

#############################################################################################.
#############################################################################################.
##########                                                                        ###########.
##########                                 The end.                               ###########.
##########                                                                        ###########.
#############################################################################################.
#############################################################################################.
