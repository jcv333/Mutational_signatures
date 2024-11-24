#########################################################################################################.
###
###   Author: Jose Camacho Valenzuela.
###   Date: XX.
###   R version: R version 4.1.3 (2022-03-10) -- "One Push-Up"
###   Script: Analyze Signature 3 using SigMA.
###             - Input files: 3 vcfs of 3 individual patient's tumors from exome-sequencing.
###
#########################################################################################################.

# 1. Load the package SigMA.
library(SigMA)

# 2. Set the correct working directory.
setwd("C:/Users/path/to/working/dir/SigMA")

# 3. Load devtools.
devtools::load_all()

# 4. Import the tumor vcfs.
data_dir <- "C:/Users/path/to/dir/Input"
list.files(data_dir)

# 5. Make matrix.
genomes_matrix <- make_matrix(data_dir, file_type = 'vcf')
genomes <- conv_snv_matrix_to_df(genomes_matrix)

# 6. Output directory and SigMA analysis.
genome_file = "C:/Users/path/to/dir/Output/patient1.csv"
write.table(genomes,
            genome_file,
            sep = ',',
            row.names = F,
            col.names = T ,
            quote = F)

message(paste0('96-dimensional matrix is saved in ', genome_file))
message('Running SigMA')

run(genome_file, 
    data = "seqcap",
    do_assign = T, 
    do_mva = T,
    tumor_type = 'breast',
    lite_format = T, add_sig3 = T,
    check_msi = T)


# ------------------------------------------------------------------------------------------- #
# ---------------               SigMA output customization                    --------------- #
# ------------------------------------------------------------------------------------------- #

# 7. Load the required libraries.
library(data.table)
library(openxlsx)
library(tidyverse)
library(openxlsx)
library(readr)
library(readxl)

# 8. Read the SigMA output file.
setwd("C:/Users/path/to/dir/Output")
list.files()
SigMA.output <- fread("patient1_output_tumortype_breast_platform_WholeExomeSequencing37Mb_cf1_lite.csv") %>% setDT()

# 9. Modify both the header and the name of the sample's ID.
colnames(SigMA.output)
setnames(SigMA.output, "tumor", "Tumor")
head(SigMA.output$Tumor)
SigMA.output$Tumor <- gsub("C:/Users/path/to/dir/Input/patient1_tumor1.vcf",
                           "Tumor1",
                           SigMA.output$Tumor)
head(SigMA.output$Tumor)

SigMA.output$Tumor <- gsub("C:/Users/path/to/dir/Input/patient1_tumor2.vcf",
                           "Tumor2",
                           SigMA.output$Tumor)
head(SigMA.output$Tumor)

SigMA.output$Tumor <- gsub("C:/Users/path/to/dir/Input/patient1_tumor3.vcf",
                           "Tumor3",
                           SigMA.output$Tumor)
head(SigMA.output$Tumor)


# 10. Write the resulting dataframe as xlsx output file.
setwd("C:/Users/path/to/dir/Output")
openxlsx::write.xlsx(SigMA.output,
                     file = "patient1_SigMA_customized.xlsx",
                     sheetName = "Sheet1",
                     colNames = T, rowNames = F, append = F)


# ########################################################################################### #.
# ########################################################################################### #.
# #########                                                                        ########## #.
# #########                                 The end.                               ########## #.
# #########                                                                        ########## #.
# ########################################################################################### #.
# ########################################################################################### #.
