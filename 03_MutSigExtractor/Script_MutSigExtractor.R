#########################################################################################.
#########################################################################################.
#####.
#####   Author: Jose Camacho Valenzuela.
#####   Date: XX.
#####   R version: R version 4.1.3 (2022-03-10) -- "One Push-Up".
#####   Script: Analyze mutational signatures with MutSigExtractor of a cancer patient.
#####           - Analyze SBS and INDELs signatures - COSMIC catalogue.
#####           - Input: 3 vcfs of 3 tumors of a breast cancer patient.
#####.
#########################################################################################.
#########################################################################################.

# 1. Load the required libraries.
library(mutSigExtractor)
library(data.table)
library(openxlsx)
library(tidyverse)

# 2. Load the input vcf.
setwd("C:/Users/path/to/file")
list.files()
vcf_snv <- 'patient1.vcf'


# ----------------------------------------------------------------------------------------------.
################################################################################################.
##### ------------------------------------------------------------------------------------ #####.
#####                          Determine SBS signatures                                    #####.
##### ------------------------------------------------------------------------------------ #####.
################################################################################################.

# 3. Extract the mutational profile of the sample in question.
contexts_snv <- extractSigsSnv(vcf.file = vcf_snv,
                               vcf.filter = 'PASS',
                               output = 'contexts')
head(contexts_snv)

# 4. Perform signature fitting.
sigs_snv_PCAWG <- fitToSignatures(mut.context.counts = contexts_snv[,1],
                                  signature.profiles = SBS_SIGNATURE_PROFILES_V3)

head(sigs_snv_PCAWG)
a <- as.data.table(sigs_snv_PCAWG, keep.rownames = T)


################################################################################################.
##### ------------------------------------------------------------------------------------ #####.
#####                          Determine INDEL signatures                                  #####.
##### ------------------------------------------------------------------------------------ #####.
################################################################################################.

# 5. Extract the mutational indel profile / contexts of the sample in question.
contexts_indel <- extractSigsIndel(vcf.file = vcf_snv,
                                   vcf.filter = 'PASS',
                                   method = 'PCAWG',
                                   output = 'contexts')
head(contexts_indel)

# 6. Perform signature fitting with the PCAWG INDELs signatures (from COSMIC).
sigs_indel_PCAWG <- fitToSignatures(mut.context.counts = contexts_indel[,1],
                                    signature.profiles = INDEL_SIGNATURE_PROFILES)

head(sigs_indel_PCAWG)
b <- as.data.table(sigs_indel_PCAWG, keep.rownames = T)

# 7. Rename properly the headers of the resulting objects.
colnames(a)
setnames(a, "rn", "SBS_PCAWG_Signatures")
setnames(a, "sigs_snv_PCAWG", "Exposure")

colnames(b)
setnames(b, "rn", "INDEL_PCAWG_Signatures")
setnames(b, "sigs_indel_PCAWG", "Exposure")

# 8. Write the final ouput dataframes as a single output xlsx.
sigs.outputs <- list("SBS_COSMIC_PCAWG_60" = a,
                     "IDS_COSMIC_PCAWG" = b)

setwd("C:/Users/path/to/dir/MutSigExtractor")
write.table(a,
            file = "patient1_MutSigExtractor_SBS_signatures.txt",
            sep = "\t",
            col.name = T,
            row.names = F,
            quote = F)

write.table(b,
            file = "patient1_MutSigExtractor_ID_signatures.txt",
            sep = "\t",
            col.name = T,
            row.names = F,
            quote = F)

openxlsx::write.xlsx(sigs.outputs,
                     "patient1_MutSigExtractor_SBS_ID_signatures.xlsx",
                     colNames = T,
                     rowNames = F,
                     append = F)


# ###########################################################################################.
# ###########################################################################################.
# #########                                                                        ##########.
# #########                                 The end.                               ##########.
# #########                                                                        ##########.
# ###########################################################################################.
# ###########################################################################################.

