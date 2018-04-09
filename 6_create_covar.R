# -------------------------------------------------------
# Making the covariate file
# -------------------------------------------------------

# Covariatess:
# - genetic PCs 1:10
# - cell counts
# - age

# Only using quantitative covariates

# -------------------------------------------------------
# Load the data needed
# -------------------------------------------------------
library(tidyverse)
library(stringr)

setwd("/panfs/bt14")

# ARIES samplesheet
samplesheet <- read.delim('epi_data/ARIES_samplesheet.txt', stringsAsFactors = FALSE)

# Cell counts
cell_counts <- read.delim('epi_data/ARIES_cellcounts.txt', stringsAsFactors = FALSE)

fam <- read.table("kinship/methyl_FOM.fam", header = F)

# -------------------------------------------------------
# Extract the covariates needed from the data
# -------------------------------------------------------

#Principal components
PC_dat <- read.table("pca/FOM_pcs.eigenvec", sep = " ", header = F, stringsAsFactors = F)
head(PC_dat)
colnames(PC_dat) <- c("FID", "IID", paste0(rep("PC", times = 20), 1:20))
PC_dat$cidB3030 <- gsub("[A-Z]", "", PC_dat[["FID"]])
PC_dat <- dplyr::select(PC_dat, -IID, -FID)

#Ensure column names correct
colnames(samplesheet)[colnames(samplesheet) == "Sample_Name"] <- "IID"
colnames(fam) <- c("FID", "IID", "PID", "MID", "SEX", "PHENO")
fam$FID <- as.character(fam$FID)

# Select PCs 1:10 
PCs <- left_join(fam, PC_dat, by = c("FID" = "cidB3030"))
PCs <- dplyr::select(PCs, -PID, -MID, -SEX, -PHENO)
PCs <- PCs[, 1:12]

nrow(PCs) - nrow(PCs[complete.cases(PCs), ])
# 109 people removed due to lack of genotype info

	
# Combine the variables into one table
dat <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "IID"), list(fam, samplesheet, cell_counts, PCs))

head(dat)
dat <- dat[!is.na(dat$PC1), ]
dim(dat)

# Sanity check
stopifnot(sum(dat$FID.x == dat$cidB3030) == nrow(dat))
stopifnot(sum(dat$FID.y == dat$cidB3030) == nrow(dat))

# -------------------------------------------------------
# Covariates for methylation GCTA
# -------------------------------------------------------

# Select the covariates for the analysis
m_covars <- c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "age", paste0("PC", 1:10))

# Extract values for individuals
meth_covars <- dplyr::select(dat, cidB3030, IID, one_of(m_covars))

head(meth_covars)

write.table(meth_covars, "phen/FOM.qcovar", quote = F, row.names = F, col.names = F)

# -------------------------------------------------------
# Covariates for genetic GCTA
# -------------------------------------------------------

g_covars <- c("age", paste0("PC", 1:10))
gen_covars <- dplyr::select(dat, cidB3030, cidB, one_of(g_covars))
head(gen_covars)


write.table(gen_covars, paste0("phen/FOM_gen.qcovar"), quote = F, row.names = F, col.names = F)
