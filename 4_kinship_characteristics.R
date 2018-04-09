rm(list = ls())
library(tidyverse)
setwd("/panfs/bt14")

source("GRM_functions.R")

# Read in the methylation/genetic kinship matrix
GRM_dat <- readGRM(rootname = "kinship/methyl_FOM")

str(GRM_dat)

# -------------------------------------------------------
# Examine the data, create histograms
# -------------------------------------------------------

off_diag <- filter(GRM_dat$grm, id1 != id2)
var(off_diag$grm)
diag <- filter(GRM_dat$grm, id1 == id2) %>% 
	mutate(cidB3030 = GRM_dat$id$V1) %>%
	mutate(IID = GRM_dat$id$V2)
var(diag$grm)

pdf(paste0("kinship/descriptives/meth_hist1.pdf"))
hist(diag$grm, xlab=paste0("IBS relatedness"), breaks = 30)
dev.off()

pdf(paste0("kinship/descriptives/meth_hist2.pdf"))
hist(off_diag$grm, xlab=paste0("IBS relatedness"), breaks = 30)
dev.off()

# -------------------------------------------------------
# Examine outliers in the data
# -------------------------------------------------------

Tukey <- function(df, column) {
	three_IQR <- 3 * IQR(df[, column], na.rm = T)
	quartiles <- summary(df[, column])[c(2, 5)]
	lo_lim <- as.numeric(quartiles[1] - three_IQR)
	up_lim <- as.numeric(quartiles[2] + three_IQR)

	large_vals <- df[df[, column] > lo_lim, ]
	vals <- large_vals[large_vals[, column] < up_lim, ]

	df <- df[rownames(df) %in% rownames(vals), ]
	return(df)
}

non_outliers <- Tukey(diag, "grm")
nrow(non_outliers)
outliers <- filter(diag, !(id1 %in% non_outliers$id1))
nrow(outliers) # 31 outliers

head(outliers)
write.table(outliers$IID, file = "outliers_IID.txt", quote = F, row.names = F, sep="")


# -------------------------------------------------------
# Genetic kinship matrix only: examine relatedness
# -------------------------------------------------------

#examine if any individuals >0.025 relatedness
max(off_diag$grm)  #0.253
length(off_diag$grm[off_diag$grm >0.025]) #30

#find off_diag values >0.025
rel <- which(off_diag$grm > 0.025)

#find which ids match to values in GRM_dat$id
rel_id <- off_diag[rel,]

#extract id of 1 individual in each related pair
rel_id1 <- rel_id$id1

#create list to exclude from analysis
ids <- GRM_dat$id[rel_id1,]
write.table(ids$V1, file = "cidlist_gen025.txt", quote = F, row.names = F, sep="")
