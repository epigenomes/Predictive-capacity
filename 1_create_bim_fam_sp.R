# -------------------------------------------------------
# Extract data
# -------------------------------------------------------

# Files needed
# ARIES betas
# ARIES samplesheet
# fdata

rm(list = ls())

library(tidyverse)
library(lme4)

setwd("/panfs/bt14")

# ARIES betas
load("mum_betas5.rdata")
# ARIES samplesheet
samplesheet <- read.delim('epi_data/ARIES_samplesheet.txt', stringsAsFactors = FALSE)
# fdata
load("epi_data/fdata_new.RData")


# Need to remove individuals who have multiple read outs at the same time point

timepoints <- c("FOM")

dat <- filter(samplesheet, time_point %in% timepoints)

fam <- dat %>%
	mutate(PID = 0) %>%
	mutate(MID = 0) %>%
	mutate(PHENO = 0) %>%
	dplyr::select(cidB3030, Sample_Name, PID, MID, Sex, PHENO)

remove <- names(which(table(fam$cidB3030) != length(timepoints)))

print(paste0("Number of people removed = ", length(remove)))

fam <- filter(fam, !(cidB3030 %in% remove))

# Sanity check 
stopifnot(as.numeric(names(which(table(fam$cidB3030) != length(timepoints)))) == 0)

unique(fam$Sex)
fam[["Sex"]] <- ifelse(fam$Sex == "F", 2, 1)
head(fam)

bim <- fdata.new %>%
	mutate(gd = 0) %>%
	mutate(a1 = "A") %>%
	mutate(a2 = "T") %>%
	select(CHR, TargetID, gd, COORDINATE_37, a1, a2) %>%
	filter(TargetID %in% rownames(beta5))

unique(bim[["CHR"]])
bim[["CHR"]] <- as.character(bim[["CHR"]])

bim[bim[["CHR"]] == "X", "CHR"] <- "23"
bim[bim[["CHR"]] == "Y", "CHR"] <- "24"

bim[["CHR"]] <- as.numeric(bim[["CHR"]])
bim[["COORDINATE_37"]] <- as.numeric(bim[["COORDINATE_37"]])


createSp <- function(bim, fam, dat)
{
	stopifnot(nrow(bim) == nrow(dat))
	stopifnot(nrow(fam) == ncol(dat))
	bim2 <- bim[order(bim$CHR, bim$COORDINATE_37), ]
	index <- match(bim2$TargetID, bim$TargetID)
	dat2 <- dat[index, ]
	return(list(sp=dat2, bim=bim2, fam=fam))
}

dat <- as.data.frame(beta5) %>%
	dplyr::select(one_of(fam[["Sample_Name"]]))


# -------------------------------------------------------
# Regress batch out
# -------------------------------------------------------
# Batch vars = Slide, BCD_plate, BCD_id, MSA4Plate_id
# From ARIES_samplesheet.txt

filt_sample <- filter(samplesheet, Sample_Name %in% colnames(dat))

write.table(filt_sample, file="epi_data/filt_samplesheet.txt")

index <- match(colnames(dat), filt_sample$Sample_Name)
filt_sample <- filt_sample[index, ]

# Sanity check
stopifnot(all(filt_sample$cidB3030 == colnames(dat)))

dat_resid <- apply(dat, 1, function(x) {resid(lmer(as.numeric(x) ~ (1 | filt_sample$BCD_plate)))})
rownames(dat_resid) <- filt_sample$cidB3030

stopifnot(colnames(dat) == fam[["cidB3030"]])


# -------------------------------------------------------
# Creating BIM, FAM and SP files
# -------------------------------------------------------

sp <- createSp(bim, fam, dat)

#change names
write.table(sp$sp, file = paste0("methyl_FOM.sp"), quote = F, col.names = F, row.names = F)
write.table(sp$bim, file = paste0("methyl_FOM.bim"), quote = F, col.names = F, row.names = F)
write.table(sp$fam, file = paste0("methyl_FOM.fam"), quote = F, col.names = F, row.names = F)


# ------------------------------------------------------------
# List of ids with methylation data (to keep for genetic GCTA)
# ------------------------------------------------------------
fam <- read.table(paste0("kinship/methyl_FOM.fam"), header = F)
keep <- data.frame(col1=fam$V1, col2=fam$V1)
write.table(keep, file = "data/keep.txt", quote = F, row.names = F, sep=" ")
