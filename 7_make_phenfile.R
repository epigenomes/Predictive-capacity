library(tidyverse)

setwd("/panfs/bt14")

# -------------------------------------------------------
# Set input/output parameters
# -------------------------------------------------------

phen_dat <- read.delim("/phen/pheno_data.txt")
#set negative (error/missing) values to NA to exclude from analysis
phen_dat$fm1ms111[phen_dat$fm1ms111 < 0] <- NA

phens <- "BMI"
parameters <- expand.grid(phen = phens)
parameters$phenfile <- "phen/pheno_FOM_BMI.txt"

# -------------------------------------------------------
# Select IDs for analysis
# -------------------------------------------------------

# Read in the fam file 
fam <-read.table("data/gen_FOM.fam", header = TRUE)

#keep only individuals with genetic and methylation data present
ids <-read.table("data/gen_FOM.grm.id")
fam <- subset(fam, (IID %in% ids$V1))

#METHYLATION ANALYSIS ONLY: remove outlier IDs
outliers <- read.delim('outliers_IID.txt', stringsAsFactors = FALSE)
fam <- subset(fam, !(IID %in% outliers$x))

#GENETIC ANALYSIS ONLY: remove related IDs
related <- read.delim('cidlist_gen025.txt', stringsAsFactors = FALSE)
fam <- subset(fam, !(IID %in% related$x))
# -------------------------------------------------------

dat <- subset(phen_dat, (phen_dat$cidB3030 %in% fam$FID))
#include only individuals with phenotype data present
dat <- dat[complete.cases(dat$fm1ms111), ]
dim(dat)

# -------------------------------------------------------
# Create phenotype file
# -------------------------------------------------------

#create file in which trait is normally distributed

row = 1
pheno_col = 9
dat = dat
dat$fm1ms111 <- as.numeric(dat$fm1ms111)

createPhen <- function(parameters, row, dat, fam, pheno_col)
{
	p <- parameters[row, ]
	if(file.exists(p$phenfile)) return(NULL)
	require(GenABEL)
	d <- dat
	value <- names(d[pheno_col])
	f <- fam[, 1:2]
	d <- merge(d, f, by.x="cidB3030", by.y="FID")
	d <- subset(d, select=c(cidB3030, IID, get(value)))
	d[[value]] <- rntransform(d[[value]])
	d <- subset(d, !duplicated(cidB3030))
	write.table(d, file=p$phenfile, row=F, col=F, qu=F)
	return(NULL)
}

createPhen(parameters, row = 1, dat = dat, fam = fam, pheno_col = 9)
