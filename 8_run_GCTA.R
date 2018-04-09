#load gcta software
#module add apps/gcta_1.91.0beta

library(tidyverse)

setwd("/panfs/bt14")

### Create text files of 2 kinship matrices to include
kinship_files <- c("kinship/methyl_FOM", "data/gen_FOM")
write.table(kinship_files, file ="kinship/methgen_FOM.txt", quote = F, row.names = F, col.names = F, sep = "\n")

# -------------------------------------------------------
# Set file input/output parameters
# -------------------------------------------------------

kinships <- "kinship/methgen_FOM.txt"
phens <- "BMI"
parameters <- expand.grid(phen = phens, kinship = kinships)

parameters$mgrm <- 2
parameters$phenfile <- "phen/pheno_FOM_BMI.txt"
parameters$kinshipfile <- "kinship/methgen_FOM.txt"
parameters$outfile <- "estimates/comb_FOM_BMI"
parameters$qcovar <- "phen/FOM.qcovar"


runGcta <- function(parameters, row, mgrm=2)
{
	p <- parameters[row,]
	if(file.exists(paste(p$outfile, ".hsq", sep=""))) return(NULL)
	if (mgrm > 1) p$kinshipfile <- p$kinshipfile
	grm <- ifelse(mgrm>1, "mgrm", "grm")
	cmd <- paste("gcta --", grm, " ", p$kinshipfile, " --reml --reml-no-constrain --reml-no-lrt --pheno ", p$phenfile, " --qcovar ", p$qcovar, " --out ", p$outfile, " --thread-num 10 --reml-maxit 500", sep="")
	system(cmd)
}

# -------------------------------------------------------
# run GCTA
# -------------------------------------------------------

runGcta(parameters, row = 1, mgrm = 2)
