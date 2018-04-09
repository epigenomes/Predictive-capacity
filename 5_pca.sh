# -------------------------------------------------------
# Generating principle components in ARIES
# -------------------------------------------------------

cd /panfs/bt14/data

# -------------------------------------------------------
# Create ID file
# -------------------------------------------------------

# Extract all the ID's from the "FOM" timepoint
grep $FOM ARIES_samplesheet.txt | awk ' { print $1$11" "$1$11 }' > FOM.txt

# -------------------------------------------------------
# Run PCA
# -------------------------------------------------------

data="/panfs/bt14/data/gen_FOM"


# Get snp list with no long range LD regions
awk -f ld.awk ${data}.bim > no_ld.txt

# Get independent SNPs excluding any long range LD regions
plink --bfile $gen_FOM --exclude no_ld.txt --indep 100 5 1.01 --out indep

# Calculate PCs
plink --bfile $data --keep $FOM.txt --extract indep.prune.in --pca 20 --out $FOM_pcs

#keep only individuals with DNA methylation data available