# Change into working directory
cd $PBS_O_WORKDIR
qsub -I -l nodes=1:ppn=16,walltime=08:00:00

# Execute code

# -------------------------------------------------------
# Create DNA methylation profiles matrix
# -------------------------------------------------------

./ldak5.linux --sp methyl_FOM \
 --calc-kins-direct methyl_FOM \
 --SNP-data NO \
 --ignore-weights YES \
 --power -0.25 


# -------------------------------------------------------
# Create genetic profiles matrix
# -------------------------------------------------------

./ldak5.linux --bfile gen_FOM \
 --calc-kins-direct gen_FOM \
 --keep FOM.txt \
 --SNP-data YES \
 --ignore-weights YES \
 --power -0.25 

#keep only individuals with DNA methylation data available