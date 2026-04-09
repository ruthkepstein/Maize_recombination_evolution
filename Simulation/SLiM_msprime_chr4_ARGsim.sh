#!/bin/bash -l

# Set up random seed
bash -c RANDOM=42;
SEED=`echo $RANDOM`

# Script paths
#put comment elsewhere
MSPRIME1=/home2/rke27/msprime_sims_variation/msprime_generate_neutral # coalescent sims
SLIM=/home2/rke27/msprime_sims_variation/slim_maize_chr4 # SLiM script
MSPRIME2=/home2/rke27/msprime_sims_variation/msprime_addmutations # simplify, add muts, output tree sequence

# First trees file name
TREES_INITIAL="$(pwd)/${SEED}_initial.trees"
TREES_SLIM="$(pwd)/${SEED}_slim.trees"

## Run msprime simulation
#module load conda/msprime
python3 $MSPRIME1.py "$TREES_INITIAL" "$SEED"

## Run SLiM
echo $TREES_INITIAL
/programs/SLiM-4.3/bin/slim -d "trees='$TREES_INITIAL'" -seed $SEED "$SLIM.slim"

## Run msprime/pyslim simplify, add mutation -
python3 $MSPRIME2.py "$TREES_SLIM" "$SEED" "${SEED}.vcf"

gzip "${SEED}.vcf" > "${SEED}.vcf.gz"
