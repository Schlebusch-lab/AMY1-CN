#!/usr/bin/env bash

#### MODULES

module load PLINK


#### DIRECTORIES

VCFs=".../project_amylase/simulations/01_simulated_VCFs_and_CN/VCFs"
MDIST=".../project_amylase/simulations/02_simulated_MDIST_PCA/EVEN_SAMPLING/MDIST"
PCA=".../project_amylase/simulations/02_simulated_MDIST_PCA/EVEN_SAMPLING/PCA"


#### MAKE PCA AND MDIST

for FILE in "$VCFs"/*.vcf; do

    [ -e "$FILE" ] || continue

    NAME=$(basename "$FILE")
    NAME="${NAME%.vcf}"

    echo "******* Processing $NAME..."

    plink \
        --vcf "$FILE" \
        --maf 0.05 \
        --pca \
        --out "$PCA/${NAME}_maf0.05"

    plink \
        --vcf "$FILE" \
        --maf 0.05 \
        --distance square 1-ibs \
        --out "$MDIST/${NAME}_maf0.05"

done

rm -f "$MDIST"/*.log "$MDIST"/*.nosex \
      "$PCA"/*.log "$PCA"/*.nosex


