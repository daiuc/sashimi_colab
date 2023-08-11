#!/bin/bash

prep_region="chr10:99989193-100989193"
top_snp="chr10:99514301:T:C"

plot_region="chr10:100300000-100400000"


# data prep
./sashimi_ingredients.py \
    -V ZOD14598_AD_GRM_WGS_2021-04-29_chr10.recalibrated_variants.leftnorm.filtered.AF.chr10_99489193_101489193.test_samples.rename.new.recode.vcf \
    -B <(realpath *.bw) \
    -C counts_clu_109523.txt \
    -P individuals.txt \
    -T template.ini \
    -S $top_snp \
    -R $prep_region \
    -E 1000 \
    -O ./plot \
    && echo "data prep done" \
    && echo "plot sashimi: $plot_region" \
    && /scratch/midway3/chaodai/miniconda3/envs/pygenometracks/bin/pyGenomeTracks \
            --tracks plot.ini \
            --region  $plot_region \
            --title "QTL:${prep_region}|TOP:${top_snp}"\
            --width 30 \
            --trackLabelFraction .1 --fontSize 12 \
            -out plot_sashimi.pdf
