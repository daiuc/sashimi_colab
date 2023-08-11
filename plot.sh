#!/bin/bash

data_pre_region="chr10:99989193-100989193"
top_snp="chr10:99514301:T:C"

plot_region="chr10:100500000-100540000"
counts=counts_clu_113382.txt

vcf=ZOD14598_AD_GRM_WGS_2021-04-29_chr10.recalibrated_variants.leftnorm.filtered.AF.chr10_99489193_101489193.test_samples.rename.new.recode.vcf
sample_id_list=individuals.txt

clu=$(echo $counts | sed -n 's/[^_]*_\(.*\)\.txt/\1/p')
out_pdf=sashimi_${clu}.pdf

# data prep
./sashimi_ingredients.py \
    -V $vcf \
    -B <(realpath *.bw) \
    -C $counts \
    -P $sample_id_list \
    -T template.ini \
    -S $top_snp \
    -R $data_pre_region \
    -E 1000 \
    -O ./plot \
    && echo "data prep done" \
    && echo "plot sashimi: $plot_region" \
    && /scratch/midway3/chaodai/miniconda3/envs/pygenometracks/bin/pyGenomeTracks \
            --tracks plot.ini \
            --region  $plot_region \
            --title "QTL:${plot_region}|TOP:${top_snp}|$clu"\
            --width 30 \
            --trackLabelFraction .1 --fontSize 10 \
            -out $out_pdf \
    && echo "plot done"
