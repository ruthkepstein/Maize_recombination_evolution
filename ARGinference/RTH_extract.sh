#!/bin/bash

#first arg is chr #
#second arg is end of chr
#third arg is name of corresponding files
export PATH=/programs/argweaver/bin:$PATH
module load python/2.7.15
export PYTHONPATH=/programs/argweaver:$PYTHONPATH

echo "$@"
echo "$#"
chr=$1
x=$2
name=$3

echo $chr
echo $x
echo $name
q=1000000

awk -v chr="$chr" '$1 == chr' /home2/rke27/meiotic_genes_v4_coordinates_sorted.bed > meiotic_genes_v4_3cols_chr${chr}.bed
awk -v chr="$chr" '$1 == chr' /home2/rke27/upstream_meiosis_genes_v4_50kb_3cols_sorted.bed > upstream_meiosis_genes_v4_50kb_3cols_sorted_chr${chr}.bed

#extracting --rth across all non-recombined regions
for i in $(seq 1 $q $x); do arg-summarize -a ${name}_$i_$(($i+$q))_ARGweaver_output.bed.gz -r ${chr}:$i-$(($i+$q)) --rth --mean --log-file ${name}_$i_$(($i+$q))_ARGweaver_output.log > ${name}_$i_$(($i+$q))_rth.txt ; done
cat ${name}_*_rth.txt > ${name}_combined_rth.txt
sort -k1,1 -k2,2n ${name}_combined_rth.txt | sed '/^#/d' - > ${name}_combined_rth_sorted.txt

#intersecting recombination genes + null gene set to see if rth measures are significantly different
bedtools intersect -a meiotic_genes_v4_3cols_chr${chr}.bed -b ${name}_combined_rth_sorted.txt -wb > meiosis_genes_chr${chr}_rth.bed
bedtools intersect -a upstream_meiosis_genes_v4_50kb_3cols_sorted_chr${chr}.bed -b ${name}_combined_rth_sorted.txt -wb > null_genes_chr${chr}_rth.bed
