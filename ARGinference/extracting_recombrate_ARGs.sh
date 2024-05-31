#!/bin/bash
echo "$@"
echo "$#"
x=$1
name=$2
directory=$3

echo $x
echo $name
echo $directory

#####normalizing recombination rate by branch length
q=1000000

export PATH=/programs/argweaver/bin:$PATH
module load python/2.7.15
export PYTHONPATH=/programs/argweaver:$PYTHONPATH

#converting to bed file
for i in $(seq 1 $q $x); do smc2bed-all -s 600 ${name}_$i_$(($i+$q))_ARGweaver_output; done

#extracting branch lengths per time & MCMC iteration
for i in $(seq 1 $q $x); do arg-summarize -a ${name}_$i_$(($i+$q))_ARGweaver_output.bed.gz --branchlen-per-time --log-file ${name}_$i_$(($i+$q))_ARGweaver_output.log > ${name}_$i_$(($i+$q))_branchlength_per_time.bed; done

#extracting ARGs
for i in $(seq 1 $q $x); do for p in {600..1000..10}; do echo "smc2arg ${name}_$i_$(($i+$q))_ARGweaver_output.$p.smc.gz ${name}_$i_$(($i+$q))_ARGweaver_output.$p.arg" ; done; done > extract_recombs.txt
parallel --jobs 30 < extract_recombs.txt

#pulling recombination breakpoints from ARG
for i in $(seq 1 $q $x); do python3 /home2/rke27/msprime_sims_argweaver/recomb_date_extract.py ${name}_$i_$(($i+$q))_ARGweaver_output.*.arg ; done

#sorting recombination events from a certain time point
for i in $(seq 1 $q $x); do sort -k1,1n -k2,2n ${name}_$i_$(($i+$q))_ARGweaver_output_122_recombs_mcmc_it.bed > ${name}_$i_$(($i+$q))_ARGweaver_output_122_recombs_mcmc_it_sorted.bed; done

#matching recombs to branch length
for i in $(seq 1 $q $x); do awk '{print $4"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24}' ${name}_$i_$(($i+$q))_branchlength_per_time.bed | tail -n +4 - | sort-bed - | bedtools intersect -a - -b ${name}_$i_$(($i+$q))_ARGweaver_output_122_recombs_mcmc_it.bed -c > ${name}_$i_$(($i+$q))_ARGweaver_output_122_recombs_mcmc_it_sorted_countreps.bed; done

##calculating recombination rate per MCMC iter and then summing over all
for i in $(seq 1 $q $x); do python3 ${directory}/mean_branch_len.py ${name}_$i_$(($i+$q))_ARGweaver_output_122_recombs_mcmc_it_sorted_countreps.bed; done

cat ${name}_*_recombrate_perMCMC_recent.bed > ${name}_combined_recombrate_perMCMC_recent.bed

echo "all done!"
