#  /bin/bash/python3
import sys
import tskit
import msprime
import pyslim
import pandas as pd
import numpy as np

trees = sys.argv[1]
seed = sys.argv[2]
vcf = sys.argv[3]

#loading in trees from SLiM
orig_ts1 = tskit.load(trees)

#adding neutral mutations to SLiM output
# Fix old SLiM format
orig_ts = pyslim.update(orig_ts1)

# Verify everything coalesces - should be 1 root
orig_max_roots = max(t.num_roots for t in orig_ts.trees())
print(f"Maximum number of roots: {orig_max_roots}\n")

# reading in chr. 4 recombination map
col_names = ['start', 'rate']
map_chr4 = pd.read_csv("Chr4_slim_recombmap_msprime_sort.bed", sep='\t', names = col_names)
positions = map_chr4['start']
rates = map_chr4['rate']
rate_map = msprime.RateMap(position=positions, rate=rates[1:])

# Get individuals that are alive at the end of the SLiM sim
alive_inds = pyslim.individuals_alive_at(orig_ts, 0)

rng = np.random.default_rng(seed = int(seed))

# Keep just 60 individuals that were alive at the end of the simulation
keep_indivs = rng.choice(alive_inds, 60, replace=False)

keep_nodes = []
for i in keep_indivs:
  keep_nodes.extend(orig_ts.individual(i).nodes)

# Simpify the ts to just have those nodes for the individuals
sts = orig_ts.simplify(keep_nodes, keep_input_roots=True)

# Add the neutral mutations on
mut_model = msprime.SLiMMutationModel(type = 1, slim_generation = 1552) # neutral mutations are type 1, slim_gen to ensure domestication ages positive
ts = msprime.sim_mutations(
           sts,
           rate=4.93e-10,
           model=mut_model,
           keep=True,
           random_seed=seed
)
# Check that it worked (there should be more mutations)
print(f"The tree sequence now has {ts.num_mutations} mutations,\n"
      f"and mean pairwise nucleotide diversity is {ts.diversity():0.3e}.")


### Check metadata and write tree sequence out
ts.dump(trees + "_neutralmutations_added_final.trees")

indivlist = []
for i in pyslim.individuals_alive_at(ts, 0):
    ind = ts.individual(i)
    if ts.node(ind.nodes[0]).is_sample():
       indivlist.append(i)
       # if one node is a sample, the other should be also:
       assert ts.node(ind.nodes[1]).is_sample()

# Write VCF
with open(vcf, "w") as vcffile:
    pyslim.convert_alleles(pyslim.generate_nucleotides(ts)).write_vcf(vcffile, individuals=indivlist)
