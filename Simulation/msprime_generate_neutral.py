#! /bin/bash/python3
# msprime coalescent simulation; generate neutral genetic variation for SLiM forward sim

import sys
import tskit
import msprime
import pyslim
import numpy as np
import pandas as pd

# command line parameters
trees = sys.argv[1]
seed = sys.argv[2]

### coalescent simulation in msprime - ancestry + neutral mutations, output tree sequence for SLiM
# simple demographic model for ancestral teosinte, scaled down 10x
demog_model = msprime.Demography()
demog_model.add_population(initial_size=12278)

# reading in chr. 4 recombination map
col_names = ['start', 'rate']
map_chr4 = pd.read_csv("Chr4_slim_recombmap_msprime_sort.bed", sep='\t', names = col_names)
positions = map_chr4['start']
rates = map_chr4['rate']
rate_map = msprime.RateMap(position=positions, rate=rates[1:])

# simulate ancestry for 656 samples (bottleneck population size/starting size for SLiM)
ots = msprime.sim_ancestry(
        samples=656,
        demography=demog_model,
        recombination_rate=rate_map,
        sequence_length=246510001,
        random_seed=seed)

# Annotate individuals, nodes, etc.
ats = pyslim.annotate(ots, model_type="WF", tick=1, stage="late")

# Add SLiM mutations
mut_model = msprime.SLiMMutationModel(type = 1, slim_generation = 1) # neutral mutations are type 1
ts = msprime.sim_mutations(
							ats,
							rate=4.93e-10, # mutation rates are same as recombination rates, so p/u = 1
							model=mut_model,
							keep=True,
							random_seed=seed
)


print(f"The tree sequence now has {ts.num_mutations} mutations, at "
      f"{ts.num_sites} distinct sites.")
tables = ts.tables
ts_metadata = tables.metadata
ts_metadata["SLiM"]["model_type"] = "WF"
tables.metadata = ts_metadata
ts = tables.tree_sequence()

# Change time units to generation
t = ts.dump_tables()
t.time_units = 'generations'
ts2 = t.tree_sequence()

# Outputs the trees
ts2.dump(trees)
