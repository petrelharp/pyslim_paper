import tskit
import msprime
import pyslim

ts = tskit.load("/Users/shyamag/Desktop/recap_example.trees")
demography = msprime.Demography.from_tree_sequence(ts, initial_size=100)
demography.add_migration_rate_change(
	time=ts.metadata['SLiM']['tick'] + 100,
	rate=0.1, source="p2", dest="p1",
)
rts = pyslim.recapitate(initial_state=ts, demography=demography,recombination_rate=1e-8)