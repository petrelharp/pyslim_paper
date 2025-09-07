import tskit, pyslim

ts = tskit.load("all.trees")

assert ts.num_populations == 6

sample_sets = [ts.samples(population=k) for k in range(6)]

x4 = ts.diversity(sample_sets[4], mode='node')
d4 = [sum(x4[s]) for s in sample_sets]
for k, d in enumerate(d4):
    print(f"Total in subpopulation {k} on paths between two samples from 4: {d}")
# check that 04 inherits from both 01 sometimes and 02 sometimes
# (so that reticulation worked)
assert d4[1] > 0
assert d4[2] > 0
# and that 04 does not inherit from 03 or 05
assert d4[3] == 0
assert d4[5] == 0
