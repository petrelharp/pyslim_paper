# Bridging forward-in-time and coalescent simulations using pyslim

by Shyamalika Gopalan, Murillo Rodrigues, Ben Haller, and Peter Ralph

Code accompanying the sections:

1. Recapitation: [notebook](code/recapitation/recapitation.ipynb), [SLiM script](code/recapitation/recap_example.slim)
2. Generating initial diversity: [notebook](code/generating_diversity/generating.ipynb), [SLiM script](code/generating_diversity/reload_annotated.slim)
3. Generating genetic data: [notebook](code/generating_genetic_data/generating.ipynb)
4. Parallelizing multiple species: [notebook](code/parallelizing_multiple_species/parallelizing.ipynb), [SLiM script](code/parallelizing_multiple_species/simulate_branch.slim), [Makefile](code/parallelizing_multiple_species/parallel_sims.make)



5. Simulating pathogen evolution:
NOTE - we occasionally run into the tskit.union() error "Shared portions of the tree sequences are not equal" when merging the final tree sequences. This happens because there are individuals present in one of the final trees, but are missing in another. This is why we think this could happen (although it is very rare): 1) an individual is marked as a founder in one lineage, 2) that individual's parent happens to be a founder in a parallel lineage, but 3) the parent's tree sequence is not directly ancestral to the child's tree sequence in the phylogeny. The error can be resolved by setting check_shared_equality=False, but we have not wanted to make this the default behaviour of the script. -SSG
