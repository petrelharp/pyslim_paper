# Bridging forward-in-time and coalescent simulations using pyslim

by Shyamalika Gopalan, Murillo Rodrigues, Ben Haller, and Peter Ralph

Code accompanying the sections:

- Section 2: **Recapitation**

    * [SLiM script](code/recapitation/recap_example.slim)
    * [notebook](code/recapitation/recapitation.ipynb) with python code

- Section 3: **Generating initial diversity**

    * [notebook](code/generating_diversity/generating.ipynb) with python code
    * [SLiM script](code/generating_diversity/reload_annotated.slim)

- Section 4: **Generating genetic data**
    
    * [notebook](code/generating_genetic_data/generating.ipynb)

- Section 5: **Parallelizing forard-in-time simulations of multiple species**

    * [SLiM script](code/parallelizing_multiple_species/simulate_branch.slim)
    * [Makefile](code/parallelizing_multiple_species/parallel_sims.make)
    * [notebook](code/parallelizing_multiple_species/parallelizing.ipynb)

- Section 5: **Metapopulation dynamics with simulation networks**

    * [README](code/metapopulation_dynamics/README.md)
    * [SLiM script](code/metapopulation_dynamics/virus_in_host.slim)
    * [Makefile](code/metapopulation_dynamics/Makefile)
    * python script to [reset the state](code/metapopulation_dynamics/reset.py)
        of a tree sequence for SLiM
    * python script to [merge different tree sequences](code/metapopulation_dynamics/merge_ts.py)
        which uses [these functions](code/metapopulation_dynamics/merge_utils.py)
