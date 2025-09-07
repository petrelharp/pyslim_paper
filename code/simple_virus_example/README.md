# Overview

1. Run SLiM script to product tree sequence.
2. If necessary, use `merge_ts.py` to merge tree sequences for a multiply-infected individual.
3. For each new host, use `reset.py` to get the tree sequence ready to start a new run.

Transmission diagram:
```
1               00
                /\
10            01  \
              /\   \
15          03  \   \
                |    \
20              |     02
                |    /
                |   /
                |  /
35              04
                |
60              05
```

## The SLiM script

In `virus_in_host.slim`:

- `FOUNDING_TIME` is a dictionary whose keys are host IDs.
- `OUTPUT_HOST_IDS` and `TRANSMISSION_DAYS` are lists of the same length that say to whom and when
    to transmit individuals.

0. When initializing a seed simulation, saves the current tick to `FOUNDING_TIME`.

1. When loading from a file, 

    * saves the `FOUNDING_TIME` to metadata so it will passed on;

2. When transmitting,

    * samples founding individuals for each output host
    * saves the current tick to `FOUNDING_TIME`, also under the output host ID key.
    * remembers these founders
    * creates a new subpopulation with id `HOST_ID`
    * transfers the founders to the new subpopulation

3. On the sampling day,

    * chooses some individuals to remember

4. At the end,

    * kills everyone
    * outputs the tree sequence.

## The reset script

Running
```
python reset.py infile.trees host_id outfile.trees
```
does:

1. Finds the founding time under `host_id` in metadata.
2. Finds the founding individuals as the only individuals in the population whose id is `host_id`.
3. Computes the (tskit) time ago from the founding time and the current tick recorded in metadata.
4. Removes individuals with duplicate pedigree IDs from the list of founders.
5. Uses `pyslim.set_slim_state` to reset to that time and with those founders alive.
6. Writes out the tree sequence.

## The merge script

Given two tree sequences, `merge_ts.py` does:

1. Shfits the node and mutation times in each input so that time matches,
    by adding to the "time ago" of each the difference between their SLiM tick and the largest SLiM tick.

2. Merges the population tables so that any populations described in the second
    are also described in the first.

3. Updates individual metadata and flags for any individuals contained in both,
    identified by unique combinations of birth population, pedigree_id,
    so that these agree between the two tree sequences for shared individuals.

4. Constructs a node mapping from the second to the first, by identifying those
    nodes sharing unique (birth population, slim_id) pairs.

5. Unions the two with this mapping.

6. Copies missing keys from `FOUNDING_TIME` in the second's metadata to the union.

7. Outputs the result.
