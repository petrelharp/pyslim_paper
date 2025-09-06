# Overview

1. Run SLiM script to product tree sequence.
2. If necessary, use `merge_ts.py` to merge tree sequences for a multiply-infected individual.
2. For each new host, use `reset.py` to get the tree sequence ready to start a new run.


## The SLiM script

In `virus_in_host.slim`:

- `FOUNDERS` and `FOUNDING_TIME` are dictionaries whose keys are host IDs.
- `OUTPUT_HOST_IDS` and `TRANSMISSION_DAYS` are lists of the same length that say to whom and when
    to transmit individuals.

0. When initializing a seed simulation, saves an empty list to `FOUNDERS`
    and the current tick to `FOUNDING_TIME`.

1. When loading from a file, 

    * saves the `FOUNDERS` and `FOUNDING_TIME` to metadata so it will passed on;

2. When transmitting,

    * samples founding individuals for each output host
    * saves the pedigree IDs of these founders to `FOUNDERS` under the output host ID key,
        and also the current tick to `FOUNDING_TIME`, also under the output host ID key.
    * remembers these founders
    * creates a new population with id `HOST_ID`
    * transfers those individuals listed under `HOST_ID` in `FOUNDERS` to the new population

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

1. Finds the founding time and list of founders under `host_id` in metadata.
2. Computes the (tskit) time ago from the founding time and the current tick recorded in metadata.
3. Uses `pyslim.set_slim_state` to reset to that time and with those founders alive.
4. Writes out the tree sequence.

## The merge script

Given two tree sequences, `merge_ts.py` does:

1. Shfits the node and mutation times in each input so that time matches,
    by adding to the "time ago" of each the difference between their SLiM tick and the largest SLiM tick.

2. Merges the population tables so that any populations described in the second
    are also described in the first.

2. Constructs a node mapping from the second to the first:

    * Finds the shared founder nodes between the two,
        by taking the intersection of their `FOUNDERS` dictionaries,
        and then finding all nodes belonging to the individuals with those SLiM IDs.
    * Finds in each the ancestors of those shared founder nodes,
        as all those reachable by edges.
    * These are identical, so should be in the same order; verifying this,
        map from one to the other.

3. Unions the two with this mapping.

4. Copies missing keys from `FOUNDERS` and `FOUNDING_TIME` in the second's metadata to the union.

5. Outputs the result.
