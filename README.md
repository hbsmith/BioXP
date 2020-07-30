# BioXP
Network expansion on biological seed compounds

## Getting started
The easiest way to use this package is to check out the example notebooks. Here you can see the functionality that motivated the consturction of this package.

## Example notebooks 
Can be found in the directory `examples`. The data used to successfully run these example notebooks is not included in the BioXP respository, but these should give you an idea how to structure your code. 
### `expansion-encel_seeds.ipynb`
Example of how to run basic network expansions over many different organisms in a nested file structure. 

### `make_random_seeds_for_minimal_seed_expansions.ipynb`
Example of how to generate ordered seed sets for many different organisms in a nested file structure. 

### `expansion-find_minimal_seeds.ipynb`
Example of how to run minimal seed set network expansions over many different organisms in a nested file structure (remember that you'll need to have generated your ordered seed sets in order to do the minimal seed set expansions).

### `expansion-find_minimal_seeds-fixed_set.ipynb`
Just like the example above, except this time there's a function in the notebook that allows us to fix certain seeds to the front of the ordered seed set lists, so that they are not allowed to be removed. This allows you to figure out what the least number of compounds you need to add to an existing seed set is. 

## Organization of src
### BioXP.jl
This is the core of the package, containing the code for all functions which "do" the network expansion, and find minimal seed sets. It also contains helper functions for this, including a simple writing function (which is terribly bloated in terms of how much space it takes up since it writes jsons). All functions/types can be access through this file because all the functions/types which are meant to be called externally are imported and exported here. 

#### Important functions defined here:
- `expand` -- Does a network expansion
- `find_minimal_seed_set` -- Finds minimal seed sets

### dg.jl
Functions that restrict which reactions are accessible based on user-inputted free energy thresholds and other options.

#### Important functions defined here:
- `filter_reactions_by_dg` -- Returns `Vector{Bool}` of length of the number of reactions you have, indicating if the reaction is allowed or not--based on your specified dg constraints. Returns forward and backward vectors.

### format.jl
Functions that format output from the expansions/minimal seed sets.

#### Important functions defined here:
- `formatbioxpoutput` -- Takes the simple outputs from the expansions/minimal seed sets and adds a lot of information you may or may not care about, to avoid some "preprocessing" when doing data analysis.

### parser.jl
Functions that coerce data from the `ecg` repo output into appropriate inputs for the `BioXP` repo. Also serves as a check on types.

#### Important functions defined here:
- `readmaster` -- coerces the "master" file from ecg
- `readcompounds` -- coerces the compound dir from ecg
- `readids` -- coerces the seeds/targets
- `readkeyedids` -- coerces the seeds/targets (if they're keyed)

### randomize.jl
Functions that randomize orderings of compounds, which is useful for doing the minimal seed set expansions. 

#### Important functions defined here:
- `randomizecompounds` -- Return a vector of run results. Each run result is a randomized list of biosystem_compounds.

### structs.jl
Definitions of structs and types used in BioXP (these are basically the equivilent of python classes).

#### Important structs defined here:
- `Compound`, `Compounds`, `Reaction`, `Reactions`, `IDs`

## Documentation
Unfortunately there is none right now (I know, lame). But I hope the Julia code is relatively easy to read and the example files help show how the package is meant to be used.

## To do...
- Add example files
- It would be nice to move all the `ecg` generated data to an AWS bucket, and to store and access it using [quilt](https://docs.quiltdata.com/). Then I could easily pull the quilt data for the example files. 


