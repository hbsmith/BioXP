### Gathering Data

- Use KEGG_REST_API to download list of all KEGG:

  - enzymes http://rest.kegg.jp/list/enzyme
  - compounds http://rest.kegg.jp/list/compound
  - reactions http://rest.kegg.jp/list/reaction
  - (and current KEGG stats/info) http://rest.kegg.jp/info/kegg
  - This can also be done for the other KEGG database entries, listed here: http://www.kegg.jp/kegg/docs/dbentry.html

- Use the TogoWS database + Biopython's TogoWS module, and loop over each database entry (from the database.txt files download above) while extracting from the TogoWS database (which is kept in sync with the KEGG database. e.g. It is updated at least every day from spot checking it).

  - http://togows.dbcls.jp/
  - See module documentation here. There are no examples, but data can be pulled simply by calling the `.entry` attribute. This returns a `handle` . Calling `handle.read()` will return the data, which can be written out to a file in the standard way python files are written out.
    - http://biopython.org/DIST/docs/api/Bio.TogoWS-module.html
  - this is done using the `generate_entry_files.py` to generate json files for `reactions`,`compounds`, and `enzymes`

- Parse the `reactions` `.json` files in order to extract and store (in a better way than how they are stored from pulling) the compound IDs associated with a reaction's substrates and products. These are stored in the `reactions_detailed` directory. These jsons are updated with 5 additonal fields if they do not contain glycans as part of the reaction. The fields are: `substrates`,`products`,`substrate_stoichiometries`,`product_stoichiometries`,`glycans`. All fields except `glycans` contains lists of strings. `glycans` is a boolean value. The lits of `substrates` and `substrate_stoichiometries` (for example) contain lists whose elements match up the compound ID to the stoichiometries. Some stoichiometries are letters (`n` or `m`). This will form the basis of the information from which to construct an updated KEGG network, and the data which can be used expanding networks from seed compounds.

  - This is done using a lot of regex in the file `add_substrates_products_stoichiometry_to_reaction_jsons.py`

  - A `.json` was created using `generate_reaction_substrate_product_edges_json.py` that contains a structure like:

  - ``````json
    {
      "products": {
        "R07881": [
          "C14101", 
          "C00010"
          ],...
      }, 
      "substrates": {
        "R07881": [
          "C14120", 
          "C00001"
        ],...
      }
    }
    ``````


- ​
- Figure out the best way to store the reaction/compound data for easily converting to a networkx graph object.
  - I'm thinking the easiest way would be to add fields directly to the existing reaction jsons:
    - substrate field
    - product field

### Creating reaction scopes

- Use objects to construct labeled graphs
- Create lists of which compounds can be reached from each other compound in x number of steps

### Notes

There appear to be duplicate reactions amongst all the reaction ids. For example, `R04193` is the same as `R06183`. I believe this only happens when a reaction involves glycans. Glycans appear to be degenerately stored as both compounds and glycans, in both the compound and glycan database. For example, `C03323 == G10497`, and this compound is known as all of the following: `"1,2-beta-D-fructan:1,2-beta-D-fructan 1F-beta-D-fructosyltransferase"`.

### Log

- Dec. 13, 2017: I've identified 3 compound IDs which are present in some reactions' substrate and product lists, but which do not have official entries in KEGG COMPOUND. e.g. they are in:

- ```python
  substrates = nx.get_node_attributes(G,"substrates")

  products = nx.get_node_attributes(G,"products")

  ```

- but not in:

  `reactions = nx.get_node_attributes(G,"reactions")`

  ​

- These are:

  - C009991 (found in R11044), labeled as Ferrocytochrome b51, no hyperlink in REACTION db. Not listed in the ENZYME db (1.14.19.6). I suspect this is a typo based on this, and this enzymes other reaction (R11043). It should likely be C009999 (Ferrocytochrome b5).
  - C03061 (found in R08652), labeled as 2-Methylhomogentisate. Broken hyperlink in REACTION db. Not listed in the ENZYME db ([1.14.14.40](http://www.kegg.jp/dbget-bin/www_bget?ec:1.14.14.40)). I do not think this is a typo.
  - C21593 (found in R11803), labeled as D-Erythronate. Broken hyperlink in REACTION db. No hyperlink in the ENZYME db. I do not think this is a typo.

- Added:

- ```python
  ## remove reactions which contain compounds that don't exist
  entriesToRemove = ['R11803', 'R11044', 'R08652']
  for k in entriesToRemove:
      reactions.pop(k, None)
  ```

  to the function `advance_expansion_1_generation`

- reformatted the removal of reactions to functions (1/18/18)

- Ran through multiple seed sets

### 9-20-18

Tasks:

- See how different network expansion is without O2 reactions
- Visualize networks in expansion (see Fig. 1 in Goldford et al.)
- Figure out AA enrichment/depletion relative to composition of enzymes used by life in an Enceladus like environment. If I assume I need all amino acids in order to be able to catalyze everything in KEGG, just need to see if I can make all amino acids with the enceladus network expansion
  - to do this for individual organisms, it requires identifying the genes which code the ECs which catalyze the reactions. And then looking at the presence/absence of AA in the enzymes coded by those genes.
- Look into annotated genomic data for organisms which can hitch rides on clean spacecraft--do network expansion on their reaction sets
  - Looks like there have been some of these completely sequenced and annotated! 
  - Check how reliant these organisms are on other organisms to survive. e.g. do people know that these are autotrophs or chemoautotrophs? how do people know that? How can i interpret the quantitative results of network expansion to figure out if this organism is viable or not? e.g. how to quantify if organism has access to the right molecules needed for growth and reproduction?
- Look into most "minimal" , "resiliant", or "autotrophic" organisms and do network expansion on their reaction sets
- Look into research on whether the molecules we've observed from enceldaus's plume, e-rings, etc are accurate predictions of the available substances in the ocean.
- Look into the metal availability on spacecraft--could these fufull the requirements necessary for metal cofactors?

Notes:

- Only L amino acids used in proteins and by cells. What about beta amino acids? These seem superfluous (via John), don't need to condider anything but L for purposes of life.
- 

Sara notes:
- Try different genomes/metagenomes on enceladus and goldford. Can you get organisms that create amino acids, or "core metabolism", or anything?
- Might try to pick the largest organisms, or most complete genomes
- Combine all clean room or rad-hardened organisms in one network expansion since individual ones don't seem to be working.

### 9-24-2018 Monday

- Use Biopython and KEGG REST API to get mappings between:
  - paths (maps) and reactions (and vice versa)
  - paths (maps) and compounds (and vice versa)
  - reactions and compounds (and vice versa)
  - And download description data and lists of each of these things.
- Goldford used KEGG map IDs to test for pathways enrichment, as well as KEGG MO module IDS to test for module enrichement, but only pathway enrichment was plotted (this is a higher level)
- e.g. `http://rest.kegg.jp/link/reaction/map00471`
- rewrite network expansion algorithms using matrices and make sure they yeild same results. Can run them with Pypy to go faster this way.
- Figure out how the F Goldford did http://www.biostathandbook.com/fishers.html fisher's exact test on the pathway enrichment. Need to talk to a bio person... Then it is straightforward to do Benjamini-Hochberg multiple comparison test to make sure False Discovery Rate (FDR) is small http://www.biostathandbook.com/multiplecomparisons.html 
- Check handwritten notes for incorporating phosphorous compounds in Enceladus network expansion (looks like best bet is planting P2O5 or H(PO3)- and assuming it came from schriberite or apatite, and will exist even if it is at low concentrations in enceladus's ocean. Steve Desch thinks weathering will be small for an ocean at pH of 9-11 but Shock thinks it will still exist, just below detection limits for Cassini.)
- Run network expansion for different seed sets to see just how much the compounds affect the scopes here (eg it seems like Goldford and Enceladus have same scope when using all of KEGG).
- Filter JGI archaea, bacteria, eukarya, or metagenomes to only look at complete, large, genomes that can exist at high pH in order to pair down omic data for test on enceladus's seed chemistry. How do scopes differ with some plausible Earth seed sets? Or differ when using all KEGG reactions? (this could be a proxy for Earth since you could assume other organisms catalyze the other reactions)
- Find more seqeunced organisms which could survive sterile clean room environments and do network expansion on these organisms.
- Where are the files that `kegg_mapping.py` generates being used? Need to recreate lookup dicts using Biopython KEGG REST methods.
- What are life's acidity limits? Looks high enough (pH 10-11)
- What about making my last chapter more of a proposed method paper with examples? Can talk about challenges but also benefits. Section on case study on Enceladus: sections on organisms most likely to hitch a ride, most likely to survive, and most likely to produce earth's catalytic repoitoire. 

### 9-25-2018 Tuesday

- Start new repository for the network expansion and biochemical data acquisition scripts. How to best organize? What are my desired inputs/outputs?
  - Inputs: 
    - mappings among compounds, reactions, enzymes, maps (pathways). 
    - Detailed descriptions of each entry of these datasets.
    - list of JGI organisms to pull data from (or better yet their EC lists) OR list of metadata properties to use to filter out organisms
- Organism properties:
  - Ecosystem: Acquatic
  - Energy source: Chemoorganotroph, lithotroph, autotroph, Chemolithoautotroph
  - Temperature range: Psychrophile
  - Salinity: Halophile
  - pH: Alkaliphile (these can grow in pH ~10)
- Enceladus proprties:
  - pH of 9-11 (Seewald 2017 Science)

### 9-26-2018 Wednesday

- cleaned up code in the bioxp repository. Everything should be ready to write the actual network expansion algorithm at this point. 
- No files are saved in a "network" or "graph" format, but I have all the relationship data in the various json files. 

### 9-27-2018 Thursday

- For individual organisms, I need to write a function which:
  -  parses the individual organisms enzymes
  -  maps them to reactions
  -  calls the main `reaction_edges.json` to find substrates/products of these reactions to include in scope and build inputs to the network expansion algorithm (e.g. the Reactant and Product matrices). Other inputs include seed compounds, master target compounds, and intersection of current scope with master target compounds to find organismal specific target compounds.
  -  When I loop through the algorithm, I will store the single step compounds, reactions, maps added, and the running total of them.
-  Need to update my pages outline with info from the two most recent papers I read, and the papers I skimmed today about algorithms for finding all minimal seed sets for specific target compounds for organisms.