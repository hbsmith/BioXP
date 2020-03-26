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


- 
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

### 10-4-2018 Thursday
- Been reading large scale reconstruction (Borenstein 2008), trying to figure out if this is a possible way to determine seed sets for the enceladus networks.
- Sara mentioned today that my paper on network projections was really important. that it means something important that in some coarse grainings that biochemistry is scale free, but not in others. i argued that it was just a consequence of the fact that you're representing different information in differnt projections, but she thinks that's important.
- Need to finish understanding this papers methods and the other papers i printed

### 10-7-2018 Sunday
- Going to run the seed set randomizations 1000 variants of 
- 10,000 exchange operations to generate one randomized list. For each network, 1000 such lists (or should I do only 100 to start? Can see how long 1000 takes for the 28 archaea)

Bacteria energy info:
```
collections.Counter([i['EnergySource'] for i in ph_bacteria_records])
Counter({u'Chemoautotroph': 1,
         u'Chemoheterotroph': 20,
         u'Chemoheterotroph, Photoautotroph': 1,
         u'Chemolithoautotroph': 1,
         u'Chemolithotroph': 1,
         u'Chemoorganoheterotroph': 2,
         u'Chemoorganotroph': 27,
         u'Chemoorganotroph, Heterotroph': 2,
         u'Heterotroph': 31,
         u'Lithoheterotroph': 1,
         u'Methylotroph': 1,
         u'Phototroph': 1,
         u'zzz': 177})
```
Bacteria habitat info (101 of 266 bacteria contain the word "Aquatic")
```
collections.Counter([i['Habitat'] for i in ph_bacteria_records])
```

Archaea energy info:
```
collections.Counter([i['EnergySource'] for i in ph_archaea_records])
Counter({u'Chemoorganotroph': 10,
         u'Heterotroph': 1,
         u'Organotroph': 2,
         u'zzz': 15})
```
Archaea EcosystemCategory info:
```
collections.Counter([i['EcosystemCategory'] for i in ph_archaea_records])
Counter({u'Aquatic': 24, u'Terrestrial': 1, u'Wastewater': 1, u'zzz': 2})
```

### 2018-10-9/10 Tues/Wed

From selfish-metabolism paper:

> The ecosystem,
> like the reductive chemo-autotrophic
> cell, takes in CO2, H2 (or other reductants),
> NH3, H2S, and H3PO4, to create
> this core which consists of twenty
> amino acids, four ribonucleotides, four
> deoxyribonucleotides, and a few sugars,
> polar lipids, and cofactors. The
> rest of biomolecular complexity, and
> most of the distinction among individuals
> and species, is defined by combinatorial
> assembly of building blocks
> from this core.

pH actually should be 11-12, not 9-11. But principles still hold

- Created new directory for organized versions of randomization jsons (`min_seeds_final`)
  - Added in all archaea; all outdirectory bacteria, and b2s, b3s

- `compounds_new[0]` in step 1 is the same as `stats['scope_seeds']`



```python
def read_formatted_jsons_streamlined(INDIR,encel):
    generation_dfs = []
    stats_dicts = []
    
    for domain in os.listdir(INDIR):
        for org in os.listdir(os.path.join(INDIR,domain)):
            for fname in glob.glob(os.path.join(INDIR,domain,org,"*.json")):
        
                d = dict()
        
                with open(fname) as f:
                    datajson = json.load(f)   
                
                d["seed"] = os.path.basename(fname).strip(".json")
                d["org_id"] = org
                d["domain"] = domain
                d["path"] = fname
                
                d["network_compounds"] = datajson["stats"]["scope_compounds"]
                d["network_reactions"] = datajson["stats"]["scope_reactions"]
                d["seed_compounds"] = datajson["stats"]["scope_seeds"]
                d["seed_compounds_on_enceladus"] = [list(set(clist) & set(encel)) for clist in d["seed_compounds"]]
                d["target_compounds"] = datajson["stats"]["scope_targets"]
                d["target_compounds_in_seeds"] = [list(set(clist) & set(d["target_compounds"])) for clist in d["seed_compounds"]]
                d["n_generations"] = max(datajson["generations"].keys())
                d["scope_compounds"] = datajson["generations"][d["n_generations"]]["compounds_cumulative"]
                d["scope_reactions"] = datajson["generations"][d["n_generations"]]["reactions_cumulative"]

                stats_dicts.append(datajson["stats"])
                generation_dfs.append(pd.DataFrame(datajson["generations"]))
            
    return generation_dfs, stats_dicts
```

## 2019-08-20 Tuesday

- Found the code used to make figure 1 in the submitted paper:

  `BioXP/jupyter/netexp-plotting-PART1.ipynb`

  I made a copied version to make edits to:

  `BioXP/jupyter/netexp-plotting-PART1-revisions.ipynb`

  so now I can mess around with showing the data in diff ways

### 2019-09-9 

- Need to add `element_conservation` field to reactions **done**
- Need to remove those reactions from the `reaction_edges.json` KEGG file (or rather duplicate that version of KEGG and remove these from the duplicated version) **done. i have explicitly renamed the old files, appending `no_element_conservation` to the name. the new conserved files have the original name without the addendum.**
- Then I can rerun the `netexp_preprocessing.py` script with this new `reaction_edges_element_conserved.json` file as an input, and I shouldn't have to change any other inputs--just need to specify a new output directory. **done**
- Rerun expansions with new `ph_edge_jsons`.
  - New results in `results/simple_2019-09-09/`
    - Starting with kegg (`kegg_edge_json_P/`)
    - Also did no P version
- Now I need to figure out why I call a function which doesn't exist (prepare_seeds and get that to work, or change how i import seeds.)

### 2019-09-10

- reran network expansions for all `ph_edge_json`s and `ph_edge_json_P`s
- reran postprocessing for all of these, plus the `kegg_edge_json` and `kegg_edge_json_P` (this needs to be done before randomization can be prepared, since randomization is based on the formatted files)
- reran `netexp_preparescoperandomization.py` on the archea/bacteria files
- I manually split all the same archaea and bacteria from before into 5 archaea files and 6 bacteria files within the `2019-09-09/ph_edge_jsons/` dir 
- now i need to find the minimal seeds for these
  - started 11 concurrent runs at 6:50pm on Tues sept 10
  - last batch of 6 finished running at ~1pm Wed sept 11

### 2019-09-12

- for some reason, water is NOT being included as a scope target even for seeds which start with water. 
  - it's because for some reason water is missing from the `targets/Freilich09.json` file.
  - I think I should just assume that everything makes water.... **update, it doesn't. check the 2019-09-12 jupyter notebook in the encxp dir**
  - 

## 2020-03-24

Trying to understand why I can't modify the temperature involved in the eQuilibrator calcuations--the official explanation is:

- > The group contribution method enables us to approximate ΔfG of compounds at a particular temperature (the temperature at which they were measured) [[JM08\]](http://equilibrator.weizmann.ac.il/static/classic_rxns/faq.html#jm08). As the change in free energy is defined as ΔG = ΔH - TΔS and we don’t know the value of ΔS in most cases, we cannot predict how changes in temperature will affect ΔfG.

- Need to understand this better by reading there paper--as there should be estimates for how to improve this.

Kristin shared some thermodynamic resources with me:

- HKF method (https://www.hindawi.com/journals/geofluids/2019/5750390/)
- orchyd asu database (no longer maintained and only one person can access at a time http://orchyd.asu.edu/)
- https://chnosz.net/vignettes/obigt.html

And some other ones I found:

- https://gitlab.com/equilibrator/component-contribution (general API)
- https://gitlab.com/equilibrator/equilibrator-api (eQuilibrator)

Today I also finished producing a new seeds file which reflects the most recent estimates of enceladus's composition. Those seeds are now under `data/seeds/encel_papers_2019.json` . The script used to create it was `jupyter_new/update_seeds.ipynb`. 

## 2020-03-25

Trying to install the local equilibrator API and it recommends installing in a virutalenv. This is allowing me to import the package into a python3 instance within the directory, but for some reason when I try to import equilibrator while in a jupyter lab notebook it's not finding the package. **ok resolved this by not using virtualenv**

### Equilibrator and changing temperature

Not possible. Because normally $\Delta _rG’^o$ (=$-RT ln(K’)$) is calculated by measuring the apparent equilibrium constant K' (that is the concentrations of species at equilibrium) at a particular temperature. For more explanation, see the `Noor et al. 2013 SI Section 1 Training data`:

> Nearly all Gibbs energy measurements, for compounds and reactions in aqueous solutions at near-room temperature, are derived from the equilibrium constants of enzyme-catalyzed reactions. Typically, an enzyme that speciﬁcally catalyzes a certain reaction is puriﬁed and added to a medium that contains the reaction substrates. After the reaction reaches equilibrium, all reactant concentrations are measured. The equilibrium constant K ′ is deﬁned as the ratio between the product of all product concentrations and the product of all substrate concentrations. Since there is no easy way to distinguish between pseudoisomers of the same compound, the concentration of every reactant is actually the sum of all its protonation 2
>
> states. Therefore, K ′ is the apparent equilibrium constant, which is related to the standard transformed Gibbs energy of reaction (∆ r G ′◦ [1]). The problem with using this data as-is lies in the fact that K ′ and ∆ r G ′◦ depend on the aqueous environment (e.g. pH, ionic strength, and pMg). The measurements listed in TECRDB span a wide range of pH and ionic strength values, and many of the reactions have only been measured in non-standard conditions.

Also, from `Jankowski et al 2008:`

>  The experimentally measured  values reported in these references were captured under a variety of temperature and pH conditions. Only data captured within one pH unit and 15 K of the chosen reference state of pH 7 and 298 K was utilized

### Ionic strength

Another question is: What ionic strength to use? Default on equilibrator is 0.1M but that's not a good reason to choose it.

Maybe .2-.6 based on the results of the paper by `Hsu et al 2015`:

> The existence of silica nanoparticles also provides strict constraints on the salinity of Enceladus’ subsurface waters because silica colloids aggregate and precipitate quickly at high ionic strength[12](https://www.nature.com/articles/nature14262#ref-CR12),[13](https://www.nature.com/articles/nature14262#ref-CR13). The critical coagulation concentration of NaCl at pH 9 is 2% or ∼0.3 M (1.5% or ∼0.2 M at pH 10, 4% or ∼0.6 M at pH 8)[13](https://www.nature.com/articles/nature14262#ref-CR13). 

`Zolotov 2007` shows 0.4-0.1 Molal maybe:

> The solution pH ranges from 8 to 11 ([Figures 1](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007GL031234#grl23670-fig-0001) and [2b](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007GL031234#grl23670-fig-0002)). The salinity (2–20 g/kg H2O) and ionic strength (0.04–0.1 molal) are higher at high‐*T* ([Figure 2c](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007GL031234#grl23670-fig-0002)) and low‐*W*/*R* conditions, but less than in Earth's seawater (∼35 g kg−1 and 0.7 molal).

### Temperature of enceladus

Fig 3 from `Hsu et al 2015` shows that temp is unknown and is probably somewhere between 90 and 0 C:

![Figure 3](41586_2015_Article_BFnature14262_Fig3_HTML.jpg)

### More compounds from enceladus...

I found another good paper on compounds on enceladus: `Khawaja et al 2019` which has actual tables (!) of compounds or estimated compounds based on how ionization (or fragmentation?) occured in laboratory experiements. See tables 1 and 2. See also the SM.

### Todo

- map equilibrator delta Gs to KEGG reactions for various combinations of pH and ionic strength 

- check to see how network expansion in my code works if i want to limit reaction directionality based on deltaG and pH

- figure out how to download some of the organisms that Shawn recommended to me

  - also redownload pH organisms, and check if i can sample phyla like i helped yoko do?

- think about  how to analyze results from my paper based on whether or not organisms are:

  - aerobic/anaerobic
  - autotroph/heterotroph
  - photo/chemotroph?
  - different metabolisms (eg. sulfur reducing or methanogens)

  