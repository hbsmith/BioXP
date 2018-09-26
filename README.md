# BioXP
Network expansion on biological seed compounds

## kegg
The kegg directory `2018-09-25` contains directories storing the KEGG entry files for pathways, enzymes, reactions, and compounds. Despite the folder name, all directories except `pathway` were retrieved from KEGG on `2017-12-01`. The `.json` files containing the ids and names of all entries across all directories were retrieved on `2018-09-25`. 

The only file meant to be called at the moment is `retrieve_kegg.py`. The functions (in order) in this file will download kegg data, get detailed entries, map them to each other, and organize/clean them.

The `kegg/date/links` directory includes mapping between all the various kegg ids. Be careful, because the `pathway` mappings contains pathways in the form of `path:map` and `path:rn` and `path:ec` depending on the mapping type.
