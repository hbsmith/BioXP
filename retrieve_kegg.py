import Bio
import Bio.TogoWS as TogoWS
from Bio.KEGG import REST, Enzyme, Compound, Map
import json
import os

def write_kegg_id_name_jsons(keggdir):
    """
    Get jsons of each of map (pathway), ec, rn, cpd ids and names

    :param keggdir: directory to write to 
    """
    
    if not os.path.exists(keggdir):
        os.makedirs(keggdir)

    kegg_types = ["pathway","enzyme","reaction","compound"]
    for kegg_type in kegg_types:
        
        ## Retreive all entry ids and names
        id_name_dict = dict()
        raw_list = REST.kegg_list(kegg_type)
        id_name_list = [s.split('\t') for s in raw_list.read().splitlines()]
        for i in id_name_list:
            id_name_dict[i[0]] = i[1]

        ## Write json of all entry ids and names
        outpath = keggdir+kegg_type+'.json'
        with open(outpath, 'w') as outfile:   
            json.dump(id_name_dict, outfile, indent=2)

def write_kegg_entry_jsons(keggdir):
    """
    Get jsons of entry for each item in map (pathway), ec, rn, cpd ids and names

    :param keggdir: directory to read master lists from, and to write to
    """

    ## Get entries on all kegg types
    kegg_types = ["pathway","enzyme","reaction","compound"]
    for kegg_type in kegg_types:

        ## Read list of all kegg ids
        inpath = keggdir+kegg_type+'.json'
        with open(inpath) as data_file:    
            data = json.load(data_file)#[0]

        ## Create dir to store entries in
        outdir = keggdir+kegg_type
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        ## Grab each entry in list
        for i,entry in enumerate(data):
            
            entry_id = entry.split(":")[1]
            entry_fname = entry+".json"

            print "Saving (verifying) %s entry %s of %s (%s)..."%(kegg_type,i+1,len(data),entry_id)

            while entry_fname not in os.listdir(outdir):
                try:
                    handle = TogoWS.entry(kegg_type, entry_id, format="json")
                    with open(outdir+'/'+entry_fname, 'a') as f:
                        f.write(handle.read())
                except:
                    pass

            break



def main():
    keggdir = "kegg/2018-09-25/"
    # write_kegg_id_name_jsons(keggdir)
    write_kegg_entry_jsons(keggdir)

if __name__ == '__main__':
    main()


## Get list of each of map, ec, rn, cpd ids and names

## Loop through lists to retrieve detailed files

## Add detailed information for reactions

## Retrieve links between all combinations of map, ec, rn, cpd

