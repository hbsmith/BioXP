import Bio
import Bio.TogoWS as TogoWS
from Bio.KEGG import REST, Enzyme, Compound, Map
import json


def write_kegg_id_name_jsons(outdir):
    """
    :param outdir: directory to write to 
    """
    ## Gather map, ec, rn, cpd detailed files
    kegg_types = ["pathway","enzyme","reaction","compound"]
    for kegg_type in kegg_types:
        
        id_name_dict = dict()
        raw_list = REST.kegg_list(kegg_type)
        id_name_list = [s.split('\t') for s in raw_list.read().splitlines()]
        for i in id_name_list:
            id_name_dict[i[0]] = i[1]

        outpath = outdir+kegg_type+'.json'
        with open(outpath, 'w') as outfile:   
            json.dump(id_name_dict, outfile, indent=2)


def main():
    outdir = "kegg/2018-09-25"
    # write_kegg_id_name_jsons(outdir)

if __name__ == '__main__':
    main()


## Get list of each of map, ec, rn, cpd ids and names

## Loop through lists to retrieve detailed files

## Add detailed information for reactions

## Retrieve links between all combinations of map, ec, rn, cpd

