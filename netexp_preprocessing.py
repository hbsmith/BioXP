import json
import glob
import os

def load_json(fname):
    """
    Wrapper to load json in single line

    :param fname: the filepath to the json file
    """
    with open(fname) as f:
        return json.load(f)

def write_json(data,fname):
    """
    Wrapper to write json in single line

    :param data: the data to write out
    :param fname: the filepath to write out
    """
    with open(fname, 'w') as outfile:
        json.dump(data, outfile)

def get_org_enzlist_from_org_json(organism_json):

    organism_json = load_json(organism_json)
    
    ec_list = list()
    for d in organism_json['records']:
        ec_list.append(d['EnzymeID'].split("EC:")[1])
    
    return ec_list

def get_org_rxnlist_from_ec_list(ec_list,ec_rxn_link_json):
    
    ec_rxn_link_json =load_json(ec_rxn_link_json)

    rxn_list = list()
    for ec in ec_list:
        ec_key = 'ec:'+ec
        if ec_key in ec_rxn_link_json:
            rxn_list+=ec_rxn_link_json[ec_key]
        
    return rxn_list

def get_org_reactant_products_from_rxnlist(rxn_list,reaction_edges_json):

    rxn_edges_json = load_json(reaction_edges_json)
    
    org_rxn_edges_dict = {'products':{},'substrates':{}}
        
    for rxn in rxn_list:
        rxn_key = rxn.strip(':rn')
        if (rxn_key in rxn_edges_json['products']) and (rxn_key in rxn_edges_json['substrates']):
            org_rxn_edges_dict['products'][rxn_key] = rxn_edges_json['products'][rxn_key]
            org_rxn_edges_dict['substrates'][rxn_key] = rxn_edges_json['substrates'][rxn_key]
            
    return org_rxn_edges_dict

def convert_single_jgi_json_to_edge_json(organism_json,ec_rxn_link_json,reaction_edges_json):

    ec_list = get_org_enzlist_from_org_json(organism_json)
    rxn_list = get_org_rxnlist_from_ec_list(ec_list,ec_rxn_link_json)
    org_rxn_edges = get_org_reactant_products_from_rxnlist(rxn_list,reaction_edges_json)

    return org_rxn_edges

def main():

    ec_rxn_link_json = 'kegg/2018-09-25/links/enzyme_reaction.json'
    reaction_edges_json = 'kegg/2018-09-25/reaction_edges.json'
    outdir = 'jgi/2018-09-29/ph_edge_jsons/'

    dirs = ['bacteria']
    for d in dirs:
       
        fnames = glob.glob('jgi/2018-09-29/ph_jsons/%s/*.json'%d)
        print("Writing %s edge_jsons..."%d)
        for fname in fnames:

            org_rxn_edges = convert_single_jgi_json_to_edge_json(fname,ec_rxn_link_json,reaction_edges_json)

            outpath = outdir+d+'/'
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            outpath = outpath+os.path.basename(fname)
            write_json(org_rxn_edges,outpath)


if __name__ == '__main__':
    main()

# organism_json = '../jgi/Bacillus_pumilus_SAFR-032.json'
# ec_rxn_link_json = '../kegg/2018-09-25/links/enzyme_reaction.json'
# reaction_edges_json = '../kegg/2018-09-25/reaction_edges.json'

# ec_list = get_org_enzlist_from_org_json(organism_json)
# rxn_list = get_org_rxnlist_from_ec_list(ec_list,ec_rxn_link_json)
# org_rxn_edges = get_org_reactant_products_from_rxnlist(rxn_list,reaction_edges_json)

# write_dir = '../jgi/edges_test/'+os.path.basename(organism_json)
# write_json(org_rxn_edges,write_dir)

    
    