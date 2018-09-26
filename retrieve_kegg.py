import Bio
import Bio.TogoWS as TogoWS
from Bio.KEGG import REST, Enzyme, Compound, Map
import json
import os
import glob
import re

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

    raise Warning("UNRELIABLE. TOGOWS CALL SOMETIMES RETURNS NULL ENTRY VALUES. NEEDS DEBUGGING.")

    ## Get entries on all kegg types
    kegg_types = ["compound"] #["pathway","enzyme","compound", "reaction"
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

def write_kegg_reaction_detail(keggdir):
    """
    Get add reaction details in convenient fields for rn jsons.
    Output in new "reaction_detailed" dir.

    :param keggdir: directory to read master lists from, and to write to
    """

    reaction_path = keggdir+'reaction/'
    reaction_detail_path = keggdir+'reaction_detailed/'
    if not os.path.exists(reaction_detail_path):
        os.makedirs(reaction_detail_path)
    
    for path in glob.glob(reaction_path+"*.json"):
        
        outpath = reaction_detail_path+os.path.basename(path)

        with open(path) as data_file:    
            data = json.load(data_file)
            
            equation = data[0]["equation"]

            # print equation

            if re.search(r'(G\d+)',equation) == None: ## Only find entries without glycans

                for i, side in enumerate(equation.split(" <=> ")):
                    # print i, side
                    # if i==0:
                    #   reactants = []
                    #   reactants_stoichiometry = []
                    # elif i==1:
                    #   products = []
                    #   products_stoichiometry = []

                    compounds = []
                    stoichiometries = []

                    ## match (n+1) C00001, (m-1) C00001 or similar
                    matches = re.findall(r'(\(\S*\) C\d+)',side)
                    # print matches
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = re.search(r'(\(\S*\))',match).group(1)
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match 23n C00001, m C00001 or similar
                    matches = re.findall(r'(\d*[n,m] C\d+)',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = re.search(r'(\d*[n,m])',match).group(1)
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match C06215(m+n), C06215(23m) or similar
                    matches = re.findall(r'(C\d+\(\S*\))',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = re.search(r'(\(\S*\))',match).group(1)
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "3 C00002" or similar (but NOT C00002 without a number)
                    matches = re.findall(r'(\d+ C\d+)',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = match.split(' '+compound)[0]# re.search(r'(\(\S*\))',match).group(1)
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "C00001 "at the start of the line (no coefficients)
                    matches = re.findall(r'(^C\d+) ',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = '1'
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "+ C00001 " (no coefficients)
                    matches = re.findall(r'(\+ C\d+ )',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = "1"
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "+ C00001" at the end of the line (no coefficients)
                    matches = re.findall(r'(\+ C\d+$)',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = "1"
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    ## match "C00001" which is at the start and end of the line
                    matches = re.findall(r'(^C\d+$)',side)
                    if len(matches) != 0:
                        for match in matches:
                            compound = re.search(r'(C\d+)',match).group(1)
                            stoichiometry = "1"
                            
                            compounds.append(compound)
                            stoichiometries.append(stoichiometry)

                    if i==0:
                        # print "Substrates!"
                        data[0]["substrates"] = compounds
                        data[0]["substrate_stoichiometries"] = stoichiometries
                    elif i==1:
                        # print "Products!"
                        data[0]["products"] = compounds
                        data[0]["product_stoichiometries"] = stoichiometries

                    # print compounds
                    # print stoichiometries
                    assert len(compounds) == len(stoichiometries)
                    # print "="*70
                    data[0]["glycans"] = False




            else:

                data[0]["glycans"] = True


        with open(outpath, 'w') as outfile:
            
            json.dump(data, outfile, indent=2)


#### TESTS--NEED TO BE MOVED TO APPROPRIATE TESTING FILES
def test_parsing(keggdir):

    dirname = keggdir+'reaction_detail/'
    with open(dirname+"rn:R00001.json") as data_file:    
        data = json.load(data_file)[0]
        assert set(data["substrates"]) == set(["cpd:C00404","cpd:C00001"])
        assert set(data["products"]) == set(["cpd:C02174"])

        assert set(data["substrate_stoichiometries"]) == set(["1","n"])
        assert set(data["product_stoichiometries"]) == set(["(n+1)"])

    with open(dirname+"rn:R00006.json") as data_file:    
        data = json.load(data_file)[0]
        assert set(data["substrates"]) == set(["cpd:C00900","cpd:C00011"])
        assert set(data["products"]) == set(["cpd:C00022"])

        assert set(data["substrate_stoichiometries"]) == set(["1"])
        assert set(data["product_stoichiometries"]) == set(["2"])

    with open(dirname+"rn:R00008.json") as data_file:    
        data = json.load(data_file)[0]
        assert set(data["substrates"]) == set(["cpd:C06033"])
        assert set(data["products"]) == set(["cpd:C00022"])

        assert set(data["substrate_stoichiometries"]) == set(["1"])
        assert set(data["product_stoichiometries"]) == set(["2"])

    with open(dirname+"rn:R00011.json") as data_file:    
        data = json.load(data_file)[0]
        assert set(data["substrates"]) == set(["cpd:C19610","cpd:C00027","cpd:C00080"])
        assert set(data["products"]) == set(["cpd:C19611","cpd:C00001"])

        assert set(data["substrate_stoichiometries"]) == set(["1","2"])
        assert set(data["product_stoichiometries"]) == set(["2"])

    with open(dirname+"rn:R05624.json") as data_file:    
        data = json.load(data_file)[0]
        assert set(data["substrates"]) == set(["cpd:C00001","cpd:C06215"])
        assert set(data["products"]) == set(["cpd:C06215"])

        assert set(data["substrate_stoichiometries"]) == set(["1","(m+n)"])
        assert set(data["product_stoichiometries"]) == set(["(m)","(n)"])

    with open(dirname+"rn:R07640.json") as data_file:    
        data = json.load(data_file)[0]
        assert set(data["substrates"]) == set(["cpd:C00002","cpd:C00046","cpd:C00046"])
        assert set(data["products"]) == set(["cpd:C00020","cpd:C00013","cpd:C00046"])

        assert set(data["substrate_stoichiometries"]) == set(["1","(m)","(n)"])
        assert set(data["product_stoichiometries"]) == set(["1","(n+m)"])

def write_kegg_link_mappings_jsons(keggdir):
    """
    Get jsons of each of mapping between each of: map (pathway), ec, rn, cpd ids and names

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



            



def main():
    keggdir = "kegg/2018-09-25/"

    # write_kegg_id_name_jsons(keggdir)
    # write_kegg_entry_jsons(keggdir)
    write_kegg_reaction_detail(keggdir)
    # test_parsing(keggdir)

if __name__ == '__main__':
    main()


## Get list of each of map, ec, rn, cpd ids and names

## Loop through lists to retrieve detailed files

## Add detailed information for reactions

## Retrieve links between all combinations of map, ec, rn, cpd

