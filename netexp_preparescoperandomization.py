# netexp_prepare_scope_randomization

import json
import glob
import os
import random
import collections
import copy
import math
import pandas as pd

def load_json(fname):
    with open(fname) as f:
        return json.load(f)

def read_formatted_jsons(INDIR):
    generation_dfs = []
    stats_dicts = []
    
    for domain in os.listdir(INDIR):
        for fname in glob.glob(INDIR+domain+"/*.json"):

            org_id = os.path.basename(fname).strip(".json")

            with open(fname) as f:
                datajson = json.load(f)

            datajson["stats"]["org_id"] = org_id
            datajson["stats"]["domain"] = domain
            datajson["stats"]["path"] = fname

            stats_dicts.append(datajson["stats"])
            generation_dfs.append(pd.DataFrame(datajson["generations"]))
            
    return generation_dfs, stats_dicts

def read_molecular_masses(cpd_dir):
    
    kegg_masses = dict()
    for fname in os.listdir(cpd_dir):
        cjson = load_json(cpd_dir+fname)
        kegg_masses[cjson[0]['entry_id']] = float(cjson[0]['exact_mass'])

    return kegg_masses ## Tuple of (cpd_id, mass)

def initialize_randomization_for_organism_scope(org_max_scope,kegg_masses):
        
    org_cpd_masses = list()
    for cpd in org_max_scope:
        org_cpd_masses.append((cpd,kegg_masses[cpd]))

    sorted_org_cpd_masses = sorted([i for i in org_cpd_masses if i[1]!=0], key=lambda x: x[1], reverse=True)
    org_cpd_zeros = [i for i in org_cpd_masses if i[1]==0]

    new_sorted_org_cpd_masses = copy.deepcopy(sorted_org_cpd_masses)
    for i,tup in enumerate(org_cpd_zeros):
        new_sorted_org_cpd_masses.insert(random.randrange(len(sorted_org_cpd_masses)+i), tup)
    
    return new_sorted_org_cpd_masses

def mix_it_up(tup_list,beta,n):
    for _ in range(n):
        i1, i2 = random.sample(range(len(tup_list)),2)  ## 2 random indices
        if swap_random(tup_list[i1],tup_list[i2],beta) == True:
            tup_list[i1], tup_list[i2] = tup_list[i2], tup_list[i1] 
        else: 
            pass
    return tup_list

def swap_random(tup1,tup2,beta):

    diff = tup1[1]-tup2[1]

    ## Determine probability of flipping
    if tup1[1]==0 or tup2[1]==0:
        p = .5 # flip with 50% chance
    
    elif diff > 0:
        p = math.e**(-diff/float(beta))

    elif diff <= 0:
        p = 1.0

    ## Flip and return
    if random.random()<p:
        return True
    else:
        return False


n_swaps = 10000
n_randomizations = 1000
beta = 20

cpd_dir = "kegg/2018-09-25/compound/"
edge_json_dir = "results/formatted/ph_edge_jsons/"

kegg_masses = read_molecular_masses(cpd_dir)
generation_dfs, stats_dicts = read_formatted_jsons(edge_json_dir)
ph_archaea_stats_dicts = [i for i in stats_dicts if i["domain"]=="archaea"]
# ph_bacteria_stats_dicts = [i for i in stats_dicts if i["domain"]=="bacteria"]
# archaea_max_scopes = [i["scope_compounds"] for i in ph_archaea_stats_dicts]
# bacteria_max_scopes = [i["scope_compounds"] for i in ph_bacteria_stats_dicts]

for org in ph_archaea_stats_dicts:

    randomized_cpd_lists = list()
    for _ in range(n_randomizations):
        new_sorted_org_cpd_masses = initialize_randomization_for_organism_scope(org["scope_compounds"],kegg_masses)
        new_sorted_org_cpd_masses = mix_it_up(new_sorted_org_cpd_masses,beta,n_swaps)
        randomized_cpd_list = [i[0] for i in new_sorted_org_cpd_masses]
        
        if randomized_cpd_list not in randomized_cpd_lists:
            randomized_cpd_lists.append(randomized_cpd_list)

    ph_archaea_stats_dicts["randomized_cpd_lists"] = randomized_cpd_lists


