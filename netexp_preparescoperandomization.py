# netexp_prepare_scope_randomization

import json
import glob
import os
import random
import collections
import copy

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

org_cpd_masses = list()
for cpd in archaea_max_scopes[0]:
    org_cpd_masses.append((cpd,kegg_masses[cpd]))

sorted_org_cpd_masses = sorted([i for i in org_cpd_masses if i[1]!=0], key=lambda x: x[1], reverse=True)
org_cpd_zeros = [i for i in org_cpd_masses if i[1]==0]

new_sorted_org_cpd_masses = copy.deepcopy(sorted_org_cpd_masses)
for i,tup in enumerate(org_cpd_zeros):
    new_sorted_org_cpd_masses.insert(random.randrange(len(sorted_org_cpd_masses)+i), tup)


cpd_dir = "kegg/2018-09-25/compound/"
edge_json_dir = "results/formatted/ph_edge_jsons/"

generation_dfs, stats_dicts = read_formatted_jsons(edge_json_dir)
ph_archaea_stats_dicts = [i for i in stats_dicts if i["domain"]=="archaea"]
ph_bacteria_stats_dicts = [i for i in stats_dicts if i["domain"]=="bacteria"]

kegg_masses = read_molecular_masses(cpd_dir)