import networkx as nx
import json
import os
import glob
import cPickle as pickle

def combine_all_compound_jsons_into_single_dict(compound_directory):

    compound_dict = dict()

    for path in glob.glob(compound_directory+"*.json"):
        # print path
        # outpath = outdirname+os.path.basename(path)

        with open(path) as data_file:    
            data = json.load(data_file)[0] #dict of single compound

            data["entry_type"] = "compound"

            compound_dict[data["entry_id"]] = data

    return compound_dict

def combine_all_reaction_jsons_into_single_dict(reaction_directory):

    reaction_dict = dict()

    for path in glob.glob(reaction_directory+"*.json"):
        # print path
        # outpath = outdirname+os.path.basename(path)

        with open(path) as data_file:    
            data = json.load(data_file)[0] #dict of single reaction

            data["entry_type"] = "reaction"

            reaction_dict[data["entry_id"]] = data

    return reaction_dict

def combine_compound_and_reaction_dicts(compound_dict,reaction_dict):
    return dict(compound_dict.items() + reaction_dict.items())

def create_networkx_graphs_from_edgelist_json(edge_filepath,attribute_dict):

    with open(edge_filepath) as data_file:
        data = json.load(data_file)

    G = nx.Graph()

    G.add_nodes_from(data["substrates"])
    G.add_nodes_from(data["products"])

    ## Add edges
    for k, v in data["substrates"].items():
        G.add_edges_from(([(k, t) for t in v]))
    for k, v in data["products"].items():
        G.add_edges_from(([(k, t) for t in v]))

    nx.set_node_attributes(G, attribute_dict)

    return G

def main():
    compound_directory = "newdata/20171201/compounds/"
    reaction_directory = "newdata/20171201/reactions_detailed/"
    edge_filepath = "newdata/20171201/reaction_edges.json"
    outpath = "newdata/20171201/kegg_reaction_compound_graph_with_attributes.pkl"

    compound_dict = combine_all_compound_jsons_into_single_dict(compound_directory)
    reaction_dict = combine_all_reaction_jsons_into_single_dict(reaction_directory)
    attribute_dict = combine_compound_and_reaction_dicts(compound_dict,reaction_dict)
    G = create_networkx_graphs_from_edgelist_json(edge_filepath,attribute_dict)

    print "Pickling networkx graph..."

    nx.write_gpickle(G, outpath)



    # create_networkx_graphs_from_edgelist_json(edge_filepath)
    # add_substrates_products_stoichiometry_to_reaction_jsons(dirname,outdirname)
    # test_parsing(outdirname)


if __name__ == '__main__':
    main()