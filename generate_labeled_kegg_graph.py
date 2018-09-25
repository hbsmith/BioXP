import re
import cPickle as pickle

class Compound(object):
    def __init__(self, entry):
        self.entry = entry
        self.names = []
        self.formula = ''
        self.reactions = []
        self.pathways = []
        self.enzymes = []
        self.chirality = None


def populate_compounds(filename):
    """
    Read in KEGG compounds.txt file, organize into dictionary of compound objects with associated info

    Parameters: 
        filename [.txt file]
            KEGG's compounds.txt

    Returns:
        compounds [dict of objects]
            Dictionary of Compound objects populated with info from compounds.txt

    """

    compounds = dict()

    more_names = False
    more_reactions = False
    more_pathways = False
    more_enzymes = False

    with open(filename) as f:
        for line in f:

            if not line.startswith(' '):
                more_names = False
                more_reactions = False
                more_pathways = False
                more_enzymes = False

                if line.startswith("ENTRY"):
                    
                    entry = line.split()[1]
                    compounds[entry] = Compound(entry)
                
                elif line.startswith("NAME"):

                    first_name = line.split('NAME')[1]
                    first_name = first_name.strip().rstrip()

                    if first_name[-1] == ';':
                        compounds[entry].names.append(first_name[:-1])
                        more_names = True
                    else:
                        compounds[entry].names.append(first_name)
                        more_names = False
                
                elif line.startswith("FORMULA"):

                    formula = line.split()[1]
                    compounds[entry].formula = formula
                
                elif line.startswith("REACTION"):

                    for i in line.split()[1:]:
                        compounds[entry].reactions.append(i) 
                        more_reactions = True
                
                elif line.startswith("PATHWAY"):

                    first_pathway = line.split()[1]
                    compounds[entry].pathways.append(first_pathway)
                    more_pathways = True
                
                elif line.startswith("ENZYME"):

                    for i in line.split()[1:]:
                        compounds[entry].enzymes.append(i) 
                        more_enzymes = True


            else: # line.startswith(' '):
                """NAME"""
                if more_names == True:
                    name = line.split('  ')[-1]
                    name = name.strip().rstrip()

                    if name[-1] == ';':
                        compounds[entry].names.append(name[:-1])
                        more_names = True
                    else:
                        compounds[entry].names.append(name)

                """REACTION"""
                if more_reactions == True:
                    
                    for i in line.split():
                        compounds[entry].reactions.append(i) 


                """PATHWAY""" 
                if more_pathways == True:
                    compounds[entry].pathways.append(line.split()[1]) 

                """ENZYME""" 
                if more_enzymes == True:
                    
                    for i in line.split():
                        compounds[entry].enzymes.append(i) 

        return compounds

def assign_chirality(compounds):
    """
    Assigns chirality to all compounds in the compounds dictionary

    Parameters:
        compounds [dict of objects]
            -without chirality attribute populated

    Returns:
        compounds [dict of objects]
            -with chirality attribute populated with True OR False

    """
    for c in compounds:
        chiral = False
        
        for name in compounds[c].names:
            matchedDSLR = re.search( r'[DSLR]\W',name) # match any name with a D,S,L,R character followed by any non-alphanumeric character
            matchedDRNA = re.search( r'[DR]NA',name)  # match any name with DNA or RNA
            
            if matchedDSLR or matchedDRNA:
                chiral = True

        compounds[c].chirality = chiral

    return compounds


def main():

    filename = "data/compounds.txt"
    compounds = populate_compounds(filename)
    compounds = assign_chirality(compounds)

    n_chiral = 0
    for c in compounds:
        if compounds[c].chirality == True: n_chiral += 1

    with open("compound_testdict.dat", "wb") as myFile:
        pickle.dump(compounds, myFile)

    print n_chiral

if __name__ == '__main__':
    main()