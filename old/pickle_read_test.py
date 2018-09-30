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

with open("compound_testdict.dat", "rb") as myFile:
    myNewPulledInDictionary = pickle.load(myFile)

print myNewPulledInDictionary["C00002"].names