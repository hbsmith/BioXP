import os

def dict_rxn_enz(kegg_reactions_file):
    input_file = open(kegg_reactions_file, 'r')
    list_rxn = []
    rxn_enz = {} #{RXN : [Enz1, Enz2, ... ]}
    newEnz = False
    for line in input_file:
        items = line.split()
        if items[0] == 'ENTRY':
            rxnName = items[1]
            list_rxn.append(rxnName)
            rxn_enz[rxnName] = []
        if items[0] == 'ENZYME':
            newEnz = True
            for s in items[1:]:
                rxn_enz[rxnName].append(s)
        if newEnz == True and items[0] != 'ENZYME':
            if line[0] == ' ':
                for s in items:
                    rxn_enz[rxnName].append(s)
            else:
                newEnz = False

    new_list_rxn = []
    new_rxn_enz = {}
    for x in list_rxn:
        if len(rxn_enz[x]) == 0:
            continue
        new_list_rxn.append(x)
        new_rxn_enz[x] = []
        for s in rxn_enz[x]:
            new_rxn_enz[x].append(s)

    return new_list_rxn, new_rxn_enz


def writing_rxn_enz(list_rxn, rxn_enz, output_file_name):
    output_file = open(output_file_name, 'w')
    output_file.write('# RXN \t ENZ')
    for x in list_rxn:
        output_file.write('\n%s'%(x))
        for s in rxn_enz[x]:
            output_file.write('\t%s'%(s))


def sort_enz_set(set_enz):
    unsorted_enz_list = []
    for s in set_enz:
        items = s.split('.')
        unsorted_enz_list.append(items)

    sorted_enz_list = sorted(unsorted_enz_list, key = lambda x: (x[0], x[1], x[2], x[3]))

    list_enz = []
    for x in sorted_enz_list:
        enz_name = '%s.%s.%s.%s'%(x[0], x[1], x[2], x[3])
        list_enz.append(enz_name)

    return list_enz


def dict_enz_rxn(rxn_enz):
    enz_rxn = {}
    set_enz = set()
    for x in rxn_enz.iterkeys():
        for s in rxn_enz[x]:
            if s not in set_enz:
                enz_rxn[s] = []
                set_enz.add(s)
            enz_rxn[s].append(x)

    list_enz = sort_enz_set(set_enz)

    return list_enz, enz_rxn

def writing_enz_rxn(list_enz, enz_rxn, output_file_name):
    EnzRxnFile = open(output_file_name, 'w')
    EnzRxnFile.write('# ENZ \t RXN')
    for s in list_enz:
        EnzRxnFile.write('\n%s'%(s))
        st_enz_rxn = sorted(enz_rxn[s])
        for x in st_enz_rxn:
            EnzRxnFile.write('\t%s'%(x))


def rxn_eq(kegg_reactions_file):
    input_file = open(kegg_reactions_file, 'r')
    list_rxn = []
    rxn_reactants = {} #{RXN : [Reac1, Reac2, ... ]}
    rxn_products = {} #{RXN : [Prod1, Prod2, ... ]}
    for line in input_file:
        items = line.split()
        if items[0] == 'ENTRY':
            rxnName = items[1]
            list_rxn.append(rxnName)
            rxn_reactants[rxnName] = []
            rxn_products[rxnName] = []

        if items[0] == 'EQUATION':
            items  = line.split('<=>')
            lefthand = items[0].split()
            righthand = items[1].split()

            elilgible_rxn = True
            for x in lefthand[1:]:
                if x[0] == 'G':
                    elilgible_rxn = False
            for x in righthand:
                if x[0] == 'G':
                    elilgible_rxn = False

            ##### list of reactants and list of products
            if elilgible_rxn:
                for x in lefthand[1:]: # [1:] ==> to exclude string EQUATION
                    if x[0] == 'C':
                        rxn_reactants[rxnName].append(x[0:6])
                for x in righthand:
                    if x[0] == 'C':
                        rxn_products[rxnName].append(x[0:6])
    return list_rxn, rxn_reactants, rxn_products

def writing_rxn_reactants(list_rxn, rxn_reactants, file_name):
    RxnReatFile = open(file_name, 'w')
    RxnReatFile.write('# Rxn \t Reactants')
    for x in list_rxn:
        if len(rxn_reactants[x]) == 0:
            continue
        RxnReatFile.write('\n%s'%(x))
        st_rxn_reactants = sorted(rxn_reactants[x])
        for w in st_rxn_reactants:
            RxnReatFile.write('\t%s'%(w))

def writing_rxn_products(list_rxn, rxn_products, file_name):
    RxnProdFile = open(file_name, 'w')
    RxnProdFile.write('# Rxn \t Products')
    for x in list_rxn:
        if len(rxn_products[x]) == 0:
            continue
        RxnProdFile.write('\n%s'%(x))
        st_rxn_products = sorted(rxn_products[x])
        for w in st_rxn_products:
            RxnProdFile.write('\t%s'%(w))


########################################################################

kegg_input_file = 'kegg/reactions.txt'

RxnEnz_file = 'data/kegg_rxn_enz.dat'
EnzRxn_file = 'data/kegg_enz_rxn.dat'
RxnReat_File = 'data/kegg_rxn_reactants.dat'
RxnProd_File = 'data/kegg_rxn_products.dat'


list_rxn, rxn_enz = dict_rxn_enz(kegg_input_file)
writing_rxn_enz(list_rxn, rxn_enz, RxnEnz_file)

list_enz, enz_rxn = dict_enz_rxn(rxn_enz)
writing_enz_rxn(list_enz, enz_rxn, EnzRxn_file)

list_rxn, rxn_reactants, rxn_products = rxn_eq(kegg_input_file)
writing_rxn_reactants(list_rxn, rxn_reactants, RxnReat_File)
writing_rxn_products(list_rxn, rxn_products, RxnProd_File)
