import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from brendapyrser import BRENDA
from itertools import filterfalse
import sys

def get_enzyme_name(reaction):
    return reaction.name

def get_substrates_and_products(reaction):
    rxnstr = reaction.reaction_str

    if rxnstr == '':
        SP = reaction.substratesAndProducts
        if len(SP) > 1:
            reactants = SP[0]['substrates']
            products = SP[0]['products']
        if len(SP) == 1:
            reactants = SP['substrates']
            products = SP['products']

    else:
        splitrxn = rxnstr.split('<=> ', 2)
        reactants = splitrxn[0].split(' + ')
        products = splitrxn[1].split(' + ')

    return(reactants, products)

def get_cofactors(reaction, species):
    cofactor_list = []

    for cofactor in reaction.cofactors:
        if species in reaction.cofactors[cofactor]['species']:
            cofactor_list.append(cofactor)

    return cofactor_list

def get_parameters(rxn, species):

    parameters_list = []

    # gets list of all KM and Kcat values for the organism
    kmvals = rxn.KMvalues.filter_by_organism(species)
    kcatvals = rxn.Kcatvalues.filter_by_organism(species)

    # pulls the substrates from both lists and combines them with no duplicates
    substrate_list_km = []
    for i in kmvals:
        substrate_list_km.append(i)

    substrate_list_kcat = []
    for i in kcatvals:
        substrate_list_kcat.append(i)

    substrate_list = list(set(substrate_list_km + substrate_list_kcat))

    for compound in substrate_list:
        KMfiltered = kmvals.filter_by_compound(compound).get_values()
        KMMfiltered = rxn.KKMvalues.filter_by_organism(species).filter_by_compound(compound).get_values()
        Kcatfiltered = kcatvals.filter_by_compound(compound).get_values()
        KMfiltered = list(filterfalse(KMMfiltered.__contains__, KMfiltered))

        KMvals_WT = []

        if len(KMfiltered) > 1:
            kval = kmvals.filter_by_compound(compound)

            metalist = []
            for k in kval.keys():
                for v in kval[k]:
                    metalist.append(v['meta'])

            for idx, meta in enumerate(metalist):
                if "wild" in meta:
                    KMvals_WT.append(kval[k][idx]['value'])

        Kcatvals_WT = []
        if len(Kcatfiltered) > 1:
            kval = kcatvals.filter_by_compound(compound)

            metalist = []
            for k in kval.keys():
                for v in kval[k]:
                    metalist.append(v['meta'])

            for idx, meta in enumerate(metalist):
                if "wild" in meta:
                    Kcatvals_WT.append(kval[k][idx]['value'])

        if len(KMvals_WT) != 0:
            KMfiltered = KMvals_WT

        if len(Kcatvals_WT) != 0:
            Kcatfiltered = Kcatvals_WT

        KM = sum(KMfiltered)/(len(KMfiltered) or 1)
        KMM = sum(KMMfiltered)/(len(KMMfiltered) or 1)
        Kcat = sum(Kcatfiltered)/(len(Kcatfiltered) or 1)

        parameters_list.append([compound, KM, Kcat])

    return parameters_list

def get_metals(reaction, species):

    LoM = []

    for metal in reaction.metals.keys():
        if species in reaction.metals[metal]['species']:
            LoM.append(metal)

    return LoM

def parser(rxn, species):

    name = get_enzyme_name(rxn)

    substrates, products = get_substrates_and_products(rxn)

    cofactors = get_cofactors(rxn, species)

    parameters = get_parameters(rxn,species)

    metals = get_metals(rxn,species)

    return(name, substrates, products, cofactors, parameters, metals)

def isNaN(string):
    return string != string

def main():
    filename = sys.argv[1]
    dataFile = sys.argv[2]

    print('this is running')

    df = pd.read_csv(filename)
    brenda = BRENDA(dataFile)

    for index, row in df.iterrows():
        ID = row['EC']
        species = row['Species']

        if ID == '-':
            print('skipping this ID')
            continue
        else:
            rxn = brenda.reactions.get_by_id(ID)

            print('now processing ', ID)

            if isNaN(species):
                species = 'Escherichia coli'

            result_df = pd.DataFrame(parser(rxn, species)).T
            values = result_df.values[0]

            df.iloc[index, 3] = values[0]

            subs = '; '.join(values[1])
            df.iloc[index, 4] = subs

            prods = '; '.join(values[2])
            df.iloc[index, 5] = prods

            cof = '; '.join(values[3])
            df.iloc[index, 6] = cof

            par = '; '.join([f'{item[0]}_Km: {item[1]}; {item[0]}_Kcat: {item[2]}; f' for item in values[4]])
            df.iloc[index, 7] = par

            met = '; '.join(values[5])
            df.iloc[index, 8] = met

    df.to_csv(filename, index=False)

if __name__ == '__main__':
    main()
