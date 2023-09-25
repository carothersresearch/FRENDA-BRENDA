import numpy as np
import pandas as pd
import math
from matplotlib import pyplot as plt
from brendapyrser import BRENDA
from itertools import filterfalse
import sys
import re

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

def replaceCofactors(common_cofactor, cof_product, substrates, products, cofactors):
    if common_cofactor in substrates:
    # add to cofactors
        if common_cofactor in cofactors:
            pass
        else:
            cofactors.append(common_cofactor)

        # add to products
        if cof_product in products:
            pass
        else:
            products.append(cof_product)

        # remove from substrates
        substrates.remove(common_cofactor)

    return substrates, products, cofactors

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

def reorderParams(parameters, substrates):
    listParams = []

    for i, item in enumerate(parameters):
        listParams.append(item[0])

    newSubsList = [None] * len(listParams)

    for substrate in substrates:
        if substrate in listParams:
            newSubsList[listParams.index(substrate)] = substrate
        else:
            pass

    return newSubsList

def get_metals(reaction, species):

    LoM = []

    for metal in reaction.metals.keys():
        if species in reaction.metals[metal]['species']:
            LoM.append(metal)

    return LoM

def average_duplicates(filtered_parameters):
    # Create a dictionary to store cumulative values for each unique element
    cumulative_values = {}

    # Iterate through the list and calculate the cumulative sum and count for each element
    for element, val1, val2 in filtered_parameters:
        if element not in cumulative_values:
            cumulative_values[element] = [0, 0, 0]  # [sum of val1, sum of val2, count of occurrences]

        # Check if the value is non-zero before adding to the cumulative sum and count
        if val1 != 0:
            cumulative_values[element][0] += val1
            cumulative_values[element][2] += 1
        if val2 != 0:
            cumulative_values[element][1] += val2

    # Calculate the averages for each element and update the list
    for element in cumulative_values:
        if cumulative_values[element][2] != 0:
            cumulative_values[element][0] /= cumulative_values[element][2]
        if cumulative_values[element][1] != 0:
            cumulative_values[element][1] /= cumulative_values[element][2]

    # Convert the dictionary back to a list
    result = [[element, val1, val2] for element, (val1, val2, _) in cumulative_values.items()]

    return result


def parser_unfilt(rxn, species):

    name = get_enzyme_name(rxn)

    substrates, products = get_substrates_and_products(rxn)

    cofactors = get_cofactors(rxn, species)

    parameters = get_parameters(rxn,species)

    metals = get_metals(rxn,species)

    return(name, substrates, products, cofactors, parameters, metals)


def parser_filt(rxn, species):

    name = get_enzyme_name(rxn)

    substrates, products = get_substrates_and_products(rxn)

    for i, substrate in enumerate(substrates):
        substrate = substrate.replace(' ','')
        substrate = substrate.replace('D-','')
        substrate = substrate.lower()
        substrate = substrate.replace('-','')
        substrate = substrate.replace('alpha','')
        substrates[i] = substrate

    for i, product in enumerate(products):
        product = product.replace(' ','')
        product = product.replace('D-','')
        product = product.lower()
        product = product.replace('-','')
        product = product.replace('alpha','')
        products[i] = product

    cofactors = get_cofactors(rxn, species)

    for i, item in enumerate(cofactors):
        item = item.replace(' ','')
        item = item.replace('D-','')
        item = item.lower()
        item = item.replace('-','')
        item = item.replace('alpha','')
        cofactors[i] = item

    # I think we could easily turn this into a function with a list of generalizable cofactors
    # I am not sure if I have all the cofactors! We might want the inverse of this too? ie. ADP to ATP?
    substrates, products, cofactors = replaceCofactors('atp', 'adp', substrates, products, cofactors)
    substrates, products, cofactors = replaceCofactors('nadh', 'nad+', substrates, products, cofactors)
    substrates, products, cofactors = replaceCofactors('nadph', 'nadp+', substrates, products, cofactors)
    substrates, products, cofactors = replaceCofactors('fadh', 'fad', substrates, products, cofactors)

    parameters = get_parameters(rxn,species)

    for listitem in parameters:
        string = listitem[0]
        string = string.replace(' ','')
        string = string.replace('D-','')
        string = string.lower()
        string = string.replace('-','')
        string = string.replace('alpha','')
        listitem[0] = string

    filtered_parameters = []
    for p in parameters:
        if(p[0] in substrates):
            filtered_parameters.append(p)

    av_filt_par = average_duplicates(filtered_parameters)

    substrates = reorderParams(av_filt_par, substrates)

    metals = get_metals(rxn,species)

    return(name, substrates, products, cofactors, av_filt_par, metals)

def isNaN(string):
    return string != string

def rmBlanks(df1, df2):
    listrm = []
    for index, row in df1.iterrows():
        ID = row['EC']
        if ID == '-':
            listrm.append(index)
    df1 = df1.drop(listrm)
    df1.reset_index(inplace=True, drop=True)
    df2 = df2.drop(listrm)
    df2.reset_index(inplace=True, drop=True)
    return df1, df2

def clean_text(text):
    if isinstance(text, str):  # Check if the value is a string
        text = re.sub(r"\([^)]*\)", "", text)  # Remove text within parentheses
        text = re.sub(r"'", "", text)  # Remove single quotes
        text = text.replace("(", "").replace(")", "")  # Remove lone parentheses
        text = re.sub(r"\[", "", text)  # Remove single quotes
        text = re.sub(r"]", "", text)  # Remove single quotes
        return text
    else:
        return text

def main():
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    dataFile = sys.argv[3]
    filtoption = sys.argv[4]

    print('this is running')

    df = pd.read_csv(filename1)
    df2 = pd.read_csv(filename2)

    df, df2 = rmBlanks(df, df2)

    brenda = BRENDA(dataFile)

    for index, row in df.iterrows():
        ID = row['EC']
        species = row['Species']

        if ID == '-':
            print('skipping this ID')
            #df = df.drop(index)
            continue
        else:
            try:
                rxn = brenda.reactions.get_by_id(ID)
            except ValueError as error:
                print(error)
                continue

            print('now processing ', ID)

            if isNaN(species):
                species = 'Escherichia coli'

            if filtoption == 'f':
                try:
                    result_df = pd.DataFrame(parser_filt(rxn, species)).T
                except Exception as error:
                    print(error)
                    print('this threw and error and should break now')
                    continue
            else:
                result_df = pd.DataFrame(parser_unfilt(rxn, species)).T

            values = result_df.values[0]

            print('the index is', index)
            print('the ID is ', ID)
            print('the enzyme name is ', values[0])

            # fourth column, adding a LABEL for now, could get rid of later if we change odbm_main ModelBuilder class
            df.iloc[index, 3] = f'R{index+1}'

            # fifth column, enzyme name
            df.iloc[index, 4] = values[0]

            # seventh column, substrates
            subs = '; '.join(values[1])
            df.iloc[index, 6] = subs

            # eigth column, cofactors
            cof = '; '.join(values[3])
            df.iloc[index, 7] = cof

            # ninth column, products
            prods = '; '.join(values[2])
            df.iloc[index, 8] = prods

            # tenth column, metals
            met = '; '.join(values[5])
            df.iloc[index, 9] = met

            if filtoption == 'f':
                if len(values[4]) == 1:
                    mech = 'MM'
                    par = '; '.join([f'kcat:{item[2]};Km:{item[1]}' for item in values[4]])

                elif len(values[4]) == 2:
                    mech = 'SOBB'

                    KM1 = values[4][0][1]
                    KM2 = values[4][1][1]
                    kcat1 = values[4][0][2]
                    kcat2= values[4][1][2]

                    if kcat1 == 0 and kcat2 == 0:
                        kcat = 0
                    elif math.isclose(kcat1, kcat2):
                        kcat = kcat1
                    elif kcat1 == 0:
                        kcat = kcat2
                    elif kcat2 == 0:
                        kcat = kcat1
                    else:
                        kcat = (kcat1 + kcat2)/2

                    par = '; '.join([f'kcat:{kcat};Km1:{KM1};Km2:{KM2}'])

                else:
                    mech = ''
                    par = '; '.join([f'{item[0]}_Km: {item[1]}; {item[0]}_Kcat: {item[2]}; ' for item in values[4]])

                df.iloc[index, 10] = par
            else:
                # eleventh column, parameters
                par = '; '.join([f'{item[0]}_Km: {item[1]}; {item[0]}_Kcat: {item[2]}; ' for item in values[4]])
                df.iloc[index, 10] = par

            # sixth column, this is mechanism, just to get it running for now. will need to fix later
            df.iloc[index, 5] = mech

            # this will be identical order to the other df, but it gets the name in Reaction.csv
            df2.iloc[index,0] = values[0]

    # Apply the functions to the DataFrame, but only to the string columns
    string_columns = df.select_dtypes(include=['object']).columns  # Select only string columns
    for column in string_columns:
        df[column] = df[column].apply(clean_text)

    # Apply the functions to the DataFrame, but only to the string columns
    string_columns = df2.select_dtypes(include=['object']).columns  # Select only string columns
    for column in string_columns:
        df2[column] = df2[column].apply(clean_text)

    df.to_csv(filename1, index=False)
    df2.to_csv(filename2, index=False)

if __name__ == '__main__':
    main()
