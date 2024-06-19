#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-07-10fabacea

import requests
from time import sleep
from ete3 import NCBITaxa
ncbi = NCBITaxa()

def is_seed_plant(org):
    try:
        tax_id = ncbi.get_name_translator([org])[org][0]
        lineage = ncbi.get_lineage(tax_id)
        if 58024 not in lineage:
            return(False)
        else:
            return(True)
    except KeyError:
        return(False)

organism = "Seed plants"

# Extract EC number list from ExPASy, which is a repository of information relative to the nomenclature of enzymes.
def eclist():
    # with open('../../Data/enzyme.dat', 'r') as outfile :
    #     lines = outfile.readlines()

    # ec_list = list()
    # for line in lines :
    #     if line.startswith('ID') :
    #         ec = line.strip().split('  ')[1]
    #         ec_list.append(ec)
    # # print(ec_list)
    # print(len(ec_list)) # 7906
    # return ec_list

    QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv'

    query = {'fields[]':['ECNumber'], 'q':'Parametertype:"kcat"'}

    request = requests.post(QUERY_URL, params = query)

    results = request.text

    ECs = [x for x in results.strip().split('\n')][1:]

    return sorted(list(set(ECs)))[1:]

def sabio_info(allEC):
    i = 0
    with open('../../Data/EC_kcat_model_' + organism + '.tsv', 'w') as ECfile :
        records = ['ECs', 'Organism', 'Uniprot IDs', 'PMID', 'Type', 'kcat', 'Temperature', 'pH', 'Substrates', 'Products', 'substrate_IDs', 'product_IDs', 'Main Substrate']
        ECfile.write('\t'.join(records) + '\n')

    with open('../../Data/max_EC_' + organism + '.tsv', 'w') as file :
        file.write('\t'.join(['EC', "max_kcat"]) + '\n')

    ids = {}
    for EC in allEC:
        QUERY_URL = 'https://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv'
        i += 1
        print('This is %s ----------------------------' %EC)
        query_dict = {"ECNumber":'%s' %EC,"Organism":'*'}
        query_string = ' AND '.join(['%s:%s' % (k,v) for k,v in query_dict.items()])

        query = {'fields[]':['Organism', 'Substrate', 'Product', 'UniprotID', 'EnzymeType', 'PubMedID', 'pH', 'Temperature', 'Parameter'], 'format':'tsv', 'q':query_string}
        
        try:
            request = requests.post(QUERY_URL, params = query)
        except requests.exceptions.SSLError:
            sleep(5)
            try:
                request = requests.post(QUERY_URL, params = query)
            except requests.exceptions.SSLError:
                continue

        results = request.text
        print('---------------------------------------------')

        if results:
            max_kcat=list()
            for line in results.strip().split('\n')[1:]:
                substrateids = []
                productids = []
                complete = True
                entry = line.split('\t')
                
                if not is_seed_plant(entry[0]):
                    continue

                if entry[8] != "kcat" or not entry[10]:
                    continue
                
                if entry[13] == 'm^(-1)':
                    entry[10] = str(float(entry[10])/60)
                    entry[13] = 's^(-1)'

                if entry[13] != 's^(-1)':
                    continue

                max_kcat.append(float(entry[10]))

                substrates=entry[1]
                products=entry[2]

                main_substrate = 'None'

                count=0
                for line2 in results.strip().split('\n')[1:]:
                    entry2 = line2.split('\t')
                    if entry2[8] == 'kcat/Km' and entry2[:7] == entry[:7]:
                        if entry[9] != main_substrate:
                                count += 1
                        main_substrate = entry2[9]
                if count == 0:
                    for line2 in results.strip().split('\n')[1:]:
                        entry2 = line2.split('\t')
                        if entry2[8] == 'Km' and entry2[:7] == entry[:7]:
                            if entry[9] != main_substrate:
                                count += 1
                            main_substrate = entry2[9]
                if count != 1 or main_substrate not in substrates:
                    complete=False

                for c in (substrates+';'+products).split(';'):
                    compound = ''
                    compoundid = ''
                    compoundinfo = ''
                    if c not in ids.keys():
                        QUERY_URL = 'https://sabiork.h-its.org/sabioRestWebServices/searchCompoundSynonyms'
                        query = {'fields[]':["SabioCompoundID"], "CompoundName":c}
                        try:
                            compound = (requests.get(QUERY_URL, params = query)).text
                        except requests.exceptions.SSLError:
                            sleep(5)
                            try:
                                compound = (requests.get(QUERY_URL, params = query)).text
                            except requests.exceptions.SSLError:
                                complete=False
                        try:
                            compoundid = compound.strip().split('\n')[1]
                        except:
                            complete = False

                        QUERY_URL = 'https://sabiork.h-its.org/sabioRestWebServices/searchCompoundDetails'
                        query = {'fields[]':["KeggCompoundID", "InChI"], "SabioCompoundID":compoundid}
                        try:
                            compoundinfo = (requests.get(QUERY_URL, params = query)).text
                        except requests.exceptions.SSLError:
                            sleep(5)
                            try:
                                compoundinfo = (requests.get(QUERY_URL, params = query)).text     
                            except requests.exceptions.SSLError:
                                complete=False                     
                        try:  
                            compoundreps = compoundinfo.strip().split('\n')[1]
                            compoundrep = compoundreps.split('\t')
                            if compoundrep[1] and compoundrep[1] != 'null':
                                ids[c] = compoundrep[1]
                            elif compoundrep[0] and compoundrep[0] != 'null':
                                ids[c] = compoundrep[0]
                            else:
                                complete=False
                        except:
                            complete = False

                if (complete == False):
                    continue
                else:
                    print(line)
                    
                    for s in substrates.split(';'):
                        substrateids.append(ids[s])

                    for p in products.split(';'):
                        productids.append(ids[p])

                    substrateids = sorted(substrateids)
                    productids = sorted(productids)

                    # substrate=''
                    # product=''

                    # for s in substrateids:
                    #     substrate += list(ids.keys())[list(ids.values()).index(s)] + ';'

                    # for p in productids:
                    #     product += list(ids.keys())[list(ids.values()).index(p)] + ';'

                    if entry[3]:
                        with open('../../Data/EC_kcat_model_' + organism + '2.tsv', 'a') as ECfile :
                            ECfile.write('\t'.join([EC, entry[0], ';'.join(list(set(entry[3].split(' ')))), entry[5], entry[4], entry[10], entry[7], entry[6],
                             substrates, products, '#'.join(substrateids), '#'.join(productids), ids[main_substrate]]) + '\n')
                        
            if max_kcat:
                print(max_kcat)
                with open('../../Data/max_EC_' + organism + '2.tsv', 'a') as ECfile :
                        ECfile.write('\t'.join([EC, str(max(max_kcat))]) + '\n')

if __name__ == '__main__' :
    allEC = eclist()
    print(allEC)
    sabio_info(allEC)


