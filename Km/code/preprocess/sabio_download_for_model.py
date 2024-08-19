#!/usr/bin/python

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

def eclist():

    QUERY_URL = 'http://sabiork.h-its.org/sabioRestWebServices/kineticlawsExportTsv'

    query = {'fields[]':['ECNumber'], 'q':'Parametertype:"Km"'}

    request = requests.post(QUERY_URL, params = query)

    results = request.text

    ECs = [x for x in results.strip().split('\n')][2:]

    return sorted(list(set(ECs)))[2:]

def sabio_info(allEC):
    i = 0
    with open('../../Data/EC_Km_model_' + organism + '.tsv', 'w') as ECfile :
        records = ['ECs', 'Organism', 'Uniprot IDs', 'PMID', 'Type', 'Km', 'Temperature', 'pH', 'Substrates', 'Products', 'substrate_IDs', 'product_IDs', 'Main Substrate']
        ECfile.write('\t'.join(records) + '\n')

    with open('../../Data/min_EC_' + organism + '.tsv', 'w') as file :
        file.write('\t'.join(['EC', "min_Km"]) + '\n')

    ids = {}
    for EC in allEC[allEC.index('1.14.14.43'):]:
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
            min_Km=list()
            for line in results.strip().split('\n')[1:]:
                substrateids = []
                productids = []
                complete = True
                entry = line.split('\t')
                
                if not is_seed_plant(entry[0]):
                    continue
                
                if entry[8] != "Km" or not entry[10]:
                    continue

                if entry[13] != 'M':
                    continue

                min_Km.append(float(entry[10]))

                print(line)

                substrates=entry[1]
                products=entry[2]
                main_substrate = "None"

                if entry[9] in substrates:
                    main_substrate = entry[9]
                else:
                    continue

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

                    for s in substrates.split(';'):
                        substrateids.append(ids[s])

                    for p in products.split(';'):
                        productids.append(ids[p])

                    if entry[3]:
                        with open('../../Data/EC_Km_model_' + organism + '.tsv', 'a') as ECfile :
                            ECfile.write('\t'.join([EC, entry[0], ';'.join(list(set(entry[3].split(' ')))), entry[5], entry[4], entry[10], entry[7], entry[6],
                             substrates, products, '#'.join(substrateids), '#'.join(productids), ids[main_substrate]]) + '\n')
                        
            if min_Km:
                print(min_Km)
                with open('../../Data/min_EC_' + organism + '.tsv', 'a') as ECfile :
                        ECfile.write('\t'.join([EC, str(min(min_Km))]) + '\n')

if __name__ == '__main__' :
    allEC = eclist()
    print(allEC)
    sabio_info(allEC)


