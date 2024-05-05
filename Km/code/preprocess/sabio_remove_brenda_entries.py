#!/usr/bin/python
################################################################################
# retrieveBRENDA
# Acces the web client and retrieves all EC data from BRENDA. Creates files with
# BRENDA output for all organisms and EC numbers for which there is data.
#
# Benjamin Sanchez. Last edited: 2018-04-10
################################################################################

# Updated by:
# Author: LE YUAN
# This code should be run under the Python 2.7 environment


#INPUTS:
#1) Path in which you wish to store all BRENDA queries:
# output_path = '/Users/.../brenda_parser/raw_data'
# output_path = '../../Data/database/brenda_ec'
output_path = '../../Data/'
#2) Last field processed (if the program was interrupted), e.g. 'KM'. If you
#   want to start from scratch, leave empty:
last_field = ''
#3) Last EC number processed (if the program was interrupted), e.g. '1.2.3.4'.
#   If you want to start from scratch, leave empty:
last_EC = ''
#4) E-mail in BRENDA:
email = 'm.van.laar@student.vu.nl'
#5) Password in BRENDA:
password = 'BRENDAenzyming'

################################################################################

def get_PMIDs_that_are_in_BRENDA(Sabio_PMIDs, df, organism):
    PMIDs_in_BRENDA = []

    for PMID in Sabio_PMIDs:

        if math.isnan(PMID):
            continue
        help_df = df.loc[df["PMID"] == PMID]
        EC = list(help_df["ECs"])[0]

        parameters = (credentials + ",ecNumber*" + EC + "#organism*" + organism + "#pubmedId*" + str(int(PMID)))
        print(parameters)
        resultString1 = client.getReference(parameters)   
        print(resultString1.decode('ascii','ignore'))

        if len(resultString1) > 0:
            resultString = resultString1.split('!')
            for string in resultString:
                reference = string.split('#')[1].split('*')[1]
                organism = string.split('#')[8].split('*')[1]
                ecNumber = string.split('#')[0].split('*')[1]
                print(reference, organism, ecNumber)

                parameters = (credentials + ",ecNumber*" + ecNumber + "#organism*" + organism + "#literature*" + reference)
                resultString2 = client.getKmValue(parameters)
                print(resultString2.decode('ascii','ignore'))
                if len(resultString2) > 0:
                    PMIDs_in_BRENDA.append(PMID)
                    break
                sleep(0.1)
        sleep(0.1)
        
    return(PMIDs_in_BRENDA)

#Main script
                    
#Change path:
import os
import math
import csv
from time import sleep
prev_path = os.getcwd()
os.chdir(output_path)

#Construct BRENDA client:
import string
import hashlib
import pandas as pd
from SOAPpy import SOAPProxy ## for usage without WSDL file
endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
client      = SOAPProxy(endpointURL)
password    = hashlib.sha256(password).hexdigest()
credentials = email + ',' + password
organism = "Arabidopsis thaliana"

df = pd.read_csv('KM_model_AT.tsv', sep='\t')

Sabio_PMIDs = list(set(df["PMID"]))
PMIDs_in_BRENDA = get_PMIDs_that_are_in_BRENDA(Sabio_PMIDs, df, organism)

droplist = []

for ind in df.index:
    if df["PMID"][ind] in PMIDs_in_BRENDA:
        droplist.append(ind)
        
df.drop(droplist, inplace = True)
sabio_df = df
#sabio_df = sabio_df.groupby(["Sequence", "substrate ID", "ECs"], as_index = False)["KM"].mean()

sabio_df.to_csv('KM_model_AT.tsv', sep='\t', index=False)


################################################################################
