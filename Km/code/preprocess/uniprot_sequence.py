#!/usr/bin/python
# coding: utf-8

# Author: LE YUAN
# Date: 2020-06-25

# This python script is to obtain protein sequence by uniprot protein id

from urllib import request
import requests
import csv
import re
from zeep import Client
import hashlib

organism = "Seed plants"

# This function is to obtain the protein sequence according to the protein id from Uniprot API
# https://www.uniprot.org/uniprot/A0A1D8PIP5.fasta 
# https://www.uniprot.org/help/api_idmapping
def uniprot_sequence(id) :
    url = "https://www.uniprot.org/uniprot/%s.fasta" % id
    IdSeq = dict()

    try :
        data = request.urlopen(url)
        respdata = data.read().decode("utf-8").strip()
        IdSeq[id] =  "".join(respdata.split("\n")[1:])
    except :
        print(id, "can not find from uniprot!")
        IdSeq[id] = None
    print(IdSeq[id])
    return IdSeq[id]

def seq_by_ec_organism(ec, organism) :
    IdSeq = dict()
    # https://www.biostars.org/p/356687/
    params = {"query": "ec:%s AND organism:%s AND reviewed:yes" % (ec, organism), "format": "fasta"}
    response = requests.get("http://www.uniprot.org/uniprot/", params=params)
    # print(type(response.text)) # <class 'str'>

    try :
        # respdata = response.text.strip()
        # # print(respdata) 
        # IdSeq[ec+'&'+organism] =  "".join(respdata.split("\n")[1:])

        respdata = response.text
        # print(respdata)
        sequence = list()
        seq = dict()
        i = 0
        for line in respdata.split('\n') :
            if line.startswith('>') :
                name=line
                seq[name] = ''
            else :
                seq[name] += line.replace('\n', '').strip()
        IdSeq[ec+'&'+organism] =  list(seq.values())

    except :
        print(ec+'&'+organism, "can not find from uniprot!")
        IdSeq[ec+'&'+organism] = None

    print(IdSeq[ec+'&'+organism])
    return IdSeq[ec+'&'+organism]

def seq_by_brenda(ec, organism) :
    # E-mail in BRENDA:
    email = 'm.van.laar@student.vu.nl'
    # Password in BRENDA:
    password = 'BRENDAenzyming'

    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    password    = hashlib.sha256(password.encode("utf-8")).hexdigest()
    client = Client(wsdl)
    # credentials = email + ',' + password

    # parameters = credentials+","+"ecNumber*%s#organism*%s" %(ec, organism)
    parameters = ( email,password,"ecNumber*%s" % ec,"organism*%s" % organism, "sequence*", "noOfAminoAcids*", "firstAccessionCode*", "source*Swiss-Prot", "id*" ) # *Swiss-Prot
    entries = client.service.getSequence(*parameters)

    sequences = list()
    # print(split_sequences)
    if entries :
        for entry in entries :
            sequences.append(entry['sequence'])

    return sequences
    
def main() :
    with open('../../Data/Km_model_' + organism + '.tsv', 'w') as outfile :
        records = ['ECs', 'Organism', 'Uniprot IDs', 'PMID', 'Type', 'Km', 'Temperature', 'pH', 'Substrates', 'Products' , 'substrate_IDs', 'product_IDs', 'Main Substrate', 'Sequence']
        outfile.write('\t'.join(records) + '\n')

    with open("../../Data/EC_Km_model_" + organism + ".tsv", "r", encoding='utf-8') as file :
        lines = file.readlines()[1:]

    for line in lines :
        if (line=="\n"):
            continue
        data = line.strip().split('\t')
        if data[2]:
            if ';' in data[2]:
                seq = list()
                for id in data[2].split(';'):
                    seq.append(uniprot_sequence(id))
                if seq:
                    with open('../../Data/Km_model_' + organism + '.tsv', 'a') as outfile :
                        outfile.write('\t'.join([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], ';'.join(seq)]) + '\n')

            else:
                with open('../../Data/Km_model_' + organism + '.tsv', 'a') as outfile :
                    seq = uniprot_sequence(data[2])
                    if 'wildtype' in data[4] and seq:
                        with open('../../Data/Km_model_' + organism + '.tsv', 'a') as outfile :
                            outfile.write('\t'.join([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], seq]) + '\n')
                    else:
                        try:
                            mutantSites = re.findall('[A-Z]\d+[A-Z]', data[4])  # re is of great use
                            mutant1_1 = [mutantSite[1:-1] for mutantSite in mutantSites]
                            mutant1_2 = [mutantSite for mutantSite in mutantSites]
                            mutant1 = [mutant1_1, mutant1_2]
                            mutant2 = set(mutant1[0])
                            if len(mutant1[0]) != len(mutant2) :
                                print(mutant1)
                                n += 1
                                print(str(n) + '---------------------------')

                            mutatedSeq = seq
                            for mutantSite in mutantSites :
                                if mutatedSeq[int(mutantSite[1:-1])-1] == mutantSite[0] :
                                    mutatedSeq = list(mutatedSeq)
                                    mutatedSeq[int(mutantSite[1:-1])-1] = mutantSite[-1]
                                    mutatedSeq = ''.join(mutatedSeq)
                                    if not mutatedSeq :
                                        print('-------------')
                                else :
                                    mutatedSeq = ''
                            
                            if mutatedSeq:
                                with open('../../Data/Km_model_' + organism + '.tsv', 'a') as outfile :
                                    outfile.write('\t'.join([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], mutatedSeq]) + '\n')

                        except:
                            continue                    

        else:
            seq = seq_by_ec_organism(data[0], data[1])
            if seq:
                with open('../../Data/Km_model_' + organism + '.tsv', 'a') as outfile :
                    outfile.write('\t'.join([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], ';'.join(seq)]) + '\n')
        #     else:
        #         seq = seq_by_brenda(data[4], organism)
        #         if seq:
        #             with open('../../Data/Km_model_' + organism + '.tsv', 'a') as outfile :
        #                 outfile.write('\t'.join(['', data[1], data[2], data[3], data[4], data[5], data[6], ';'.join(seq), data[7]]) + '\n')


if __name__ == "__main__" :
    main()


