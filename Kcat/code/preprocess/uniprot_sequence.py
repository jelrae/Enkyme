#!/usr/bin/python
# coding: utf-8

from urllib import request
import requests
import re
from zeep import Client
import hashlib

organism = "Seed plants"

def uniprot_sequence(id) :
    url = "https://www.uniprot.org/uniprot/%s.fasta" % id
    IdSeq = dict()
    try :
        data = request.urlopen(url)
        respdata = data.read().decode("utf-8").strip()
        IdSeq[id] =  "".join(respdata.split("\n")[1:])
    except :
        print(id, "Cannot find from uniprot!")
        IdSeq[id] = None
    print(IdSeq[id])
    return IdSeq[id]

def seq_by_ec_organism(ec, organism) :
    IdSeq = dict()
    params = {"query": "ec:%s AND organism:%s AND reviewed:yes" % (ec, organism), "format": "fasta"}
    response = requests.get("http://www.uniprot.org/uniprot/", params=params)

    try :
        respdata = response.text
        seq = dict()
        for line in respdata.split('\n') :
            if line.startswith('>') :
                name=line
                seq[name] = ''
            else :
                seq[name] += line.replace('\n', '').strip()
        IdSeq[ec+'&'+organism] =  list(seq.values())

    except :
        print(ec+'&'+organism, "Cannot find from uniprot!")
        IdSeq[ec+'&'+organism] = None

    print(IdSeq[ec+'&'+organism])
    return IdSeq[ec+'&'+organism]

def seq_by_brenda(ec, organism) :
    email = ''
    password = ''

    wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
    password = hashlib.sha256(password.encode("utf-8")).hexdigest()
    client = Client(wsdl)

    parameters = (email,password,"ecNumber*%s" % ec,"organism*%s" % organism, "sequence*", "noOfAminoAcids*", "firstAccessionCode*", "source*Swiss-Prot", "id*" )
    entries = client.service.getSequence(*parameters)

    sequences = list()
    if entries :
        for entry in entries :
            sequences.append(entry['sequence'])

    return sequences
    
def main() :
    with open('../../data/kcat_model_' + organism + '.tsv', 'w') as outfile :
        records = ['ECs', 'Organism', 'Uniprot IDs', 'PMID', 'Type', 'kcat', 'Temperature', 'pH', 'Substrates', 'Products' , 'substrate_IDs', 'product_IDs', 'Main Substrate', 'Sequence']
        outfile.write('\t'.join(records) + '\n')

    with open("../../data/EC_kcat_model_" + organism + ".tsv", "r", encoding='utf-8') as file :
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
                with open('../../data/kcat_model_' + organism + '.tsv', 'a') as outfile :
                    outfile.write('\t'.join([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], ';'.join(seq)]) + '\n')
            else:
                with open('../../data/kcat_model_' + organism + '.tsv', 'a') as outfile :
                    seq = uniprot_sequence(data[2])
                    if 'wildtype' in data[4] and seq:
                        with open('../../data/kcat_model_' + organism + '.tsv', 'a') as outfile :
                            outfile.write('\t'.join([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], seq]) + '\n')
                    else:
                        try:
                            mutantSites = re.findall('[A-Z]\d+[A-Z]', data[4])
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
                                with open('../../data/kcat_model_' + organism + '.tsv', 'a') as outfile :
                                    outfile.write('\t'.join([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], mutatedSeq]) + '\n')   

                        except:
                            continue               

        else:
            seq = seq_by_ec_organism(data[0], data[1])
            if seq:
                with open('../../data/kcat_model_' + organism + '.tsv', 'a') as outfile :
                    outfile.write('\t'.join([data[0], data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10], data[11], data[12], ';'.join(seq)]) + '\n')
        #     else:
        #         seq = seq_by_brenda(data[4], organism)
        #         if seq:
        #             with open('../../data/kcat_model_' + organism + '.tsv', 'a') as outfile :
        #                 outfile.write('\t'.join(['', data[1], data[2], data[3], data[4], data[5], data[6], ';'.join(seq), data[7]]) + '\n')


if __name__ == "__main__" :
    main()


