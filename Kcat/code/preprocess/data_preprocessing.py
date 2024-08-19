from zeep import Client
import hashlib
import requests
from urllib.request import urlopen, Request
import pandas as pd
import numpy as np
from os.path import join
import os
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
import pickle


def get_max_for_EC_number(EC):
    df = pd.DataFrame(columns = ["EC", "kcat VALUE"])
    df = add_kcat_for_EC_number(brenda_df = df, EC = EC)
    for ind in df.index:
        try:
            df["kcat VALUE"][ind] = float(df["kcat VALUE"][ind])
        except ValueError:
            df["kcat VALUE"][ind] = 0
    return(np.max(df["kcat VALUE"]))    


def mw_mets(metabolites):
    mw = 0
    for met in metabolites:
        if met != "":
            if met[0] == "C":
                try:
                    mol = Chem.MolFromMolFile(join("..", "..", "data", "metabolite_data",
                                           "mol-files", met + '.mol'))
                except:
                    mw = np.nan
                    break
            else:
                mol = Chem.inchi.MolFromInchi(met)
            mw = mw + Descriptors.MolWt(mol)
        
    return(mw)

def split_dataframe_enzyme(frac, df):
    df1 = pd.DataFrame(columns = list(df.columns))
    df2 = pd.DataFrame(columns = list(df.columns))
    
    #n_training_samples = int(cutoff * len(df))
    
    df.reset_index(inplace = True, drop = True)
    
    #frac = int(1/(1- cutoff))
    
    train_indices = []
    test_indices = []
    ind = 0
    while len(train_indices) +len(test_indices) < len(df):
        if ind not in train_indices and ind not in test_indices:
            if ind % frac != 0:
                n_old = len(train_indices)
                train_indices.append(ind)
                train_indices = list(set(train_indices))

                while n_old != len(train_indices):
                    n_old = len(train_indices)

                    training_seqs= list(set(df["Sequence"].loc[train_indices]))

                    train_indices = train_indices + (list(df.loc[df["Sequence"].isin(training_seqs)].index))
                    train_indices = list(set(train_indices))
                
            else:
                n_old = len(test_indices)
                test_indices.append(ind)
                test_indices = list(set(test_indices)) 

                while n_old != len(test_indices):
                    n_old = len(test_indices)

                    testing_seqs= list(set(df["Sequence"].loc[test_indices]))

                    test_indices = test_indices + (list(df.loc[df["Sequence"].isin(testing_seqs)].index))
                    test_indices = list(set(test_indices))
                
        ind +=1
    
    
    df1 = df.loc[train_indices]
    df2 = df.loc[test_indices]
    
    return(df1, df2)

def split_dataframe(frac, df):
    df1 = pd.DataFrame(columns = list(df.columns))
    df2 = pd.DataFrame(columns = list(df.columns))
    
    #n_training_samples = int(cutoff * len(df))
    
    df.reset_index(inplace = True, drop = True)
    
    #frac = int(1/(1- cutoff))
    
    train_indices = []
    test_indices = []
    ind = 0
    while len(train_indices) +len(test_indices) < len(df):
        if ind not in train_indices and ind not in test_indices:
            if ind % frac != 0:
                n_old = len(train_indices)
                train_indices.append(ind)
                train_indices = list(set(train_indices))

                while n_old != len(train_indices):
                    n_old = len(train_indices)

                    training_seqs= list(set(df["Sequence"].loc[train_indices]))
                    training_fps = list(set(df["structural_fp"].loc[train_indices]))

                    train_indices = train_indices + (list(df.loc[df["Sequence"].isin(training_seqs)].index) +
                                                           list(df.loc[df["structural_fp"].isin(training_fps)].index))
                    train_indices = list(set(train_indices))
                
            else:
                n_old = len(test_indices)
                test_indices.append(ind)
                test_indices = list(set(test_indices))

                while n_old != len(test_indices):
                    n_old = len(test_indices)

                    testing_seqs= list(set(df["Sequence"].loc[test_indices]))
                    testing_fps = list(set(df["structural_fp"].loc[test_indices]))

                    test_indices = test_indices + (list(df.loc[df["Sequence"].isin(testing_seqs)].index) +
                                                           list(df.loc[df["structural_fp"].isin(testing_fps)].index))
                    test_indices = list(set(test_indices))
                
        ind +=1
    
    
    df1 = df.loc[train_indices]
    df2 = df.loc[test_indices]
    
    return(df1, df2)
                