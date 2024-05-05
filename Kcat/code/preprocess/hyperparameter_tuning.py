import os
import numpy as np
from os.path import join
import pandas as pd
from build_GNN import *
from tqdm import tqdm
import multiprocessing
import random
random.seed(10)
np.random.seed(84)
tf.random.set_seed(84)
tf.config.threading.set_intra_op_parallelism_threads(84)
tf.random.set_seed(84)

organism = "Seed plants"
datasets_dir = "../../data"
train_df = pd.read_pickle(join(datasets_dir, "train_df_kcat_%s.pkl" %organism))
CV_train_indices = np.load(join(datasets_dir, "CV_train_indices_%s" %organism))
CV_test_indices = np.load(join(datasets_dir, "CV_test_indices_%s" %organism))

train_indices = os.listdir(join(datasets_dir, "GNN_input_data"))
train_indices = [index[:index.rfind("_")] for index in train_indices]
train_indices = list(set([index for index in train_indices if "train" in index]))

test_indices = os.listdir(join(datasets_dir, "GNN_input_data"))
test_indices = [index[:index.rfind("_")] for index in test_indices]
test_indices = list(set([index for index in test_indices if "test" in index]))

param_grid = {'batch_size': [32,64,96],
                'D': [50,100],
                'learning_rate': [0.01, 0.1],
                'epochs': [30,50,80],
                'l2_reg_fc' : [0.02, 0.05, 0.08],
                'l2_reg_conv': [0.02, 0.05, 0.08],
                'rho': [0.9, 0.95, 0.99]}

params_list = [(batch_size, D, learning_rate, epochs, l2_reg_fc, l2_reg_conv, rho) for batch_size in param_grid['batch_size'] for D in param_grid["D"] for learning_rate in param_grid['learning_rate']
                for epochs in param_grid['epochs'] for l2_reg_fc in param_grid['l2_reg_conv'] for l2_reg_conv in param_grid['l2_reg_conv'] for rho in param_grid["rho"]]

def hyperparameter_tuning(params):

    batch_size, D, learning_rate, epochs, l2_reg_fc, l2_reg_conv, rho = params
    MSE = []

    for i in range(5):
        train_index, test_index  = CV_train_indices[i], CV_test_indices[i]
        train_index = [ind for ind in train_indices if int(ind.split("_")[1]) in train_index]
        test_index = [ind for ind in train_indices if int(ind.split("_")[1]) in test_index]

        train_params = {'batch_size': batch_size,
                'folder' :join(datasets_dir, "GNN_input_data"),
                'list_IDs' : np.array(train_index),
                'shuffle': True}

        test_params = {'batch_size': len(test_index),
                'folder' : join(datasets_dir, "GNN_input_data"),
                'list_IDs' : np.array(test_index),
                'shuffle': False}

        training_generator = DataGenerator(**train_params)
        test_generator = DataGenerator(**test_params)


        model = DMPNN_without_extra_features(l2_reg_conv = l2_reg_conv, l2_reg_fc = l2_reg_fc, learning_rate = learning_rate,
                    D = D, N = N, F1 = F1, F2 = F2, F= F, drop_rate = 0.0, ada_rho = rho)
        model.fit(training_generator, epochs= epochs, shuffle = True, verbose = 1)

        #get test_y:
        test_indices_y = [int(ind.split("_")[1]) for ind in train_indices if ind in test_index]
        test_y = np.array([train_df["kcat"][ind] for ind in test_indices_y])

        pred_test = model.predict(test_generator)
        print(np.mean(abs(pred_test - np.reshape(test_y[:len(pred_test)], (-1,1)))**2))
        MSE.append(np.mean(abs(pred_test - np.reshape(test_y[:len(pred_test)], (-1,1)))**2))

    result = {"batch_size" : batch_size, "D" : D , "learning_rate" : learning_rate, "epochs" : epochs,
                    "l2_reg_fc" : l2_reg_fc, "l2_reg_conv" : l2_reg_conv, "rho" : rho, "cv_mse" : np.mean(MSE)}
    
    return(result)

if __name__ == '__main__':

    results=[]
    params_list = random.sample(params_list, 100)

    with multiprocessing.Pool(multiprocessing.cpu_count()-1) as pool:
        for result in pool.imap(hyperparameter_tuning, params_list):
            with open('tuning.txt','a') as data:  
                data.write(str(result) + '\n')
            results.append(result)

    params = min(results, key=lambda d: d['cv_mse'])

    print(params)
