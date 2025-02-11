Hi there!

This is the codebase for the Enkyme pipeline developed to predict Km and Kcat for enzymes involved in glucosinolate metabolism.
This pipeline is based on the work of Kroll et.al: https://github.com/AlexanderKroll/KM_prediction_function/tree/main and https://github.com/AlexanderKroll/kcat_prediction_function/tree/main

To run this codebase, you need the following installed:
-Python 3.7
-ete3==3.1.3
-zeep==4.2.1
-pandas==1.1.3
-bioservices==1.11.2
-colorlog==6.8.0
-easydev==0.12.1
-ete3==3.1.3
-keras==2.10.0
-numpy==1.21.5
-pandas==1.1.3
-tensorflow==2.10.0
-torch==1.13.1
-rdkit==2020.09.1.0
-fair-esm==2.0.0, 0.4.0 to run the baseline for kcat
-hyperopt==0.2.7
-matplotlib==3.5.3
-plotnine==0.8.0
-scipy==1.7.3
-xgboost==1.5.0

Because part of the hyperparameter tuning was run on online computing services with a different version of XGBoost, it might be possible that the results deviate from what is presented in the paper.

To run the pipeline, do the following:

1. Navigate to the preprocess folder in either the kcat or Km folder under the code folder
2. > python sabio_download_for_model.py
   > python uniprot sequence.py
3. Run the 01 - Data preprocessing notebook
4. Navigate to the modelling folder
5. Run the 01 Training xgboost models notebook
6. Run the 02 Analyzing results notebook

If you want to run the baselines, do the following
1. Navigate to the baseline folder in the code folder
2. Run the Tutorial kcat prediction or KM prediction notebook
