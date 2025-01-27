To run this codebase, you need the following installed 

Setup the base for the env (Python 3.7) including rdkit using:
 
conda create -c conda-forge -n ee python=3.7 rdkit

Specific packages and version:

pip install ete3==3.1.3 zeep==4.2.1 pandas==1.1.3 bioservices==1.11.2 colorlog==6.8.0 easydev==0.12.1 ete3==3.1.3 keras==2.10.0 numpy==1.21.5 pandas==1.1.3 tensorflow==2.10.0 torch==1.13.1 fair-esm==2.0.0

For some reason there is also this after fair-esm , 0.4.0

to run the baseline for kcat: 

pip install hyperopt==0.2.7 matplotlib==3.5.3 plotnine==0.8.0 scipy==1.7.3 xgboost==1.5.0