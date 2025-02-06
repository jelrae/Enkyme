#!/bin/bash
#SBATCH --job-name=Kcat_pipeline
#SBATCH --output=./server_output_files/kcat_pipe-%a.out
#SBATCH --error=./server_output_files/kcat_pipe-%a.err
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks-per-node=1           # Number of tasks per node
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=j.a.earle@uva.nl    # Email address for notifications
#SBATCH --mem-per-cpu=8000

eval "$(conda shell.bash hook)"
conda activate ee
cd /home/jearle/personal/Enkyme/Kcat/code/preprocess/
python sabio_download_for_model.py python uniprot sequence.py