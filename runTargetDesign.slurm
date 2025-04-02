#!/bin/bash
#SBATCH --job-name=targetDesign
#SBATCH --output=logs/%j_output.log
#SBATCH --gres=gpu:4          # Request 4 GPUs
#SBATCH --cpus-per-task=8      # Adjust based on workload
#SBATCH --mem=72G              # Adjust memory
#SBATCH --time=96:00:00        # Set job duration
#SBATCH --partition=gpu-rtx6k        # Adjust if needed for your cluster
#SBATCH --account=ml4ml

# Activate virtual environment if needed
source ~/miniconda3/bin/activate proteingenerator

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/mmfs1/gscratch/ml4ml/safwand/designSoftware/colabfold/localcolabfold/conda/lib"

# Run the Bash script that executes Python scripts sequentially
python targetDesign.py -inputjson PipelineOutput/TargetDesign/exp_033025_1/exp_033025_1_args.json -batchsize 3
