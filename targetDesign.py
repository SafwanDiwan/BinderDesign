from subprocess import call, Popen
import os
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser
import argparse 
import json
import glob
import numpy as np
import shutil
from pathlib import Path

def read_config(file_path):
    config = {}
    with open(file_path, 'r') as file:
        for line in file:
            if '=' in line:
                key, value = line.strip().split('=', 1)
                config[key] = value
    return config

def obtainArguments():
    parser = argparse.ArgumentParser(description="Given a target PDB file, this script will run the necessary software to design and score binders for that target.")
    parser.add_argument("-inputjson", help="Input JSON config file", required=True)
    parser.add_argument("-batchsize", help="ColabFold batch size for concurrent folding", required=True)

    configArgs = read_config('binderDesignConfig.txt')

    args = parser.parse_args()
    return args, configArgs

def pdb_to_fasta(pdb_file, output_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    chain = next(structure.get_chains())  # Get the first (and assumed only) chain
    sequence = "".join([seq1(residue.resname) for residue in chain if residue.id[0] == " "])

    output_name = os.path.splitext(os.path.basename(output_file))[0]
    fasta_content = f'>{output_name}\n{sequence}\n'

    with open(f"{output_file}.fasta", "w") as file:
        file.write(fasta_content)

    return f"{output_file}.fasta"


def runProteinGenerator(proteinGeneratorConda, proteinGeneratorBasePath, inputJson, colabFoldConda, colabFoldBasePath, inputFasta, outputFolder):
    proteinGeneratorCmd = [
        proteinGeneratorConda,
        f'{proteinGeneratorBasePath}/inference.py',
        '--input_json',
        inputJson,
    ]
    call(proteinGeneratorCmd) # Run protein generator - wait for finish so that we can evaluate via AF2

def createFastaFiles(outFolder):
    # Get FASTA from the output of protein generator
    fastaPathList = []
    pdbFiles = glob.glob(f'{outFolder}/*.pdb')
    if not pdbFiles:
        print("No PDB files found in the output folder.")
        return
    for file in pdbFiles:
        fileName = os.path.splitext(os.path.basename(file))[0]
        fastaPath = f'{os.path.splitext(file)[0]}/'
        os.makedirs(fastaPath, exist_ok=True)
        fastaPath = f'{fastaPath}{fileName}'

        pathWithExtension = pdb_to_fasta(file, fastaPath) # Convert PDB to FASTA
        fastaPathList.append(pathWithExtension)
    
    return fastaPathList

def runColabFold(fastaPaths, colabFoldConda, colabFoldBasePath, batchSize):
    foldDirectories = []
    processes = []  # To keep track of the started processes

    for i, fastaPath in enumerate(fastaPaths):
        # Prepare the command for colabFold
        colabFoldCmd = [
            colabFoldConda,
            colabFoldBasePath,
            '--templates',
            '--amber',
            fastaPath,
            os.path.dirname(fastaPath)
        ]
        call(colabFoldCmd)
        foldDirectories.append(os.path.dirname(fastaPath))
        
        # Start the process
      #  process = Popen(colabFoldCmd)
      #  processes.append(process)
        
        # Wait for processes to complete when batch size is reached
      #  if (i + 1) % batchSize == 0 or (i + 1) == len(fastaPaths):
            # Wait for all processes in the batch to finish
       #     for p in processes:
        #        p.wait()
            # Collect results after batch is done
         #   foldDirectories.extend([os.path.dirname(fastaPaths[j]) for j in range(i + 1 - batchSize, i + 1)])

            # Clear the process list for the next batch
          #  processes = []
   # print(foldDirectories)
    return foldDirectories

def filterGoodTargets(foldDirectories, goodTargetDirectory):
    goodDirectories = []
    for directory in foldDirectories:
        print(f'Searching through directory {directory}')
        for filename in os.listdir(directory):
            # Check if the pattern exists in the file name
            if 'scores_rank_001' in filename:
                scorePath = os.path.join(directory, filename)

                with open(scorePath, "r") as f:
                    data = json.load(f)
        
                # Extract the pLDDT matrix
                pLDDT = np.array(data["plddt"])
                avgpLDDT = np.mean(pLDDT)
                if avgpLDDT >= 90:
                    print(f'Found a good target in {scorePath}')
                    goodDirectories.append(f'{directory}/{os.path.basename(os.path.normpath(directory))}.fasta')
                    for filename in os.listdir(directory):
                        # Move to folder with dates that has good binder targets in it
                        if '_relaxed_rank_001' in filename:
                            pdbPath = os.path.join(directory, filename)
                            os.makedirs(goodTargetDirectory, exist_ok=True)
                            shutil.copy(pdbPath, goodTargetDirectory)
                            break
                    break
    return goodDirectories

args, configArgs = obtainArguments()

with open(args.inputjson, "r") as file:
    jsonData = json.load(file)
    
outputFolder = jsonData['out']

if not os.path.exists(outputFolder):
    os.makedirs(outputFolder)


runProteinGenerator(configArgs['proteinGeneratorConda'], configArgs['proteinGeneratorBasePath'], args.inputjson, configArgs['colabFoldConda'], configArgs['colabFoldBasePath'], 'inputFasta', 'outputFolder')

fastaPaths = createFastaFiles(outputFolder)

foldDirectories = runColabFold(fastaPaths, configArgs['colabFoldConda'], configArgs['colabFoldBasePath'], int(args.batchsize))

goodDirectories = filterGoodTargets(foldDirectories, f'{root_path}/ViableTargets')

# python targetDesign.py -inputjson designSoftware/protein_generator/examples/out/design_000000_args.json
