from subprocess import call, Popen
import os
import argparse
import numpy as np
from pathlib import Path
import shutil
from Bio import SeqIO
from Bio.SeqUtils import seq1
from Bio.Seq import Seq
from Bio.PDB import PDBParser
import glob
import json

def read_config(file_path):
    config = {}
    with open(file_path, 'r') as file:
        for line in file:
            if '=' in line:
                key, value = line.strip().split('=', 1)
                config[key] = value
    return config

def obtainArguments():
    parser = argparse.ArgumentParser(description="Given a target PDB file, this script will run the necessary software compare selectivity of already created and scored binder files.")
    parser.add_argument("-outSCFile", help="out.sc file containing scores of binder design process", required=True)
    parser.add_argument("-dgrBinderSilentFile", help="Silent file containing the binder/target pairs from the binder design process", required=True)
    parser.add_argument("-originalMotif", help="The original motif to be replaced", required=True)
    parser.add_argument("-replacementMotif", help="The motif to replace the original motif", required=True)
    parser.add_argument("-target", help="The target PDB file", required=True)

    configArgs = read_config('binderDesignConfig.txt')

    args = parser.parse_args()
    return args, configArgs

def pdb_to_fasta(pdb_file, outputFile, replacementMotif):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    chains_sequences = {}

    for model in structure:
        for chain in model:
            chain_id = chain.id
            sequence = "".join([seq1(residue.resname) for residue in chain if residue.id[0] == " "])
            chains_sequences[chain_id] = sequence

    binder_len = len(chains_sequences["A"])
    outputName = os.path.splitext(os.path.basename(outputFile))[0]
    with open(f'{outputFile}.fasta', 'w') as file:
        file.write(f'>{outputName}_{replacementMotif}\n{chains_sequences["A"]}:{chains_sequences["B"]}')

    return binder_len

def find_and_replace_motif(fasta_path, target_seq, replacement_seq):
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    modified = False
    
    sequences = list(SeqIO.parse(fasta_path, "fasta"))
    modified = False
    
    for record in sequences:
        seq_str = str(record.seq).replace("\n", "").replace("\r", "").strip()  # Remove all newlines and whitespace
        seq_parts = seq_str.split(":", 1)
        if len(seq_parts) == 2 and target_seq in seq_parts[1]:
            seq_parts[1] = seq_parts[1].replace(target_seq, replacement_seq)
            record.seq = Seq(seq_parts[0] + ":" + seq_parts[1])
            modified = True
    
    if modified:
        with open(fasta_path, "w") as output_file:
            for record in sequences:
                output_file.write(f">{record.id}\n{record.seq}\n")  # Write sequence in single line
    
    return modified  # Returns True if any replacements were made

'''
Filters good binders when given an out.sc file from the binder design process.
Parses the silent file according to the out.sc file and extracts the good binders.

@param outputFile: The out.sc file from the binder design process
@param silentOutputFile: The silent file from the binder design process
@param silentToolsPath: The path to the silent tools
@param colabFoldOutputPath: The path to the output folder
@param pdbName: The name of the PDB file
@return: The path to the output folder containing the good binders
'''

def filterGoodBinders(outputFile, silentOutputFile, silentToolsPath, colabFoldOutputPath, pdbName, originalMotif, replacementMotif):
    # Initialize a list to hold the dictionaries
    dict_list = []
    # Parse data into dictionaries
    with open(outputFile, 'r') as file:
        
        lines = file.readlines()
        
        # Obtain the keys from the first line
        keys = lines[0].split()
        keys = keys[1:]
        
        # Process the remaining lines
        for line in lines[1:]:
            values = line.split()
            values = values[1:]
            row_dict = {keys[i]: values[i] for i in range(len(keys))}
            dict_list.append(row_dict)

    # Create a temporary file that contains the pdb files of the binder/target pairs from AF2_Initial_Guess
    tempFile = f'{colabFoldOutputPath}/{pdbName}BinderOutput/temp'
    Path(tempFile).mkdir(parents=True, exist_ok=True)
    fastaDirectoriesBinderLengths = {}
    # Obtain the good binders via metric filtering
    for binder in dict_list:
        pae_interaction = float(binder['pae_interaction'])
        plddt_binder = float(binder['plddt_binder'])
        plddt_target = float(binder['plddt_target'])
        if pae_interaction < 5 and plddt_binder >= 90 and plddt_target >= 90:
            description = binder['description']
            description = description.replace('_af2pred', '')

            # Extract the pdb file from the silent file
            silentExtractCmd = [
                '/bin/bash',
                '-c',  # Use -c to pass the command as a string
                f'{silentToolsPath}/silentextractspecific {silentOutputFile} {description}'
            ]
            print(f'Extracting {description}...')
            call(silentExtractCmd)

            # Move the pdb file to the temp folder
            pdbFile = f'{tempFile}/{description}.pdb'
            os.rename(f'{description}.pdb', pdbFile)

            # Convert the pdb file to fasta and store the fasta file in the colabFoldOutput folder passed in
            fastaFilePath = f'{colabFoldOutputPath}/{pdbName}BinderOutput/{description}'
            binderLen = pdb_to_fasta(pdbFile, fastaFilePath, replacementMotif)
            fastaDirectoriesBinderLengths[f'{fastaFilePath}.fasta'] = binderLen

            # Replace the DGR motif with NGR
            find_and_replace_motif(f'{fastaFilePath}.fasta', originalMotif, replacementMotif)

    # Remove the temporary folder
    shutil.rmtree(tempFile)

    return fastaDirectoriesBinderLengths

def runComparison(colabFoldConda, colabFoldBasePath, fastaDirectories):
    processes = []
    for inputFasta in fastaDirectories:
        outputFolder = inputFasta.replace('.fasta', '')
        colabFoldCmd = [
            colabFoldConda,
            colabFoldBasePath,
            '--templates',
            '--amber',
            inputFasta,
            outputFolder
        ]
        process = Popen(colabFoldCmd)  # Start the process
        processes.append(process)  # Store the process in the list

    # Wait for all processes to complete
    for process in processes:
        process.wait()  # Ensures all processes finish before moving on


def scoreBinders(fastaDirectories):
    for directory in fastaDirectories.keys():
        outputFolder = directory.replace('.fasta', '')
        searchPattern = os.path.join(outputFolder, f"*scores_rank_001*.json")
        filePath = glob.glob(searchPattern)[0]
        print(filePath)
        
        # Load the JSON file containing PAE data (adjust the filename)
        with open(filePath, "r") as f:
            data = json.load(f)

        # Extract the PAE matrix
        pae = np.array(data["pae"])  # Convert list of lists to NumPy array

        # Get binder length
        binderLength = fastaDirectories.get(directory)

        # Compute PAE metrics
        pae_interaction1 = np.mean(pae[:binderLength, binderLength:])
        pae_interaction2 = np.mean(pae[binderLength:, :binderLength])
        pae_binder = np.mean(pae[:binderLength, :binderLength])
        pae_target = np.mean(pae[binderLength:, binderLength:])
        pae_interaction_total = (pae_interaction1 + pae_interaction2) / 2

        # Print results
        with open(f"{os.path.dirname(outputFolder)}/out.sc", "w") as sc_file:
            sc_file.write("Directory\tPAE Interaction 1\tPAE Interaction 2\tPAE Binder\tPAE Target\tPAE Interaction Total\n")
            sc_file.write(f"{directory}\t{pae_interaction1:.4f}\t{pae_interaction2:.4f}\t{pae_binder:.4f}\t{pae_target:.4f}\t{pae_interaction_total:.4f}\n")

    return

args, configArgs = obtainArguments()

pdbName = os.path.splitext(os.path.basename(args.target))[0]

print("Filtering good binders...")
fastaDirectories = filterGoodBinders(args.outSCFile, args.dgrBinderSilentFile, configArgs['silentToolsBasePath'], configArgs['colabFoldOutputPath'], pdbName, args.originalMotif, args.replacementMotif)

print("Running ColabFold comparison...")
runComparison(configArgs['colabFoldConda'], configArgs['colabFoldBasePath'], fastaDirectories.keys())

print("Scoring Binders...")
scoreBinders(fastaDirectories)

# python selectivityEvaluation2.py -outSCFile dl_binder_design/af2_initial_guess/output/short_dual_helix_target_relaxed_rank_001_alphafold2_ptm_model_3_seed_000BinderOutput/short_dual_helix_target_relaxed_rank_001_alphafold2_ptm_model_3_seed_000Out.sc -dgrBinderSilentFile dl_binder_design/mpnn_fr/output/short_dual_helix_target_relaxed_rank_001_alphafold2_ptm_model_3_seed_000/short_dual_helix_target_relaxed_rank_001_alphafold2_ptm_model_3_seed_000SequencedBinders.silent -originalMotif DGR -replacementMotif NGR -target colabfold/experiments/short_dual_helix_target/short_dual_helix_target_relaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb