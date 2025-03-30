from subprocess import call
import os
import argparse
import re
from pathlib import Path
import shutil

def combine_and_renumber_blueprints(blueprint_file1, blueprint_file2, output_file="combined_blueprint.txt"):
    """
    Combines two Blueprint files, renumbers the residues, and writes to a new Blueprint file.
    """
    # Function to parse a Blueprint file and extract residues info
    def parse_blueprint(blueprint_file):
        residues = []
        with open(blueprint_file, 'r') as f:
            for line in f:
                line = line.strip()  # Remove leading/trailing whitespace
                if line:  # Process non-empty lines
                    parts = line.split()
                    residue_seq = int(parts[0])  # Original residue sequence number
                    residue_name = parts[1]  # Residue name (e.g., A, L, G)
                    secondary_structure = parts[2] if len(parts) > 2 else '.'  # Secondary structure (default is '.')
                    residues.append({
                        'residue_seq': residue_seq,
                        'residue_name': residue_name,
                        'secondary_structure': secondary_structure
                    })
        return residues

    # Parse both Blueprint files
    residues1 = parse_blueprint(blueprint_file1)
    residues2 = parse_blueprint(blueprint_file2)

    # Combine residues from both files
    combined_residues = residues1 + residues2

    # Renumber the residues and write the combined Blueprint
    next_residue_number = 1  # Start renumbering from 1

    with open(output_file, 'w') as f:
        for residue in combined_residues:
            # Write renumbered residue info to the output file
            f.write(f"{next_residue_number} {residue['residue_name']} {residue['secondary_structure']}\n")
            next_residue_number += 1

    print(f"Combined and renumbered Blueprint file written to {output_file}")

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
    parser.add_argument("-target", help="Target PDB file path", required=True)
    parser.add_argument("-hotspot", help="Hotspots on the target to focus binding e.g: \'[A23, A24, A25]\'", required=True)
    parser.add_argument("-bindersize", help="Binder size range (e.g. 150-200)", required=True)
    parser.add_argument("-chainrange", help="Target chain range (e.g. A1-86)", required=True)
    parser.add_argument("-designcount", help="Number of RFDiffusion designs wanted", required=True)
    parser.add_argument("-mpnnsequences", help="Number of ProteinMPNN sequences wanted")

    configArgs = read_config('binderDesignConfig.txt')

    args = parser.parse_args()
    return args, configArgs

def find_and_replace_motif(filename, find_seq, replace_seq):
    with open(filename, 'r') as file:
        lines = file.readlines()
        
    find_len = len(find_seq)
    replace_len = len(replace_seq)
    
    if find_len != replace_len:
        print("Error: The find and replace sequences must be of the same length.")
        return
    
    sequence_indices = []
    
    for i in range(len(lines) - find_len + 1):
        match = True
        for j in range(find_len):
            current_line = lines[i + j].strip().split()
            if len(current_line) < 2 or current_line[1] != find_seq[j]:
                match = False
                break
        if match:
            sequence_indices.append(i)
    
    if len(sequence_indices) == 1:
        index = sequence_indices[0]
        for j in range(find_len):
            current_line = lines[index + j].strip().split()
            if current_line[1] != replace_seq[j]:
                current_line.append('PIKAA')
                current_line.append(replace_seq[j])
            lines[index + j] = ' '.join(current_line) + '\n'
        with open(filename, 'w') as file:
            file.writelines(lines)
        print(f"Formatted sequence '{find_seq}' with '{replace_seq}' starting at row {index + 1}")
    elif len(sequence_indices) > 1:
        print("Error: More than one sequence found. No changes made.")
    else:
        print("No such sequence found.")

def filterGoodBinders(outputFile, silentOutputFile, silentToolsPath, alphaFoldOutputPath, rosettaBasePath, pdbName):
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
    tempFile = f'{alphaFoldOutputPath}/{pdbName}BinderOutput/temp'
    Path(tempFile).mkdir(parents=True, exist_ok=True)
    Path(f'{tempFile}/testy').mkdir(parents=True, exist_ok=True)

    # Obtain the good binders via metric filtering
    for binder in dict_list:
        pae_interaction = float(binder['pae_interaction'])
        plddt_binder = float(binder['plddt_binder'])
        plddt_target = float(binder['plddt_target'])
        if pae_interaction < 5 and plddt_binder >= 90 and plddt_target >= 90:
            description = binder['description']
            description = description.replace('_af2pred', '')
            silentExtractCmd = [
                '/bin/bash',
                '-c',  # Use -c to pass the command as a string
                f'{silentToolsPath}/silentextractspecific {silentOutputFile} {description}'
            ]
            print(f'Extracting {description}...')
            call(silentExtractCmd)
            pdbFile = f'{tempFile}/{description}.pdb'
            os.rename(f'{description}.pdb', pdbFile)

            blueprintCreateCmdChainA = [
                '/bin/bash',
                '-c',  # Use -c to pass the command as a string
                f'{rosettaBasePath}/tools/remodel/getBluePrintFromCoords.pl -pdbfile {pdbFile} -chain A > {description}BlueprintChainA.txt'
            ]

            blueprintCreateCmdChainB = [
                '/bin/bash',
                '-c',  # Use -c to pass the command as a string
                f'{rosettaBasePath}/tools/remodel/getBluePrintFromCoords.pl -pdbfile {pdbFile} -chain B > {description}BlueprintChainB.txt'
            ]
            call(blueprintCreateCmdChainA)
            call(blueprintCreateCmdChainB)
            combine_and_renumber_blueprints(f'{description}BlueprintChainA.txt', f'{description}BlueprintChainB.txt', f'{description}Blueprint.txt')
            os.rename(f'{description}BlueprintChainA.txt', f'{tempFile}/{description}BlueprintChainA.txt')
            os.rename(f'{description}BlueprintChainB.txt', f'{tempFile}/{description}BlueprintChainB.txt')
            os.rename(f'{description}Blueprint.txt', f'{tempFile}/{description}Blueprint.txt')

            find_and_replace_motif(f'{tempFile}/{description}Blueprint.txt', 'DGR', 'NGR')

            rosettaRemodelCmd = [
                '/bin/bash',
                '-c',  # Use -c to pass the command as a string
                f'{rosettaBasePath}/source/bin/remodel.linuxgccrelease -s {pdbFile} -remodel:blueprint {tempFile}/{description}Blueprint.txt -chain B -no_optH false -ex1 -ex2 -use_input_sc -linmem_ig 10 -num_trajectory 5 -save_top 0'
            ]
            call(rosettaRemodelCmd)
            os.rename('score.sc', f'{tempFile}/score.sc')
            os.rename(f'{description}_0001.pdb', f'{tempFile}/testy/{description}.pdb')

    goodBinderOutputFile = f'{alphaFoldOutputPath}/{pdbName}BinderOutput/{pdbName}GoodBinders.silent'
    silentCreateCmd = [
        '/bin/bash',
        '-c',  # Use -c to pass the command as a string
        f'{silentToolsPath}/silentfrompdbs {tempFile}/testy/*.pdb > {goodBinderOutputFile}'
    ]
    call(silentCreateCmd)

    shutil.rmtree(tempFile)
    return goodBinderOutputFile

def runComparison(alphaFoldConda, alphaFoldPath, alphaFoldOutputPath, pdbName, binderSequenceSilentPath):
    Path(f'{alphaFoldOutputPath}/{pdbName}BinderOutput/NGRComparison').mkdir(parents=True, exist_ok=True)
    
    alphaFoldCmd = [
        alphaFoldConda,
        f'{alphaFoldPath}/predict.py',
        '-silent',
        f'{binderSequenceSilentPath}',
        '-outsilent',
        f'{alphaFoldOutputPath}/{pdbName}BinderOutput/NGRComparison/{pdbName}FoldOutputComparison.silent'
    ]
    call(alphaFoldCmd)

    outputFile = f'{alphaFoldOutputPath}/{pdbName}BinderOutput/NGRComparison/{pdbName}OutComparison.sc'
    os.rename(f'check.point', f'{alphaFoldOutputPath}/{pdbName}BinderOutput/NGRComparison/{pdbName}CheckComparison.point')
    os.rename(f'out.sc', outputFile)

    return

args, configArgs = obtainArguments()

print(configArgs)

pdbName = os.path.splitext(os.path.basename(args.target))[0]

print("Filtering good binders...")
goodBinderOutputFile = filterGoodBinders('dl_binder_design/af2_initial_guess/output/short_dual_helix_target_relaxed_rank_001_alphafold2_ptm_model_3_seed_000BinderOutput/short_dual_helix_target_relaxed_rank_001_alphafold2_ptm_model_3_seed_000Out.sc', 'dl_binder_design/mpnn_fr/output/short_dual_helix_target_relaxed_rank_001_alphafold2_ptm_model_3_seed_000/short_dual_helix_target_relaxed_rank_001_alphafold2_ptm_model_3_seed_000SequencedBinders.silent', configArgs['silentToolsBasePath'], configArgs['alphaFoldOutputPath'], configArgs['rosettaBasePath'], pdbName)

print("Running AlphaFold comparison...")
runComparison(configArgs['alphaFoldConda'], configArgs['alphaFoldBasePath'], configArgs['alphaFoldOutputPath'], pdbName, goodBinderOutputFile)
