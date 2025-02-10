from subprocess import call
import os
import argparse
import re
from pathlib import Path
# 1-33

# Example usage:
# python binderDesign.py -target dl_binder_design/Rosetta/rosettaInput/11_14_24_Novel_Target_1_Testy/output/spliced_seq50_from_110124Target_0001Testy.pdb -hotspot "[A17, A18, A19]" -bindersize "50-150" -chainrange "A1-33" -designcount 25 -mpnnsequences 8

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

def silentToFasta(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        count = 1
        for line in infile:
            if line.startswith('ANNOTATED_SEQUENCE:'):
                # Extract the sequence part
                sequence = line.split(':', 1)[1].strip()
                
                # Stop parsing after the first space character
                sequence = sequence.split(' ')[0]
                
                # Find the position of the first '[' and split one character before it
                split_pos = sequence.find('[') - 1
                if split_pos > 0:
                    left_side = sequence[:split_pos + 1]  # Include the character before '['
                    right_side = sequence[split_pos + 1:]
                    
                    # Remove all content within brackets from the right side
                    right_side = re.sub(r'\[.*?\]', '', right_side)
                    
                    processed_sequence = f"{left_side}:{sequence[split_pos]}{right_side}"
                else:
                    processed_sequence = sequence
                
                # Write to the output file
                outfile.write(f">sequence{count}\n")
                outfile.write(f"{processed_sequence}\n")
                count += 1

def runRFDiffusion(rfDiffusionConda, rfDiffusionPath, rfDiffusionOutputPath, targetPath, pdbName, chainRange, binderSize, designCount, hotSpots):
    rfDiffusionCmd = [
        rfDiffusionConda,
        f'{rfDiffusionPath}/scripts/run_inference.py',
        f'inference.output_prefix={rfDiffusionOutputPath}/{pdbName}/{pdbName}',
        f'inference.input_pdb={targetPath}',
        f'contigmap.contigs=[{chainRange}/0 {binderSize}]',
        f'ppi.hotspot_res={hotSpots}',
        f'inference.num_designs={designCount}',
        f'contigmap.inpaint_str=[{chainRange}]',
        'denoiser.noise_scale_ca=0',
        'denoiser.noise_scale_frame=0']    
    print(rfDiffusionCmd)
    call(rfDiffusionCmd)

    outputFiles = f'{rfDiffusionOutputPath}/{pdbName}/*.pdb'

    return outputFiles

def runProteinMPNN(proteinMPNNConda, proteinMPNNPath, proteinMPNNOutputPath, silentToolsPath, diffusedBinderPaths, pdbName, sequenceCount=1):
    if sequenceCount != 1:
        relaxationCount = 0
    else:
        relaxationCount = 1
    
    pdbToSilentConversionCmd = [
        '/bin/bash',
        '-c',  # Use -c to pass the command as a string
        f'{silentToolsPath}/silentfrompdbs {diffusedBinderPaths} > {pdbName}Binders.silent'
    ]

    outputSilentPath = f'{proteinMPNNOutputPath}/{pdbName}/{pdbName}SequencedBinders.silent'

    proteinMPNNCmd = [
        proteinMPNNConda,
        f'{proteinMPNNPath}/dl_interface_design.py',
        '-silent',
        f'{pdbName}Binders.silent',
        '-outsilent',
        outputSilentPath,
        '-relax_cycles',
        f'{relaxationCount}',
        '-seqs_per_struct',
        f'{sequenceCount}',
        '-debug'
    ]
    print(proteinMPNNCmd)
    Path(f'{proteinMPNNOutputPath}/{pdbName}/').mkdir(parents=True, exist_ok=True)

    print('Converting pdb files to silent...')
    call(pdbToSilentConversionCmd)
    print("Running protein MPNN...")
    call(proteinMPNNCmd)
    
    os.rename('check.point', f'{proteinMPNNOutputPath}/{pdbName}/check.point')
    os.rename(f'{pdbName}Binders.silent', f'{proteinMPNNOutputPath}/{pdbName}/{pdbName}Binders.silent')
    os.rename(f'{pdbName}Binders.silent.idx', f'{proteinMPNNOutputPath}/{pdbName}/{pdbName}Binders.silent.idx')


    return outputSilentPath

def runAlphaFoldInitialGuess(alphaFoldConda, alphaFoldPath, alphaFoldOutputPath, pdbName, binderSequenceSilentPath):
    
    Path(f'{alphaFoldOutputPath}/{pdbName}BinderOutput/').mkdir(parents=True, exist_ok=True)
    
    alphaFoldCmd = [
        alphaFoldConda,
        f'{alphaFoldPath}/predict.py',
        '-silent',
        f'{binderSequenceSilentPath}',
        '-outsilent',
        f'{alphaFoldOutputPath}/{pdbName}BinderOutput/{pdbName}FoldOutput.silent'
    ]
    call(alphaFoldCmd)

    outputFile = f'{alphaFoldOutputPath}/{pdbName}BinderOutput/{pdbName}Out.sc'
    os.rename(f'check.point', f'{alphaFoldOutputPath}/{pdbName}BinderOutput/{pdbName}Check.point')
    os.rename(f'out.sc', outputFile)

    return outputFile

args, configArgs = obtainArguments()

print(configArgs)

pdbName = os.path.splitext(os.path.basename(args.target))[0]

print("Running RFDiffusion...")
diffusedOutputFiles = runRFDiffusion(configArgs['rfDiffusionConda'], configArgs['rfDiffusionBasePath'], configArgs['rfDiffusionOutputPath'], args.target, pdbName, args.chainrange, args.bindersize, args.designcount, args.hotspot)

print("Running Protein MPNN...")
silentOutputFile = runProteinMPNN(configArgs['proteinMPNNConda'], configArgs['proteinMPNNBasePath'], configArgs['proteinMPNNOutputPath'], configArgs['silentToolsBasePath'], diffusedOutputFiles, pdbName, args.mpnnsequences)

print("Running AlphaFold...")
alphaFoldOutputFile = runAlphaFoldInitialGuess(configArgs['alphaFoldConda'], configArgs['alphaFoldBasePath'], configArgs['alphaFoldOutputPath'], pdbName, silentOutputFile)

# Steps for use (with simplest file organization):
# Install conda
# Clone https://github.com/nrbennet/dl_binder_design/
# If you want to refactor proteins, install Rosetta from https://www.rosettacommons.org/software/ and add to path
    # Create a folder "Rosetta" inside the dl_binder_design directory. Inside, create two folders rosettaApp and rosettaInput. Install rosetta inside the rosettaApp folder
    # Create a conda environment for this application
# Install RF Diffusion via this source: https://github.com/RosettaCommons/RFdiffusion
    # Create a folder "RFdiffusion" inside the dl_binder_design directory. Inside, create a folder RFDiffusionOutput. Install RFDiffusion inside the RFDiffusion folder
    # Create a conda environment for this application
# Install ProteinMPNN and the AlphaFold initial prediction via the README.md instructions in the dl_binder_design folder
# Create a binderDesignConfig.txt file using the following example:
    # alphaFoldConda=/ssd1/home/safwand/miniconda3/envs/af2_binder_design/bin/python
    # alphaFoldBasePath=dl_binder_design/af2_initial_guess
    # alphaFoldOutputPath=dl_binder_design/af2_initial_guess/output

    # proteinMPNNConda=/ssd1/home/safwand/miniconda3/envs/proteinmpnn_binder_design/bin/python
    # silentToolsBasePath=dl_binder_design/include/silent_tools
    # proteinMPNNBasePath=dl_binder_design/mpnn_fr
    # proteinMPNNOutputPath=dl_binder_design/mpnn_fr/output

    # rfDiffusionConda=/ssd1/home/safwand/miniconda3/envs/SE3nv/bin/python
    # rfDiffusionBasePath=dl_binder_design/RFDiffusion/RFdiffusion
    # rfDiffusionOutputPath=dl_binder_design/RFDiffusion/RFDiffusionOutput