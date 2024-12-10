# Binder Design

This repository leverages a wide assortment of computational tools to provide an easy to use pipeline to develop binder targets (via Rosetta Remodel) and binders (via RFDiffusion, ProteinMPNN FastRelax, and AlphaFold 2).

Setup can be tricky, so be sure to complete every step of the following instructions:
- Clone this repository. 
- Make sure you have conda installed.
- Inside it, clone the dl_binder_design repository found at https://github.com/nrbennet/dl_binder_design
- Within this repository create two folders, 'RFDiffusion' and 'Rosetta'. Inside the 'RFDiffusion' folder create a 'RFDiffusionOutput' folder. Inside the 'Rosetta' folder create 'rosettaApp' and 'rosettaInput' folders.
- Inside the 'RFDiffusion' folder, clone the RFDiffusion repository found here: https://github.com/RosettaCommons/RFdiffusion
- Inside the 'rosettaApp' folder, download and install Rosetta from here (make sure to follow the directions on how to install the suite): https://docs.rosettacommons.org/demos/latest/tutorials/install_build/install_build
- You will need four conda environments one for rosetta (should automatically be created when you install rosetta), another for RFDiffusion (the instructions to create are included in the RFDiffusion GitHub README), a third for proteinmpnn_fast_relax, and a fourth for af2_initial_guess. .yml files for each environment can be found in the environmentConfigs folder.
- There are also a couple of changes to three of the silent_tools files that you'll have to make in order to ensure that they run with the correct conda environment. 
    - In the silentls and silentrename files, change 'sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))' to 'sys.path.append('dl_binder_design/include/silent_tools/silent_tools.py')'
    - In the silentfrompdbs file there are two changes that must be made:
        - Update 'jd2_program=pyjd2' to contain your protein_mpnn conda environment path as such: 'jd2_program="/ssd1/home/safwand/miniconda3/envs/proteinmpnn_binder_design/bin/python dl_binder_design/include/silent_tools/pyjd2"'
        - Change 'silentls \$tmp_file | sed 's/_0001$//g' | silentrename \$tmp_file' to 'dl_binder_design/include/silent_tools/silentls \$tmp_file | sed 's/_0001$//g' | dl_binder_design/include/silent_tools/silentrename \$tmp_file'
- Now that everything is installed, you can update the binderDesignConfig.txt with your preferred paths for things and the paths to the conda environments you've now created. 
- You should now be set up!

# Run Instructions
## Binder Target
- To design a binder target using RosettaRemodel by updating an existing PDB with a certain motif, you must first generate a blueprint by running the following command: 'dl_binder_design/Rosetta/rosettaApp/rosetta.source.release-371/main/tools/remodel/getBlueprintFromCoords.pl -pdbfile [starting pdb] > [blueprint file]'
- Once you have a blueprint file, you can update individual residues with the type you'd like. For example:
    - '45 S . PIKAA R' updates the 45th residue (serine) to an arginine
    - '46 I . APOLAR' updates the 46th residue (isoleucine) to a non polar amino acid
- Once you've created such a file run Rosetta Remodel, via the following 'dl_binder_design/Rosetta/rosettaApp/rosetta.source.release-371/main/source/bin/remodel.linuxgccrelease -s YOUR_STARTING_PDB_TO_REMODEL.pdb -remodel:blueprint YOUR_BLUEPRINT_FILE -no_optH false -ex1 -ex2 -use_input_sc -linmem_ig 10 -num_trajectory 5 -save_top 5'
- This will result in a pdb being saved with the designed changes.

## Binder Design
- Now that everything is set up, the binder design is quite simple. Simply run the following command: 'python binderDesign.py -target YOUR_BINDER_TARGET.pdb -hotspot "[A43, A44, A45]" -bindersize "50-150" -chainrange "A1-64" -designcount 25 -mpnnsequences 2
    - You can update any of the flags to match your specifications.
