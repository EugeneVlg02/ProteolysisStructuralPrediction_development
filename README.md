
The tool for predicting probabilities of proteolytic events in proteins. It based on ML algorithms and data from CutDB database.

## Scripts. Brief Description.

- ***run.py*** is the main script for predicting scores of proteolytic cleavage sites. It runs the scripts located in the "scripts" directory and outputs log information in response.

- ***1_download_structure_with_sep_chain.py*** is the script for downloading structure data and extracting separate protein chain data from it. Structure will be downloading from [PDB database.](https://www.rcsb.org/)

- ***2_del_ligands_pyv2.py*** is the script for removing ligands from structure data using [Chimera tool.](http://www.cgl.ucsf.edu/chimera)

- ***3_parse_pdb.py*** is the script for mapping proteolytic cleavage sites onto the structure and extracting feature information of b-factor.

- ***4_create_dssp.py*** is the script for generating DSSP files from PDB files of structure data.

- ***5_extract_features.py*** is the script for extracting feature information from DSSP files and generating ad-hoc feature information as well as for mapping structure information from DSSP files and PDB files.

- ***6_norm_data.py*** is the script for normalisation of features with float types and generating dummy variables from secondary structure information.

- ***7_get_structural_score.py*** is the script for predicting scores of proteolytic cleavage sites using our structural model located in the "models" directory.

- ***map_scores_pyv2.py*** is the optional script for visualisation of proteolytic cleavage site scores predicted.
