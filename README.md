# Brief description of scripts

- ***1_download_structure_with_sep_chain.py***: a) to download experimentally verified structures from [RCSB PDB](https://www.rcsb.org/) or modelled structures from [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/); b) to extract specific protein chain from pdb-file with the whole structure into separate pdb-file (only for experimental structures).

- ***2_del_ligands_pyv2.py***: a) to remove ligands from pdb-file with specific protein chain using [Chimera](http://www.cgl.ucsf.edu/chimera) tool.

- ***3_parse_pdb.py***: a) to extract "b-factor" values (one of the structural features); b) to map polypeptide sequence and structure positions.

- ***3_parse_pdb_with_cuts.py***: a) to extract "b-factor" values (one of the structural features); b) to map polypeptide sequence and structure positions; c) to map proteolytic cleavage sites from sequence into structure.

- ***4_create_dssp.py***: a) to generate dssp-files locally or remotely using pdb-files. 

- ***5_extract_features.py***: a) to extract feature information from DSSP files - initial type of secondary structure, solvent accessibility and others; b) to generate ad-hoc features - length of loop, type of secondary structure, terminal regions - based on information of initial type of secondary structure; c) to map structure information between DSSP files and PDB files.

- ***6_norm_data.py***: a) to apply normalisation of features for variables with float type; b) to generate dummy variables from secondary structure information. As a rule, we only apply normalisation within the protein structure chain. While creating of training dataset we applied two mode of normalisation: within the protein structure chain and within the whole dataset. 

- ***7_get_structural_score.py***: a) to predict structural scores of proteolytic cleavage sites using our structural model.

- ***ROC_AUC.py***: a) to generate ROC-curves and ROC AUC scores on testing dataset; b) to compare our results with ProCleave results.

- ***corr_proba.py***: a) to visualise PWM and structural scores on plot for specific protein substrate.

- ***decision_boundary.py***: a) to visualise decision boundary plot while getting total score.  

- ***dist_proba.py***: a) to visualise distribution of scores predicted.

- ***evaluate_StrModel.py***: a) to estimate our model on training dataset; b) to visualise performance score of our model.

- ***evaluate_features.py***: a) to estimate predictive performance of specific structural features; b) to visualise estimates.

- ***map_cuts_pyv2.py***: a) to visualise (map) proteolytic cleavage sites onto structures using Chimera.  

- ***map_scores_pyv2.py***: a) to visualise (map) proteolytic cleavage scores onto structures using Chimera.

- ***preprocess1_CutDB.py***: a) to aggregate information about proteolytic cleavage sites from CutDB and sequences of protein substrates.

- ***preprocess2_CutDB.py***: a) to present information about proteolytic cleavage sites and sequence for each protein substrate in the table form.

- ***save_model.py***: a) to save ML model.

- ***statistics*.py***: a) to get summary statistics - the number of unique protein substrate ID, the number of unique structure ID, the number of proteolytic cleavage sites, the number of proteases - on the different step of creating training dataset.

- 
