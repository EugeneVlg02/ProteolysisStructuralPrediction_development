import pandas as pd
import os

basic_df = pd.read_csv("/data2/Protease/TEST/input/pe_info/derived/cuts_seqs_CutDB.csv", dtype={"PROTEIN_ID":str})
basic_df = basic_df.iloc[:,:-1]
basic_df["P1_count"] = basic_df["P1"].apply(lambda x: len(set(x.split('|'))))
basic_df["MEROPS_CODE_unique"] = basic_df["MEROPS_CODE"].apply(lambda x: set(x.split('|')))
print(basic_df)
structures_df = pd.read_csv("/data2/Protease/TEST/input/AFPDB_structures_CutDB.csv")
print(structures_df)

current_path = os.getcwd()
proteins_path = os.path.join("AFPDB_CutDB")
structures = [i for i in os.listdir(proteins_path) if '.' not in i]

AF_structures = 0
AF_substrates = 0
AF_cleavages = 0
AF_proteases = set()

PDB_structures = 0
PDB_substrates = 0
PDB_cleavages = 0
PDB_proteases = set()

for structure in structures:
    structure_path = os.path.join(proteins_path, structure)
    
    if 'AF-' in structure:
        structure_files = [i for i in os.listdir(structure_path) if '.pdb' in i]
    else:
        structure_files = [i for i in os.listdir(structure_path) if '.pdb' in i and '_' in i]
    
    for structure_file in structure_files:
        structure_name = structure_file.split('.pdb')[0]
        
        if 'AF-' in structure_name:
            AF_structures += 1
        else:
            PDB_structures += 1
        
        substrates = structures_df.loc[structures_df["structure"] == structure_name, "query"].tolist()[0]
        for substrate in substrates.split('|'):
            if 'AF-' in structure_name:
                AF_substrates += 1
            else:
                PDB_substrates += 1
            
            basic_subset = basic_df.loc[basic_df["PROTEIN_ID"] == substrate]
            if 'AF-' in structure_name:
                AF_cleavages += basic_subset["P1_count"].values[0]
                for protease_set in basic_subset["MEROPS_CODE_unique"].values:
                    for protease in protease_set:
                        AF_proteases.add(protease)
            else:
                PDB_cleavages += basic_subset["P1_count"].values[0]
                for protease_set in basic_subset["MEROPS_CODE_unique"].values:
                    for protease in protease_set:
                        PDB_proteases.add(protease)
print(AF_structures)
print(AF_substrates)
print(AF_cleavages)
print(len(AF_proteases))
print()
print(PDB_structures)
print(PDB_substrates)
print(PDB_cleavages)
print(len(PDB_proteases))
print()