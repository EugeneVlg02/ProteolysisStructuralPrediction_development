import pandas as pd
import os

current_path = os.getcwd()
proteins_path = os.path.join(current_path, "AFPDB_CutDB")
structures = [i for i in os.listdir(proteins_path) if '.' not in i]

AF_structures = 0
AF_substrates = []
AF_cleavages = 0
AF_proteases = []

PDB_structures = 0
PDB_substrates = []
PDB_cleavages = 0
PDB_proteases = []

for structure in structures:
    structure_path = os.path.join(proteins_path, structure)
    feature_path = os.path.join(structure_path, "features")
    
    feature_files = [i for i in os.listdir(feature_path) if '.csv' in i and 'init' not in i and 'final' not in i]
    for feature_file in feature_files:
        structure_name = feature_file.split('.pdb')[0]
        if 'AF-' in structure_name:
            AF_structures += 1
        else:
            PDB_structures += 1
        
        df = pd.read_csv(os.path.join(feature_path, feature_file))
        
        if len(df["is_cut"].unique()) == 1:
            continue
        
        for substrate_name in df["SUBSTRATE_name"].astype(str).unique():
            for substrate_temp in substrate_name.split('|'):
                for substrate in substrate_temp.split('&'):
                    if 'AF-' in structure_name:
                        if substrate not in AF_substrates and substrate != '-':
                            AF_substrates.append(substrate)
                    else:
                        if substrate not in PDB_substrates and substrate != '-':
                            PDB_substrates.append(substrate)
        
        if 'AF-' in structure_name:
            AF_cleavages += df["is_cut"].sum()
        else:
            PDB_cleavages += df["is_cut"].sum()
        
        for merops_code in df["MEROPS_CODE"].astype(str).unique():
            for code_temp in merops_code.split('|'):
                for code in code_temp.split('&'):
                    if 'AF-' in structure_name:
                        if code not in AF_proteases and code != '-':
                            AF_proteases.append(code)
                    else:
                        if code not in PDB_proteases and code != '-':
                            PDB_proteases.append(code)

print(AF_structures)
print(len(AF_substrates))
print(AF_cleavages)
print(len(AF_proteases))
print()
print(PDB_structures)
print(len(PDB_substrates))
print(PDB_cleavages)
print(len(PDB_proteases))
print()
print(PDB_substrates)       