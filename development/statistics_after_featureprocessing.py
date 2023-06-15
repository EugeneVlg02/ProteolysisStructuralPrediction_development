import pandas as pd

df = pd.read_csv("AFPDB_CutDB.csv")

AF = df.loc[df["structure"].str.contains("AF-")]
PDB = df.loc[df["structure"].str.contains('_')]

def get_unique(dataset, feature_name):
    values = set()
    for temp1 in dataset[feature_name].astype(str).unique():
        for temp2 in temp1.split('&'):
            for substrate in temp2.split('|'):
                if substrate != '-':
                    values.add(substrate)
    return values

PDB_substrates = get_unique(PDB, "SUBSTRATE_name")
AF_substrates = get_unique(AF, "SUBSTRATE_name")

PDB_proteases = get_unique(PDB, "MEROPS_CODE")
AF_proteases = get_unique(AF, "MEROPS_CODE")

print(AF["structure"].nunique())
print(len(AF_substrates))
print(AF["is_cut"].sum())
print(len(AF_proteases))
print()
print(PDB["structure"].nunique())
print(len(PDB_substrates))
print(PDB["is_cut"].sum())
print(len(PDB_proteases))