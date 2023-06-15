import pandas as pd

basic_df = pd.read_csv("/data2/Protease/TEST/input/pe_info/derived/cuts_seqs_CutDB.csv", dtype={"PROTEIN_ID":str})
basic_df = basic_df.iloc[:,:-1]
basic_df["P1_count"] = basic_df["P1"].apply(lambda x: len(set(x.split('|'))))
basic_df["MEROPS_CODE_unique"] = basic_df["MEROPS_CODE"].apply(lambda x: set(x.split('|')))

def get_statistics(dataset):
    c_substrates = 0
    c_cleavages = 0
    proteases = set()
    
    for substrates in dataset["query"].unique():
        substrates = substrates.split('|')
        for substrate in substrates:
            c_substrates += 1
            
            basic_subset = basic_df.loc[basic_df["PROTEIN_ID"] == substrate]
            c_cleavages += basic_subset["P1_count"].values[0]
            for protease_set in basic_subset["MEROPS_CODE_unique"].values:
                for protease in protease_set:
                    proteases.add(protease)
            
    
    print(c_substrates)
    print(c_cleavages)
    print(len(proteases))

            

df = pd.read_csv("AFPDB_structures_CutDB.csv")
PDB = df.loc[df["structure"].str.contains('_')]
AF = df.loc[df["structure"].str.contains('AF-')]

print(PDB)
print(AF)
get_statistics(PDB)
get_statistics(AF)
