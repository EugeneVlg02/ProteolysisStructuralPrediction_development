import os
import pandas as pd

script_path = os.getcwd()
derived_path = os.path.join(script_path, "derived")
blastp_path = "/data2/Protease/TEST/input/blastp"
name_of_set = input(f"Enter the name of structure set: {', '.join([i for i in os.listdir(blastp_path) if '.' not in i])}\n")
proteolytic_path = os.path.join(script_path, "proteolytic_events", name_of_set)

if not os.path.exists(proteolytic_path):
    os.mkdir(proteolytic_path)

blastp_df = pd.read_csv(os.path.join(blastp_path, name_of_set, f"{name_of_set}_blastpF.csv"), dtype={"query":str, "subject":str})
derived_df = pd.read_csv(os.path.join(derived_path, "cuts_seqs_CutDB.csv"), dtype={"PROTEIN_ID":str})

c = 1
No_structure_list = []
for item in derived_df.groupby("PROTEIN_ID"):
    PROTEIN_ID = item[0]
    print(f"{c} --- {PROTEIN_ID}")
    c += 1
    subset = item[1]
    seq = list(subset["sequence"].tolist()[0])
    len_seq = range(1, len(seq) + 1)
    
    if len(blastp_df.loc[blastp_df["query"] == PROTEIN_ID, "subject"].tolist()) == 0: 
        No_structure_list.append(PROTEIN_ID)
        continue
    structure = blastp_df.loc[blastp_df["query"] == PROTEIN_ID, "subject"].tolist()[0]
    
    substrate_df = pd.DataFrame({"SUBSTRATE_AA":seq,
                                "SUBSTRATE_num_AA":len_seq,
                                "is_cut":0,
                                "MEROPS_CODE":'-',
                                "SUBSTRATE_name":PROTEIN_ID,
                                "structure":structure
                                })
    
    for p1, code in zip(subset["P1"].tolist()[0].split('|'), subset["MEROPS_CODE"].tolist()[0].split('|')):
        if substrate_df.loc[substrate_df["SUBSTRATE_num_AA"] == int(p1), "is_cut"].tolist()[0] == 0:
            substrate_df.loc[substrate_df["SUBSTRATE_num_AA"] == int(p1), "is_cut"] = 1
            substrate_df.loc[substrate_df["SUBSTRATE_num_AA"] == int(p1), "MEROPS_CODE"] = code
        else:
            substrate_df.loc[substrate_df["SUBSTRATE_num_AA"] == int(p1), "MEROPS_CODE"] = substrate_df.loc[substrate_df["SUBSTRATE_num_AA"] == int(p1), "MEROPS_CODE"].values[0] + '|' + code
    substrate_df.to_csv(os.path.join(proteolytic_path, f"{PROTEIN_ID}.csv"), index=False)

with open(os.path.join(proteolytic_path, "NoBLASTpStructure.txt"), 'w') as file:
    file.write('\n'.join(No_structure_list))