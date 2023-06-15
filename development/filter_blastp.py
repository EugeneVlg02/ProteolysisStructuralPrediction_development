import pandas as pd
import os

script_path = os.getcwd()
dir_names = [i for i in os.listdir(script_path) if '.' not in i]
name_of_set = input(f"Enter the name of dataset: {', '.join(dir_names)}\n")
set_path = os.path.join(script_path, name_of_set)

name_columns = ["query", "subject", "%_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score", "qcovs"]
df = pd.read_csv(os.path.join(set_path, f"{name_of_set}_blastp.csv"), names=name_columns, dtype={"query":str, "subject":str})
df["subject"] = df["subject"].apply(lambda x: x.split(':')[1] if ':' in x else x.split('_')[0].upper() + '_' + x.split('_')[1])

result_df = pd.DataFrame()
for item in df.groupby("query"):
    q = item[0]
    subset = item[1]
    filter_subset = subset.loc[(subset["%_identity"] >= 90) & (subset["qcovs"] >= 67), ["query", "subject", "%_identity", "qcovs"]]
    if len(filter_subset) == 0:
        filter_subset = pd.DataFrame({"query":[q], "subject":'-', "%_identity":'-', "qcovs":'-'})
    else:
        filter_subset = filter_subset.iloc[0:1, :]
    result_df = pd.concat([result_df, filter_subset], ignore_index=True)

result_df.to_csv(os.path.join(set_path, f"{name_of_set}_blastpF.csv"), sep=',', index=False)

result_df = result_df.groupby("subject").agg({"query":'|'.join}).reset_index().rename(columns={"subject":"structure"})
result_df = result_df.loc[result_df["structure"] != '-']
result_df.to_csv(os.path.join(os.path.split(script_path)[0], f"{name_of_set.split('_')[0]}_structures_{name_of_set.split('_')[1]}.csv"), sep=',', index=False)