import os
import pandas as pd

script_path = os.getcwd()
basic_path = os.path.join(script_path, "basic")
derived_path = os.path.join(script_path, "derived")

''' Cuts info '''
cuts_df = pd.read_csv(os.path.join(basic_path, "cuts_CutDB.txt"), sep='\t', dtype={"PROTEIN_ID":str, "P1":str, "MEROPS_CODE":str})
agg_cuts_df = cuts_df.groupby("PROTEIN_ID").agg({"P1":'|'.join, "MEROPS_CODE":'|'.join}).reset_index()
agg_cuts_df["sequence"] = '-'

''' Seq info '''
with open(os.path.join(basic_path, "seqs_CutDB.fsa")) as file:
    data = file.read()

for protein in data.strip('\n').lstrip('>').split('>'):
    ID = protein.strip('\n').split('\n')[0]
    seq = ''.join(protein.strip('\n').split('\n')[1:])
    agg_cuts_df.loc[agg_cuts_df["PROTEIN_ID"] == ID, "sequence"] = seq

agg_cuts_df.to_csv(os.path.join(derived_path, "cuts_seqs_CutDB.csv"), sep=',', index=False)