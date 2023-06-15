import os
import time
import glob

import pandas as pd
from chimera import runCommand as run

''' Set paths '''
main_path = "/data2/Protease/TEST"
dataset_path = os.path.join(main_path, "datasets")
name_of_set = raw_input("Enter the name of structure set:\n{}\n".format(', '.join([i for i in os.listdir(dataset_path) if '.' not in i])))
StructureSet_path = os.path.join(dataset_path, name_of_set)
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

output_path = os.path.join("/data2/Protease/TEST/output", name_of_set)
if not os.path.exists(output_path):
    os.mkdir(output_path)
chimera_path = os.path.join(output_path, "chimera")
if not os.path.exists(chimera_path):
    os.mkdir(chimera_path)
save_path = os.path.join(chimera_path, "map_cuts")
if not os.path.exists(save_path):
        os.mkdir(save_path)

def main():
    num = 1
    for structure in structure_list:        
        structure_path = os.path.join(StructureSet_path, structure)
        feature_path = os.path.join(structure_path, "features")
        score_path = os.path.join(structure_path, "scores")
        
        score_files = glob.glob(os.path.join(score_path, '*.csv'))
        if len(score_files) == 0: continue
        feature_files = glob.glob(os.path.join(feature_path, '*_norm.csv'))
        if len(feature_files) == 0: continue
        
        for feature_file in feature_files:
        
            feature_df = pd.read_csv(feature_file, dtype={"structure":str, "chain":str, "num_AA":str, "AA":str, "is_cut":int})
            structure_name = feature_df["structure"].unique()[0] # AF-P12228-F1 or 1AUC_A #
            chain = feature_df["chain"].unique()[0]
            cuts_info = feature_df.loc[feature_df["is_cut"] == 1, ["num_AA", "MEROPS_CODE"]].values
            cuts = [":{0}.{1}".format(i[0], chain) for i in cuts_info]
            proteases = [i[1] for i in cuts_info]
            
            print "{0} --- {1}".format(num, structure_name)
            num += 1
            
            run("open {0}/{1}.pdb".format(structure_path, structure_name))
            run("color blue")
            run("color yellow {0}".format(' '.join(cuts)))
            for i in range(len(proteases)):
                run("labelopt resinfo {0}".format(proteases[i]))
                run("rlabel {0}".format(cuts[i]))
            run("save {0}/{1}.py".format(save_path, structure_name))
            run("del")
            
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print("Start: " + str(start) + '\n' + "Finish: " + str(finish) + '\n')