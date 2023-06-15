''' Import libraries '''
import os
import sys
import glob
import time

import pandas as pd
import numpy as np

''' Set paths '''
main_path = "/data2/Protease/TEST"
dataset_path = os.path.join(main_path, "datasets")
name_of_set = input("Enter the name of structure set:\n{}\n".format(', '.join([i for i in os.listdir(dataset_path) if '.' not in i])))
StructureSet_path = os.path.join(dataset_path, name_of_set)
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

get_common_dataset = input("Do you want to form common dataset? YES or NO\n")

def normalise(dataset, name_input_vector, name_output_vector):
    min_value = np.min(dataset[name_input_vector])
    max_value = np.max(dataset[name_input_vector])
    if max_value - min_value == 0:
        dataset[name_output_vector] = np.nan
    else:
        dataset[name_output_vector] = dataset[name_input_vector].apply(lambda x: round( (x - min_value) / (max_value - min_value), 3))

def main():
    warnings = []
    common_dataset = pd.DataFrame()
    
    num = 1
    for structure in structure_list:        
        structure_path = os.path.join(StructureSet_path, structure)
        feature_path = os.path.join(structure_path, "features")
        
        feature_files = glob.glob(os.path.join(feature_path, '*_final.csv'))
        for feature_file in feature_files:
            structure_name = '_'.join(feature_file.split('/')[-1].split('.')[0].split('_')[:-1])
            print(f"{num} --- {structure_name}")
            num += 1
            
            feature_df = pd.read_csv(feature_file)
            if ('AF-' in structure_name) and ('is_cut' in feature_df.columns):
                feature_df = feature_df.astype({"structure":str, "chain":str, "AA":str, "num_AA":str, "is_cut":int, "MEROPS_CODE":str, "SUBSTRATE_name":str, "init_SS_type":str, "SS_type":str, "bp1":int, "bp2":int, "ACC":float, "AF_score":float, "len_loop":int, "is_terminus":int})
            elif ('AF-' not in structure_name) and ('is_cut' in feature_df.columns):
                feature_df = feature_df.astype({"structure":str, "chain":str, "AA":str, "num_AA":str, "is_cut":int, "MEROPS_CODE":str, "SUBSTRATE_name":str, "init_SS_type":str, "SS_type":str, "bp1":int, "bp2":int, "ACC":float, "bfac":float, "len_loop":int, "is_terminus":int})
            elif ('AF-' in structure_name) and ('is_cut' not in feature_df.columns):
                feature_df = feature_df.astype({"structure":str, "chain":str, "AA":str, "num_AA":str, "init_SS_type":str, "SS_type":str, "bp1":int, "bp2":int, "ACC":float, "AF_score":float, "len_loop":int, "is_terminus":int})
            else:
                feature_df = feature_df.astype({"structure":str, "chain":str, "AA":str, "num_AA":str, "init_SS_type":str, "SS_type":str, "bp1":int, "bp2":int, "ACC":float, "bfac":float, "len_loop":int, "is_terminus":int})
            
            ''' ACC (structure) '''
            normalise(feature_df, "ACC", "ACC_N_chain")
            if True in np.unique(feature_df["ACC_N_chain"].isna()):
                warnings.append(f"{structure_name} - uncorrected normalisation of Solvent accessibility (ACC)!")
                feature_df["ACC_N_chain"].fillna(0, inplace=True)
                #continue
            ''' B-factor / AlphaFold prediction score '''
            if 'AF-' in structure_name:
                normalise(feature_df, "AF_score", "AF_score_N_chain")
                if True in np.unique(feature_df["AF_score_N_chain"].isna()):
                    warnings.append(f"{structure_name} - uncorrected normalisation of AlphaFold prediction score!")
                    feature_df["AF_score_N_chain"].fillna(0, inplace=True)
                    #continue
            else:
                normalise(feature_df, "bfac", "bfac_N_chain")
                if True in np.unique(feature_df["bfac_N_chain"].isna()):
                    warnings.append(f"{structure_name} - uncorrected normalisation of B-factor!")
                    feature_df["bfac_N_chain"].fillna(0, inplace=True)
                    #continue
            ''' Loop length '''
            normalise(feature_df, "len_loop", "len_loop_N_chain")
            if True in np.unique(feature_df["len_loop_N_chain"].isna()):
                warnings.append(f"{structure_name} - uncorrected normalisation of Loop length!")
                feature_df["len_loop_N_chain"].fillna(0, inplace=True)
                #continue

            ''' Secondary structure type coding '''
            feature_df = pd.concat([pd.get_dummies(feature_df, prefix='SS', columns=['SS_type']), feature_df["SS_type"]], axis=1)
            if 'SS_O' not in feature_df.columns:
                feature_df["SS_O"] = 0
            if 'SS_B' not in feature_df.columns:
                feature_df["SS_B"] = 0
            if 'SS_H' not in feature_df.columns:
                feature_df["SS_H"] = 0
            
            ''' Save results '''
            #feature_df.to_csv(os.path.join(feature_path, feature_file.split('.')[0] + '_norm.csv'), index=False)
            if get_common_dataset == "YES" and len(feature_df["is_cut"].unique()) == 2:
                common_dataset = pd.concat([common_dataset, feature_df])
    
    if len(common_dataset) != 0:
        ''' ACC (common) '''
        normalise(common_dataset, "ACC", "ACC_N_com")
        common_dataset.to_csv(os.path.join(dataset_path, f"{name_of_set}2.csv"), index=False)
        print(common_dataset["is_cut"].value_counts())
    
    if len(warnings) > 0:
        with open(os.path.join(StructureSet_path, f"{name_of_set}_UncorrectedNormalisation.txt"), 'w') as file:
            file.write('\n'.join(warnings))
                
''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")