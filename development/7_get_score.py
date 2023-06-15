''' Import libraries '''
import os
import sys
import glob
import time

import pandas as pd
import numpy as np
np.random.seed(8)

''' ML algorithms '''
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

methods = {
           'LinearDA':LinearDiscriminantAnalysis(),
           }


''' Set paths '''
main_path = "/data2/Protease/TEST"
dataset_path = os.path.join(main_path, "datasets")
name_of_set = input("Enter the name of structure set:\n{}\n".format(', '.join([i for i in os.listdir(dataset_path) if '.' not in i])))
StructureSet_path = os.path.join(dataset_path, name_of_set)
structure_list = [i for i in os.listdir(StructureSet_path) if '.' not in i]

''' PWM and proteases '''
PWM_path = "/data2/Protease/PWM/MEROPS"
#name_of_protease = input("Enter the name of protease used:\n{}\n".format(', '.join([i.split('_')[0] for i in os.listdir(PWM_path) if '_PWM' in i])))
#PWM = pd.read_csv(os.path.join(PWM_path, f"{name_of_protease}_PWM.txt"), sep='\t')

''' Features '''
features = {"ACC_N_chain":float, "bfac_N_chain":float, "SS_B":int, "SS_H":int, "SS_O":int, "len_loop_N_chain":float, "is_terminus":int}

''' Train dataset '''
train_dataset_file = "/data2/Protease/create_dataset_step/datasets/CutDB_PDB/CutDB_PDB.csv"
train_df = pd.read_csv(train_dataset_file, dtype={"structure":str, "chain":str, "num_aa":str, "is_cut":int, "substrate":str, "SS_B":int, "SS_H":int, "SS_O":int})

train_df["structure_chain"] = train_df["structure"] + '_' + train_df["chain"]
train_structures_list = list(train_df["structure_chain"].unique())

X = train_df[["ACC_N_com", "bfac_N_chain", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_edge"]].values
y = train_df["is_cut"].values
positive_indices = np.where(y == 1)[0]
negative_indices = np.where(y == 0)[0]

''' Train sample '''
def form_sample(p_indices, n_indices, ratio):
    
    choose_ratio = {'OneToOne':len(p_indices),
                    'OneToTwo':len(p_indices)*2,
                    'OneToFive':len(p_indices)*5,
     			    'OneToTen':len(p_indices)*10,
                    'OneToFifty':len(p_indices)*50,
                    'Dataset':len(n_indices)}
    np.random.shuffle(n_indices)

    part_n_indices = n_indices[:choose_ratio[ratio]]
    sample_indices = np.concatenate([p_indices, part_n_indices])
    np.random.shuffle(sample_indices)

    return sample_indices

sample_indices = form_sample(positive_indices, negative_indices, "OneToOne")
X_train = X[sample_indices]
y_train = y[sample_indices]

''' Train data from E.coli experiment to get total score '''
TrainData_Ecoli = pd.read_csv("/data2/Protease/create_model_step/test_data_model/test_structures/TrainData_Ecoli.csv")
train_Ecoli_structures_list = list(TrainData_Ecoli["structure_chain"].unique())
X_TrainData_Ecoli = TrainData_Ecoli[["MMP9_MMP25_PMAP", "LinearDA_Scikit_score"]].values
y_TrainData_Ecoli = TrainData_Ecoli["is_cut"].values

positive_indices = np.where(y_TrainData_Ecoli == 1)[0]
negative_indices = np.where(y_TrainData_Ecoli == 0)[0]

sample_indices = form_sample(positive_indices, negative_indices, "OneToOne")
X_TrainData_Ecoli = X_TrainData_Ecoli[sample_indices]
y_TrainData_Ecoli = y_TrainData_Ecoli[sample_indices]

''' PWM score function '''
def get_PWM_score(sequence, PWM):
    
    processed_sequence = f"--{sequence}--"
    PWM_scores = []
    for i in range(2, len(sequence) + 1):
        seq_frame = processed_sequence[i - 2: i + 4]
        local_score = 0
        for aa, pos in zip(seq_frame, ["P3", "P2", "P1", "P1'", "P2'", "P3'"]):
            if aa == '-': continue
            local_score += PWM.loc[PWM["AA"] == aa, pos].values[0]
        PWM_scores.append(local_score)  
    
    PWM_scores.append(np.nan)
    return PWM_scores

def main():
    
    num = 1
    for structure in structure_list:        
        structure_path = os.path.join(StructureSet_path, structure)
        feature_path = os.path.join(structure_path, "features")
        
        feature_files = glob.glob(os.path.join(feature_path, '*_norm.csv'))
        for feature_file in feature_files:
            structure_name = '_'.join(feature_file.split('/')[-1].split('.')[0].split('_')[:-2])
            print(f"{num} --- {structure_name}")
            num += 1
            
            feature_df = pd.read_csv(feature_file, dtype=features)
            if len(feature_df["is_cut"].unique()) != 2:
                print(f"{structure_name} without cuts!")
                continue
            if feature_df["structure"].unique()[0] in train_structures_list: 
                print(f"{structure_name} in training set!")
                continue
            if feature_df["structure"].unique()[0] in train_Ecoli_structures_list: 
                print(f"{structure_name} in Ecoli set!")
                continue
            
            score_path = os.path.join(structure_path, "scores")
            if not os.path.exists(score_path):
                os.mkdir(score_path)
            
            ''' Structural score '''
            
            X_test = feature_df[features.keys()].values
            for method in methods:
                clf = methods[method]
                clf.fit(X_train, y_train)
                structural_scores = clf.predict_proba(X_test)[:, 1]
                
                structural_score_col = f"{method}_Scikit"
                feature_df[structural_score_col] = structural_scores
                feature_df[structural_score_col] = feature_df[structural_score_col].map(lambda x: round(x, 4))
            
            ''' PWM (Sequence) score '''
            sequence = ''.join(feature_df["AA"].tolist())
            if 'X' in sequence:
                print('Unusual amino acid in sequence')
                continue
                #sys.exit()
            
            names_of_protease = [i for i in feature_df["MEROPS_CODE"].unique() if i != '-']
            
            for name_of_protease in names_of_protease:
                PWM = pd.read_csv(os.path.join(PWM_path, f"{name_of_protease}_PWM.txt"), sep='\t')
 
                PWM_scores = get_PWM_score(sequence, PWM)
                PWM_score_col = f"{name_of_protease}_PWM"
                
                copy_df = feature_df.copy()
                copy_df[PWM_score_col] = PWM_scores    
                copy_df = copy_df.iloc[:-1, :]
                
                ''' Only one protease '''
                copy_df.loc[~((copy_df["MEROPS_CODE"].str.contains(name_of_protease)) | (copy_df["MEROPS_CODE"].str.contains('-'))), "is_cut"] = 0
                
                ''' Total score '''
                score_cols = [i for i in copy_df.columns if 'PWM' in i] + [i for i in copy_df.columns if 'Scikit' in i]
                scores_array = copy_df[score_cols].values
               
                #total_clf = {
                #            "RandomForest":RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1),
                #            "AdaBoost":AdaBoostClassifier(),
                #            "GaussianNaiveBayes":GaussianNB(),
                #            "QuadraticDA":QuadraticDiscriminantAnalysis()
                #            }
                #for name_clf in total_clf:
                #    clf = total_clf[name_clf]
                #    clf.fit(X_TrainData_Ecoli, y_TrainData_Ecoli)
                #    total_score = clf.predict_proba(scores_array)[:, 1]
                #    copy_df[f"{name_clf}_Total_score"] = total_score
                #    copy_df[f"{name_clf}_Total_score"] = copy_df[f"{name_clf}_Total_score"].map(lambda x: round(x, 4))
                
                ''' Save results '''
                #copy_df.to_csv(os.path.join(score_path, f"{structure_name}_{name_of_protease}.csv"), index=False)

''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")