''' Import libraries '''
import _pickle as cPickle
import os
import sys
import glob
import time

import pandas as pd
import numpy as np
np.random.seed(8)

''' ML algorithms '''
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB

''' Set paths '''
script_path = os.getcwd()
main_path = os.path.split(script_path)[0]
model_path = os.path.join(main_path, "models")
if not os.path.exists(model_path):
    os.mkdir(model_path)

''' Sample '''
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

''' Train dataset for structural model '''
train_dataset_file = "/data2/Protease/create_dataset_step/datasets/CutDB_PDB/CutDB_PDB.csv"
train_df = pd.read_csv(train_dataset_file, dtype={"structure":str, "chain":str, "num_aa":str, "is_cut":int, "substrate":str, "SS_B":int, "SS_H":int, "SS_O":int})

X = train_df[["ACC_N_com", "bfac_N_chain", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_edge"]].values
y = train_df["is_cut"].values

positive_indices = np.where(y == 1)[0]
negative_indices = np.where(y == 0)[0]
sample_indices = form_sample(positive_indices, negative_indices, "OneToOne")

X_train = X[sample_indices]
y_train = y[sample_indices]

clf = LinearDiscriminantAnalysis()
clf.fit(X_train, y_train)

with open(os.path.join(model_path, 'structural_model.pkl'), 'wb') as file:
    cPickle.dump(clf, file)

''' Train data from E.coli experiment to get total score '''
TrainData_Ecoli = pd.read_csv("/data2/Protease/create_model_step/test_data_model/test_structures/TrainData_Ecoli.csv")

X_TrainData_Ecoli = TrainData_Ecoli[["MMP9_MMP25_PMAP", "LinearDA_Scikit_score"]].values
y_TrainData_Ecoli = TrainData_Ecoli["is_cut"].values

positive_indices = np.where(y_TrainData_Ecoli == 1)[0]
negative_indices = np.where(y_TrainData_Ecoli == 0)[0]
sample_indices = form_sample(positive_indices, negative_indices, "OneToOne")

X_TrainData_Ecoli = X_TrainData_Ecoli[sample_indices]
y_TrainData_Ecoli = y_TrainData_Ecoli[sample_indices]

clf = GaussianNB()
clf.fit(X_TrainData_Ecoli, y_TrainData_Ecoli)

with open(os.path.join(model_path, 'total_model.pkl'), 'wb') as file:
    cPickle.dump(clf, file)