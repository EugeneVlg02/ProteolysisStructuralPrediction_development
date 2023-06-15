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
