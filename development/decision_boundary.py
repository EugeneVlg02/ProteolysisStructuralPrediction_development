import matplotlib.pyplot as plt
cm = plt.cm.RdBu

from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB

import pandas as pd
import numpy as np
np.random.seed(8)

''' Train data from E.coli experiment to get total score '''
TrainData_Ecoli = pd.read_csv("/data2/Protease/create_model_step/test_data_model/test_structures/TrainData_Ecoli.csv")

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

X_TrainData_Ecoli = TrainData_Ecoli[["MMP9_MMP25_PMAP", "LinearDA_Scikit_score"]].values
y_TrainData_Ecoli = TrainData_Ecoli["is_cut"].values

positive_indices = np.where(y_TrainData_Ecoli == 1)[0]
negative_indices = np.where(y_TrainData_Ecoli == 0)[0]

sample_indices = form_sample(positive_indices, negative_indices, "OneToOne")
X_TrainData_Ecoli = X_TrainData_Ecoli[sample_indices]
y_TrainData_Ecoli = y_TrainData_Ecoli[sample_indices]

clf = RandomForestClassifier(max_depth=5, n_estimators=10, max_features=1)
clf.fit(X_TrainData_Ecoli, y_TrainData_Ecoli)
disp = DecisionBoundaryDisplay.from_estimator(clf, X_TrainData_Ecoli, response_method="predict_proba", xlabel="PWM score", ylabel="Structural score", cmap=cm, alpha=0.8)
plt.ylim([-0.1, 1.1])
plt.savefig('decision_boundary_RF.jpg')