""" Import libraries """
import os
import time
import argparse
import _pickle as cPickle

import pandas as pd
import numpy as np
np.random.seed(8)

""" Set paths """
current_path = os.getcwd()
fitted_models_path = os.path.join(current_path, "fitted_models")
pwmData_path = os.path.join(current_path, "PWM")
script_path = os.path.join(current_path, "scripts")

parser = argparse.ArgumentParser(description='### Here should be the description! ###')
parser.add_argument("-input", help="The name of your PDB ID with separate chain, for example '4GAW_A'", type=str)
parser.add_argument("-model", help="The name of model used, for example 'GaussianNB'", type=str)
parser.add_argument("-protease", help="The name of protease used, for example 'MMP9'", type=str)
parser.add_argument("-pwm", help="The name of PWM used, for example 'PMAP'", type=str)

args = parser.parse_args()

structure_chain = args.input
chosen_model = args.model + ".pkl"
protease_name = args.protease
pwm_name = args.pwm

feature_path = os.path.join(current_path, structure_chain)

col_type = {"structure":str, "chain":str, "num_aa":str, "substrate":str}

features = ["bfac_N_chain", "SS_O", "SS_B", "SS_H", "len_loop_N_chain", "is_edge"]

def get_pwm_score(sequence, pwm):

    global_score = []
    for index_aa in range(len(sequence)):
        local_score = []
        if (index_aa == 0):
            seq_part = sequence[:index_aa+4]
            local_score.append(np.mean(pwm["P3"]))
            local_score.append(np.mean(pwm["P2"]))
            for aa, pos in zip(seq_part, ["P1", "P1'", "P2'", "P3'"]):
                if aa == "X":
                    aa = "M"
                local_score.append(pwm.loc[pwm["AA"] == aa, pos].values[0])
        elif index_aa == 1:
            seq_part = sequence[:index_aa+4]
            local_score.append(np.mean(pwm["P3"]))
            for aa, pos in zip(seq_part, ["P2", "P1", "P1'", "P2'", "P3'"]):
                if aa == "X":
                    aa = "M"
                local_score.append(pwm.loc[pwm["AA"] == aa, pos].values[0])
        elif index_aa == len(sequence)-3:
            seq_part = sequence[index_aa-2:]
            local_score.append(np.mean(pwm["P3'"]))
            for aa, pos in zip(seq_part, ["P3", "P2", "P1", "P1'", "P2'"]):
                if aa == "X":
                    aa = "M"
                local_score.append(pwm.loc[pwm["AA"] == aa, pos].values[0])
        elif index_aa == len(sequence)-2:
            seq_part = sequence[index_aa-2:]
            local_score.append(np.mean(pwm["P2'"]))
            local_score.append(np.mean(pwm["P3'"]))
            for aa, pos in zip(seq_part, ["P3", "P2", "P1", "P1'"]):
                if aa == "X":
                    aa = "M"
                local_score.append(pwm.loc[pwm["AA"] == aa, pos].values[0])
        elif index_aa == len(sequence)-1:
            continue
        else:
            seq_part = sequence[index_aa-2:index_aa+4]
            for aa, pos in zip(seq_part, ["P3", "P2", "P1", "P1'", "P2'", "P3'"]):
                if aa == "X":
                    aa = "M"
                local_score.append(pwm.loc[pwm["AA"] == aa, pos].values[0])

        #print(seq_part)
        #print(local_score)
        pwm_score = sum(local_score)
        global_score.append(round(pwm_score, 4))
    global_score.append(np.nan)
    return global_score

def get_max_min(pwm):
    list = []
    for p3 in pwm["P3"].values:
        for p2 in pwm["P2"].values:
            for p1 in pwm["P1"].values:
                for p1_ in pwm["P1'"].values:
                    for p2_ in pwm["P2'"].values:
                        for p3_ in pwm["P3'"].values:
                            list.append(p3 + p2 + p1 + p1_ + p2_ + p3_)
    print(max(list), min(list))

def get_structural_score(method, X_test):
    structural_scores = {}
    with open(os.path.join(fitted_models_path, method), 'rb') as file:
        clf = cPickle.load(file)

    structural_score = clf.predict_proba(X_test)[:, 1]
    structural_score_name = f"{method.split('.')[0]}_score"
    structural_scores[structural_score_name] = list(map(lambda x: round(x, 5), list(structural_score)))

    return pd.DataFrame(structural_scores)

def main():
    # Prepare PWM #
    pwm = pd.read_csv(os.path.join(pwmData_path, f"{protease_name}_{pwm_name}.txt"), sep='\t')

    # Work with every structure #
    feature_files = [i for i in os.listdir(feature_path) if (i == f'{structure_chain}.csv')]
    for index, feature_file in enumerate(feature_files):

        test_data = pd.read_csv(os.path.join(feature_path, feature_file), dtype=col_type)
        if "SS_B" not in list(test_data.columns):
            test_data["SS_B"] = 0
        X_test = test_data[features].values

        # Get PWM score #
        test_data[f"{protease_name}_{pwm_name}_score"] = get_pwm_score(test_data["AA"].values, pwm)

        # Get structural score #
        structural_scores_df = get_structural_score(chosen_model, X_test)
        test_data = pd.concat([test_data, structural_scores_df], axis=1)
        #print(test_data)
        # Save results #
        test_data.to_csv(os.path.join(feature_path, f"{feature_file.split('.')[0]}_{protease_name}_{pwm_name}_{chosen_model.split('.')[0]}.csv"), index=False)

""" Launch script """
if __name__ == "__main__":
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
