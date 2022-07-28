""" Import libraries """
import os
import time

import pandas as pd
import numpy as np

""" Set paths """
current_path = os.getcwd()
script_path = os.path.join(current_path, "scripts")

col_type = {"structure":str,
            "chain":str,
            "num_aa":str,
            "substrate":str}

""" Set auxiliary functions """
def normalise(dataset, name_vector, name_normalised_vector):
    min_value = np.min(dataset[name_vector])
    max_value = np.max(dataset[name_vector])
    dataset[name_normalised_vector] = dataset[name_vector].apply(lambda x: round( (x - min_value) / (max_value - min_value) , 3))

""" Set main function """
def main():
    feature_files = [i for i in os.listdir(script_path) if '.csv' in i]
    for feature_file in feature_files:
        df = pd.read_csv(os.path.join(script_path, feature_file), dtype=col_type)

        # ACC bfac len_loop dist_loop #
        normalise(df, "bfac", "bfac_N_chain")
        normalise(df, "len_loop", "len_loop_N_chain")
        """
        normalise(df_chain, "ACC", "ACC_N_chain")
        """
        # ACC bfac len_loop  #
        """
        normalise(new_df_structure, "ACC", "ACC_N_str")
        normalise(new_df_structure, "len_loop", "len_loop_N_str")
        """

        ''' COMMON NORMALISATION '''
        # ACC bfac len_loop #
        normalise(df, "ACC", "ACC_N_com")
        #normalise(df, "len_loop", "len_loop_N_com")

        ''' CHECK IN NAN AND INF '''
        #print("Check bfac")
        if True in np.unique(df["bfac_N_chain"].isna()):
            #bfac_Nan_df = new_df.loc[new_df["bfac_N_chain"].isna(), ["structure", "chain", "num_aa", "AA", "SS_type", "bfac", "bfac_N_chain"]]
            #input("bfac - NA\nPlease, press on 'Enter' ...\n")
            df["bfac_N_chain"] = df["bfac_N_chain"].fillna(0.0)
        if True in np.unique(np.isinf(df["bfac_N_chain"].values)):
            #input("bfac - inf\nPlease, press on 'Enter' ...\n")
            df.loc[np.isinf(df["bfac_N_chain"].values) == True, "bfac_N_chain"] = 0.0

        #print("Check ACC")
        if True in np.unique(df["ACC_N_com"].isna()):
            #input("ACC - NA\nPlease, press on 'Enter' ...\n")
            df["ACC_N_com"] = df["ACC_N_com"].fillna(0.0)
        if True in np.unique(np.isinf(df["ACC_N_com"].values)):
            #input("ACC - inf\nPlease, press on 'Enter' ...\n")
            df.loc[np.isinf(df["ACC_N_com"].values) == True, "ACC_N_com"] = 0.0

        #print("Check len_loop")
        if True in np.unique(df["len_loop_N_chain"].isna()):
            #LenLoop_Nan_df = new_df.loc[new_df["len_loop_N_chain"].isna(), ["structure", "chain", "num_aa", "AA", "len_loop", "len_loop_N_chain"]]
            #input("len_loop - NA\nPlease, press on 'Enter' ...\n")
            df["len_loop_N_chain"] = df["len_loop_N_chain"].fillna(0.0)
        if True in np.unique(np.isinf(df["len_loop_N_chain"].values)):
            #input("len_loop - inf\nPlease, press on 'Enter' ...\n")
            df.loc[np.isinf(df["len_loop_N_chain"].values) == True, "len_loop_N_chain"] = 0.0
        '''
        if True in np.unique(df["len_loop_N_com"].isna()):
            #LenLoop_Nan_df = new_df.loc[new_df["len_loop_N_chain"].isna(), ["structure", "chain", "num_aa", "AA", "len_loop", "len_loop_N_chain"]]
            input("len_loop - NA\nPlease, press on 'Enter' ...\n")
            df["len_loop_N_com"] = df["len_loop_N_com"].fillna(0.0)
        if True in np.unique(np.isinf(df["len_loop_N_com"].values)):
            input("len_loop - inf\nPlease, press on 'Enter' ...\n")
            df.loc[np.isinf(df["len_loop_N_com"].values) == True, "len_loop_N_com"] = 0.0
        '''

        ''' ONE-HOT ENCODING '''
        df = pd.concat([pd.get_dummies(data=df, prefix="SS", columns=["SS_type"]), df["SS_type"]], axis=1)

        ''' SAVING '''
        columns = list(df.columns)
        df = df[columns[:4] + [columns[-1]] + columns[4:-1]]

        feature_path = os.path.join(current_path, feature_file.split('.')[0])
        if not os.path.exists(feature_path):
            os.mkdir(feature_path)
        df.to_csv(os.path.join(feature_path, feature_file), index=False)
        os.remove(os.path.join(script_path, feature_file))

""" Launch script """
if __name__ == '__main__':
    start = time.ctime()
    main()
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")
