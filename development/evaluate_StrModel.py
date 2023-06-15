import pandas as pd
import numpy as np
np.random.seed(8)

from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

from xgboost.sklearn import XGBClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression

import os
import time

import plotly.graph_objects as go
import plotly.express as px

script_path = os.getcwd()
datasets = {"PDB_CutDB":"/data2/Protease/create_dataset_step/datasets/CutDB_PDB/CutDB_PDB.csv",
            "AFPDB_CutDB":"/data2/Protease/TEST/datasets/AFPDB_CutDB.csv"
            }
output_path = os.path.join("/data2/Protease/TEST/output", "Structural_model")
if not os.path.exists(output_path):
    os.mkdir(output_path)

methods = {
           'DecisionTree':DecisionTreeClassifier(), #636efa
           'RandomForest':RandomForestClassifier(), #636efa
           'XGBoost':XGBClassifier(), #00cc96
		   'LinearDA':LinearDiscriminantAnalysis(), #aa64fa
		   'QuadraticDA':QuadraticDiscriminantAnalysis(), #ffa25a
           'GaussianNB':GaussianNB(), #19d3f4
		   'KNeighbors':KNeighborsClassifier(), #ff6691
           'SVM':SVC(kernel='rbf', probability=True), #b5e784
		   'LogisticRegression':LogisticRegression() #ff97ff
           }
           
features = ["ACC_N_com", "bfac_N_chain", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_edge"]
approaches = ["OneToOne", "OneToTwo", "OneToFive", "OneToTen", "OneToFifty", "Dataset"]

# Auxiliary function 1 #
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

# Auxiliary function 2 #
def cv_auc(clf, X_cv, y_cv, X_ost, y_ost, cv):
    
    common_AUC = np.array([], dtype="float")
    
    for train_index, test_index in cv.split(X_cv, y_cv):
        X_train = X_cv[train_index]
        y_train = y_cv[train_index]
        X_part2_test = X_cv[test_index]
        y_part2_test = y_cv[test_index]
        X_test = np.concatenate([X_part2_test, X_ost])
        y_test = np.concatenate([y_part2_test, y_ost])
        
        clf.fit(X_train, y_train)
        y_proba = clf.predict_proba(X_test)[:, 1]
        
        common_AUC_score = roc_auc_score(y_test, y_proba)
        common_AUC = np.append(common_AUC, common_AUC_score)
        
    return common_AUC

def main(dataset, name_dataset, features=features, approaches=approaches, methods=methods):
    df = pd.read_csv(dataset)
    #3.3); 3.6) - 
    #df["AF_score_N_chain"] = df["AF_score_N_chain"].fillna(round(df["AF_score_N_chain"].mean(), 3))
    
    X = df[features].values
    y = df["is_cut"].values
    
    skf = StratifiedKFold(n_splits=10)
    
    indices = np.arange(0, len(df), 1)
    positive_indices = np.where(y == 1)[0]
    negative_indices = np.where(y == 0)[0]
    
    result_df = pd.DataFrame()
    for approach in approaches:
        print(f"Approach --- {approach}")
        cv_indices = form_sample(positive_indices, negative_indices, approach)
        
        ost_indices = np.array([], dtype='int')
        for i in indices:
            if i not in cv_indices:
                ost_indices = np.append(ost_indices, i)
                
        X_cv = X[cv_indices]
        y_cv = y[cv_indices]
    
        X_part1_test = X[ost_indices]
        y_part1_test = y[ost_indices]
              
        for method in methods:
            clf = methods[method]
            print(f"Method: {method}")
            common_AUC = cv_auc(clf, X_cv, y_cv, X_part1_test, y_part1_test, skf)
            temp_df = pd.DataFrame({"Approach":[approach]*10,
                                    "Method":[method]*10,
                                    "AUC":common_AUC}
            )
            result_df = pd.concat([result_df, temp_df], ignore_index=True)
    if (len(approaches) == 1) and (len(methods) == 1):        
        result_df.to_csv(os.path.join(output_path, f"Scikit_{name_dataset}_Evaluating{approaches[0]}{list(methods.keys())[0]}.csv"), sep=',', index=False)
    elif (len(approaches) == 1) and (len(methods) != 1):
        result_df.to_csv(os.path.join(output_path, f"Scikit_{name_dataset}_Evaluating{approaches[0]}.csv"), sep=',', index=False)
    elif (len(approaches) != 1) and (len(methods) == 1):
        result_df.to_csv(os.path.join(output_path, f"Scikit_{name_dataset}_Evaluating{list(methods.keys())[0]}.csv"), sep=',', index=False)
    else:
        result_df.to_csv(os.path.join(output_path, f"Scikit_{name_dataset}_Evaluating.csv"), sep=',', index=False)

def common_barplot():

    approaches = {"Dataset":"Dataset", "OneToOne":"1:1", "OneToTwo":"1:2", "OneToFive":"1:5", "OneToTen":"1:10", "OneToFifty":"1:50"}
    
    plot_df = pd.read_csv(os.path.join(output_path, "Scikit_PDB_CutDB_Evaluating.csv"))
    plot_df["Approach"] = plot_df["Approach"].replace(approaches)
    plot_df["Approach"] = pd.Categorical(plot_df["Approach"], ["1:1", "1:2", "1:5", "1:10", "1:50", "Dataset"])
    plot_df["Method"] = pd.Categorical(plot_df["Method"], methods.keys())
    
    statistics_plot_df = plot_df.groupby(["Approach", "Method"]).agg(AUC_mean=("AUC", np.mean), AUC_std=("AUC", np.std)).reset_index()
    #print(statistics_plot_df)
    #input()
    ''' Common plot '''
    fig = px.bar(statistics_plot_df, x = "Approach", y = "AUC_mean", color = "Method", error_y = "AUC_std", barmode="group")
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=False, gridwidth=1, gridcolor='lightgray')
    fig.update_yaxes(range=[0, 1], dtick=0.1, showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray', ticks="outside")
    fig.update_layout(title=dict(text="<b>CutDB PDB</b>", y=0.95, x=0.465, xanchor="center", yanchor="top"),
                xaxis_title="Approach",
                yaxis_title="AUC score",
                paper_bgcolor="white",
                plot_bgcolor="white",
                bargroupgap=0.2
                )
    #fig.write_html(os.path.join(output_path, f"Scikit_PDB_CutDB_EvaluatingCommon.html"))
    fig.write_image(os.path.join(output_path, f"Scikit_PDB_CutDB_EvaluatingCommon.jpg"), width = 650, height = 550, scale = 4)
    
    ''' 1:1 plot '''
    subset_statistics_plot_df = statistics_plot_df.loc[statistics_plot_df["Approach"] == "1:1"]
    fig = px.bar(subset_statistics_plot_df, x = "Approach", y = "AUC_mean", color = "Method", error_y = "AUC_std", barmode="group", width=750)
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=False, gridwidth=1, gridcolor='lightgray')
    fig.update_yaxes(range=[0, 1], dtick=0.1, showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray', ticks="outside")
    fig.update_layout(title=dict(text="<b>CutDB PDB (1:1)</b>", y=0.95, x=0.465, xanchor="center", yanchor="top"),
                xaxis_title="Approach",
                yaxis_title="AUC score",
                paper_bgcolor="white",
                plot_bgcolor="white",
                bargap=0.05
                )
                
    fig.update_traces(width=[0.07]*9)
    #fig.write_html(os.path.join(output_path, f"Scikit_PDB_CutDB_EvaluatingOneToOne.html"))
    fig.write_image(os.path.join(output_path, f"Scikit_PDB_CutDB_EvaluatingOneToOne.jpg"), width = 650, height = 550, scale = 4)
    
    ''' Only LinearDA '''
    subset_statistics_plot_df = statistics_plot_df.loc[statistics_plot_df["Method"] == "LinearDA"]
    fig = px.bar(subset_statistics_plot_df, x = "Approach", y = "AUC_mean", color = "Method", error_y = "AUC_std", barmode="group", width=650)
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=False, gridwidth=1, gridcolor='lightgray')
    fig.update_yaxes(range=[0, 1], dtick=0.1, showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray', ticks="outside")
    fig.update_layout(title=dict(text="<b>CutDB PDB (1:1)</b>", y=0.95, x=0.465, xanchor="center", yanchor="top"),
                xaxis_title="Approach",
                yaxis_title="AUC score",
                paper_bgcolor="white",
                plot_bgcolor="white",
                #bargap=0.001,
                #bargroupgap=0.5
                )
    fig.update_traces(marker=dict(color="#aa64fa"), width=[0.5]*6)
    
    #fig.write_html(os.path.join(output_path, f"Scikit_PDB_CutDB_EvaluatingLinearDA.html"))
    #fig.write_image(os.path.join(output_path, f"Scikit_PDB_CutDB_EvaluatingLinearDA.jpg"), width = 650, height = 550, scale = 4)
    
    
def PDB_AFPDB_comparison_plot():
    approaches = {"OneToOne":"1:1"}
    
    plot_PDB_df = pd.read_csv(os.path.join(output_path, "Scikit_PDB_CutDB_EvaluatingOneToOneXGBoost2.csv"))
    plot_PDB_df["Approach"] = plot_PDB_df["Approach"].replace(approaches)
    plot_PDB_df["Type"] = "PDB"
    plot_AFPDB_df = pd.read_csv(os.path.join(output_path, "Scikit_AFPDB_CutDB_EvaluatingOneToOneXGBoost2.csv"))
    plot_AFPDB_df["Approach"] = plot_AFPDB_df["Approach"].replace(approaches)
    plot_AFPDB_df["Type"] = "AlphaFold + PDB"
    concat_plot_df = pd.concat([plot_PDB_df, plot_AFPDB_df], ignore_index=True)
    print(concat_plot_df)
    print(concat_plot_df.groupby("Type").agg({"AUC":"median"}))
    input()
    fig = px.box(concat_plot_df, x="Type", y="AUC", width=650)
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=False, gridwidth=1, gridcolor='lightgray')
    fig.update_yaxes(range=[0, 1], dtick=0.1, showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray', ticks="outside")
    fig.update_layout(title=dict(text="<b>XGBoost (1:1)</b>", y=0.95, x=0.465, xanchor="center", yanchor="top"),
                xaxis_title="Dataset",
                yaxis_title="AUC score",
                paper_bgcolor="white",
                plot_bgcolor="white",
                boxgap=0.5
                )
    
    fig.update_traces(marker=dict(color="#00cc96"))

    #fig.write_html(os.path.join(output_path, f"Scikit_PDB_AFPDB_CutDB_EvaluatingOneToOneXGBoost.html"))
    #fig.write_image(os.path.join(output_path, f"Scikit_PDB_AFPDB_CutDB_EvaluatingOneToOneXGBoost2.jpg"), width = 650, height = 550, scale = 4)
    
''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    #1) main()
    #2) 
    common_barplot()
    
    ''' LinearDA '''
    #3.0) main(datasets["PDB_CutDB"], "PDB_CutDB", ["ACC_N_com", "bfac_N_chain", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_edge"], ["OneToOne"], {"LinearDA":LinearDiscriminantAnalysis()})
    #3.1) main(datasets["PDB_CutDB"], "PDB_CutDB", ["ACC_N_com", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_edge"], ["OneToOne"], {"LinearDA":LinearDiscriminantAnalysis()})
    #3.2) main(datasets["AFPDB_CutDB"], "AFPDB_CutDB", ["ACC_N_com", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_terminus"], ["OneToOne"], {"LinearDA":LinearDiscriminantAnalysis()})
    #3.3) main(datasets["AFPDB_CutDB"], "AFPDB_CutDB", ["ACC_N_com", "AF_score_N_chain", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_terminus"], ["OneToOne"], {"LinearDA":LinearDiscriminantAnalysis()})
    ''' XGBoost '''
    #3.4) main(datasets["PDB_CutDB"], "PDB_CutDB", ["ACC_N_com", "bfac_N_chain", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_edge"], ["OneToOne"], {"XGBoost":XGBClassifier()})
    #3.5) main(datasets["PDB_CutDB"], "PDB_CutDB", ["ACC_N_com", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_edge"], ["OneToOne"], {"XGBoost":XGBClassifier()})
    #3.6) main(datasets["AFPDB_CutDB"], "AFPDB_CutDB", ["ACC_N_com", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_terminus"], ["OneToOne"], {"XGBoost":XGBClassifier()})
    #3.7) main(datasets["AFPDB_CutDB"], "AFPDB_CutDB", ["ACC_N_com", "AF_score_N_chain", "SS_B", "SS_H", "SS_O", "len_loop_N_chain", "is_terminus"], ["OneToOne"], {"XGBoost":XGBClassifier()})
    
    #4) 
    PDB_AFPDB_comparison_plot() 
    
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")