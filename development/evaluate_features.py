import pandas as pd
import numpy as np
np.random.seed(8)

from sklearn.metrics import roc_auc_score

import os
import time

import plotly.graph_objects as go
import plotly.express as px

script_path = os.getcwd()
dataset = "/data2/Protease/TEST/datasets/AFPDB_CutDB.csv"
output_path = os.path.join("/data2/Protease/TEST/output", "Structural_model")
if not os.path.exists(output_path):
    os.mkdir(output_path)

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
    
def main():
    df = pd.read_csv(dataset)
    AF_df = df.loc[df["structure"].str.contains("AF-")]
    str_df = df.loc[~df["structure"].str.contains("AF-")]

    features = {"AF_score_N_chain":"AlphaFold confidence", "ACC_N_com":"Solvent accessibility", "len_loop_N_chain":"Loop length"}
    
    ''' Common '''
    plot_df = pd.DataFrame()
    for feature in features:    
        feature_name = features[feature]
        
        if feature_name == "B-factor":
            y_true = str_df["is_cut"].tolist()
            y_score = str_df[feature].tolist()
        else:
            y_true = AF_df["is_cut"].tolist()
            y_score = AF_df[feature].tolist()
        
        AUC = round(roc_auc_score(y_true, y_score), 3)
        if AUC < 0.5:
            AUC = 1 - AUC
            
        temp_df = pd.DataFrame({"Feature":feature_name, "AUC":[AUC]})
        plot_df = pd.concat([plot_df, temp_df], ignore_index=True)
    
    print(plot_df)
    input()
    fig = px.bar(plot_df, x="Feature", y="AUC", color="Feature", width=650, color_discrete_map={"AlphaFold confidence": "limegreen", "Solvent accessibility":"orange", "Loop length":"lightblue"})
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=False, gridwidth=1, gridcolor='lightgray', showticklabels=False)
    fig.update_yaxes(range=[0, 1], dtick=0.1, showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray', ticks="outside")
    fig.update_layout(title=dict(text="<b>Feature performance</b>", y=0.95, x=0.465, xanchor="center", yanchor="top"),
                xaxis_title="",
                #xaxis_position=0,
                yaxis_title="AUC score",
                paper_bgcolor="white",
                plot_bgcolor="white",
                #bargroupgap=0.1,
                )
    fig.update_traces(width=[0.45]*3)
    
    #fig.write_html(os.path.join(output_path, f"Scikit_AFPDB_EvaluatingFeatures.html"))
    #fig.write_image(os.path.join(output_path, f"Scikit_AFPDB_EvaluatingFeatures.jpg"), width = 650, height = 550, scale = 4)
    
    ''' OneToOne '''
    y = AF_df["is_cut"].values
    positive_indices = np.where(y == 1)[0]
    negative_indices = np.where(y == 0)[0]
    sample_indices = form_sample(positive_indices, negative_indices, "OneToOne")
    AF_subset = AF_df.iloc[sample_indices]
    
    plot_df = pd.DataFrame()
    for feature in features:
        feature_name = features[feature]
        
        y_true = AF_subset["is_cut"].tolist()
        y_score = AF_subset[feature].tolist()
        
        AUC = round(roc_auc_score(y_true, y_score), 3)
        if AUC < 0.5:
            AUC = 1 - AUC
            
        temp_df = pd.DataFrame({"Feature":feature_name, "AUC":[AUC]})
        plot_df = pd.concat([plot_df, temp_df], ignore_index=True)
    
    fig = px.bar(plot_df, x="Feature", y="AUC", color="Feature", width=650, color_discrete_map={"AlphaFold confidence": "limegreen", "Solvent accessibility":"orange", "Loop length":"lightblue"})
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', showgrid=False, gridwidth=1, gridcolor='lightgray', showticklabels=False)
    fig.update_yaxes(range=[0, 1], dtick=0.1, showline=True, linewidth=2, linecolor='black', showgrid=True, gridwidth=1, gridcolor='lightgray', ticks="outside")
    fig.update_layout(title=dict(text="<b>Feature performance (1:1)</b>", y=0.95, x=0.465, xanchor="center", yanchor="top"),
                xaxis_title="",
                yaxis_title="AUC score",
                paper_bgcolor="white",
                plot_bgcolor="white",
                )
    fig.update_traces(width=[0.45]*3)
    #fig.write_html(os.path.join(output_path, f"Scikit_AFPDB_EvaluatingFeaturesOneToOne.html"))
    #fig.write_image(os.path.join(output_path, f"Scikit_AFPDB_EvaluatingFeaturesOneToOne.jpg"), width = 650, height = 550, scale = 4)
    
''' Launch script '''                
if __name__ == "__main__":
    start = time.ctime()
    main() 
    finish = time.ctime()
    print(f"\nScript: {__file__}\nStart: {start}\nFinish: {finish}")